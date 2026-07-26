# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

r"""BGJ-Sieve-AMX backend (``svp_tool``): Zhao–Ding–Yang, *Sieving with Streaming Memory Access*.

Sliding-window passes (params.sieve_msd_*), then one full ``-f`` pass **without**
``-amx`` (the bgj3-amx kernel stalls on the final window — a known svp_tool bug).
"""

from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

from hawk_recover.types import WarmBasis

_REF_IMPL = Path(__file__).parent.parent.parent / "ref_impl"
_LOCAL_BUILD = _REF_IMPL / "bgj_sieve_src" / "BGJ-Sieve-AMX-main" / "app" / "svp_tool"


def _log(verbose: bool, t0: float, msg: str) -> None:
    if verbose:
        print(f"[bgj {time.time() - t0:7.1f}s] {msg}", flush=True)


def _cpu_flags() -> set[str]:
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return set()
    for line in text.splitlines():
        if line.startswith("flags"):
            return set(line.split(":", 1)[1].split())
    return set()


def write_ntl(Bint, path: str | Path) -> None:
    """svp_tool input: NTL mat_ZZ format."""
    n = Bint.nrows()
    with open(path, "w") as f:
        f.write("[")
        for i in range(n):
            f.write("[" + " ".join(str(Bint[i, j]) for j in range(n)) + "]\n")
        f.write("]\n")


# svp_tool has a compiled-in MAX_NTHREADS (include/bgj_epi8.h); passing more
# threads overflows its fixed-size per-thread arrays and the stack protector
# kills the process ("stack smashing detected", no output). The Nix build
# raises the upstream limit of 112 to 192 (postPatch in nix/sieve-src.nix);
# this constant must stay equal to that build-time value.
_BGJ_MAX_THREADS = 192


@dataclass
class BGJBackend:
    """Sieve via the svp_tool binary (subprocess, .ntl files, restart-safe)."""

    name = "bgj"

    threads: int = 50
    ssd: int = 48  # sieve start dimension (research-run default)
    # Pass schedule: None → use the parameter set's (params.sieve_msd_*).
    msd_slide: tuple[int, ...] | None = None
    msd_final: int | None = None
    esd_final: int | None = None

    @classmethod
    def availability(cls) -> tuple[bool, str]:
        flags = _cpu_flags()
        local_ok = (
            _LOCAL_BUILD.exists() and (_LOCAL_BUILD.parent / "run_local.sh").exists()
        )
        if not local_ok:
            return (
                False,
                "no svp_tool binary (run scripts/setup.sh to build it with Nix)",
            )
        if "avx512_vnni" not in flags:
            return False, "CPU lacks required ISA (svp_tool needs avx512_vnni)"
        return True, ""

    def _launcher(self) -> Path:
        """The Nix-built binary's wrapper. availability() has already
        validated that it exists and that the CPU can run it."""
        return (_LOCAL_BUILD.parent / "run_local.sh").resolve()

    def _effective_threads(self) -> int:
        """self.threads clamped to the binary's compiled-in thread limit."""
        return min(self.threads, _BGJ_MAX_THREADS)

    def _run_pass(self, ntl: Path, args: list[str], log_fh, log_path: Path) -> None:
        """One svp_tool pass; raises RuntimeError if no output basis appears.

        All paths handed to the subprocess must be absolute: it runs with a
        different cwd, and svp_tool's error for an unopenable file is the
        misleading 'incorrect file format'."""
        ntl = ntl.resolve()
        threads = self._effective_threads()
        cmd = [str(self._launcher()), str(ntl), *args, "-t", str(threads)]
        if "-v" not in args:
            cmd += ["-v", "1"]
        env = {"OMP_NUM_THREADS": str(threads), "PATH": "/usr/bin:/bin"}
        t0 = time.time()
        subprocess.run(
            cmd,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
            cwd=str(ntl.parent),
        )
        log_fh.flush()
        out = ntl.with_name(ntl.name + "r")
        if not out.exists():
            tail = ""
            try:
                tail = log_path.read_text()[-500:]
            except OSError:
                pass
            raise RuntimeError(
                f"svp_tool pass produced no output after {time.time() - t0:.1f}s "
                f"(cmd: {' '.join(cmd[:4])}…). Log tail:\n{tail}"
            )
        out.replace(ntl)

    def sieve_to_target(
        self, warm: WarmBasis, *, workdir: str | Path, verbose: bool = False
    ) -> list:
        """Run the pass schedule, then harvest exact λ₁² vectors from the
        final pass's output log."""
        from hawk_recover.svp import harvest_log

        t0 = time.time()
        p = warm.params
        msd_slide = self.msd_slide if self.msd_slide is not None else p.sieve_msd_slide
        msd_final = self.msd_final if self.msd_final is not None else p.sieve_msd_final
        esd_final = self.esd_final if self.esd_final is not None else p.sieve_esd_final
        workdir = Path(workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)
        n = p.rank_tau
        ntl = workdir / "basis.ntl"
        progress_file = workdir / "passes_done.txt"
        if not ntl.exists():
            write_ntl(warm.Bint, ntl)
            progress_file.write_text("0")
        done = int(progress_file.read_text() or "0") if progress_file.exists() else 0

        # -amx accelerates the sliding passes on AMX machines; never on the
        # final -f pass (bgj3-amx stall bug).
        amx_args = ["-amx"] if "amx_tile" in _cpu_flags() else []
        slide_log = workdir / "slide.log"
        pass_idx = 0
        with open(slide_log, "a") as fh:
            for msd in msd_slide:
                for lo in range(0, n - msd - 19, 20):
                    pass_idx += 1
                    if pass_idx <= done:
                        continue  # completed before a restart
                    hi = min(lo + msd + 20, n)
                    self._run_pass(
                        ntl,
                        [
                            "--msd",
                            str(msd),
                            "--ind_l",
                            str(lo),
                            "--ind_r",
                            str(hi),
                            "--ssd",
                            str(self.ssd),
                            *amx_args,
                        ],
                        fh,
                        slide_log,
                    )
                    # Progress recorded only for passes that verifiably
                    # produced output (_run_pass raises otherwise).
                    progress_file.write_text(str(pass_idx))
                _log(verbose, t0, f"slide msd={msd} done ({pass_idx} passes)")

        final_log = workdir / "final.log"
        with open(final_log, "w") as fh:
            self._run_pass(
                ntl,
                ["--msd", str(msd_final), "--esd", str(esd_final), "-f", "-v", "3"],
                fh,
                final_log,
            )
        _log(verbose, t0, f"final -f msd={msd_final} done")

        return harvest_log(final_log.read_text(), warm)
