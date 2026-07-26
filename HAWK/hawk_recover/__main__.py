# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""CLI: HAWK public key in, verified secret key out (paper §2.1).

Stages (decode → cocycle → svp → intertwiner → tower → encode) checkpoint
into --workdir; re-running resumes from the last completed stage.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="hawk_recover",
        description="HAWK key recovery via the τ-cocycle lattice attack.",
    )
    ap.add_argument("pk", help="public key: hex string, or path to a .hex/.bin file")
    ap.add_argument(
        "-o",
        "--out",
        default="recovered_sk.bin",
        help="output path for the recovered private key, in the reference "
        "implementation's decoded format (default: %(default)s)",
    )
    ap.add_argument(
        "-w",
        "--workdir",
        default="hawk_recover_work",
        help="checkpoint/scratch directory (default: %(default)s)",
    )
    ap.add_argument(
        "--ell",
        type=int,
        default=None,
        help="cyclotomic conductor exponent (default: inferred from the pk "
        "length; end-to-end recovery supports HAWK-256, ell=9, only)",
    )
    ap.add_argument(
        "--svp",
        choices=("auto", "bgj"),
        default="auto",
        help="sieve backend: bgj (needs an AVX512-VNNI CPU) or auto-select (default)",
    )
    ap.add_argument(
        "--skip-verify",
        action="store_true",
        help="skip the final reference-implementation sign/verify check",
    )
    args = ap.parse_args(argv)
    verbose = True
    t0 = time.time()

    def log(msg: str) -> None:
        if verbose:
            print(f"[hawk_recover {time.time() - t0:8.1f}s] {msg}", flush=True)

    p = Path(args.pk)
    if p.exists():
        data = p.read_bytes()
        try:
            pk_bytes = bytes.fromhex(data.decode("ascii").strip())
        except (UnicodeDecodeError, ValueError):
            pk_bytes = data
    else:
        pk_bytes = bytes.fromhex(args.pk.strip())

    # Resolve to absolute: stage code hands paths under this directory to
    # subprocesses that run with a different cwd.
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    from hawk_recover.cocycle import cocycle_lattice
    from hawk_recover.decode import decode
    from hawk_recover.encode import encode
    from hawk_recover.infra import shim
    from hawk_recover.infra.checkpoint import CheckpointDir
    from hawk_recover.infra.sieve_backends import select_backend
    from hawk_recover.intertwiner import intertwiner
    from hawk_recover.svp import find_shell_vectors
    from hawk_recover.tower import recover_basis
    from hawk_recover.types import CocycleLattice, Intertwiner, ShellVectors

    ckpt = CheckpointDir(workdir / "checkpoints")

    # 1. decode
    if args.ell:
        from hawk_recover.params import HawkParams

        pk = decode(pk_bytes, params=HawkParams.for_ell(args.ell))
    else:
        pk = decode(pk_bytes)
    log(f"decoded HAWK-{pk.params.n} public key ({len(pk_bytes)} bytes)")
    if pk.params.ell != 9:
        # The SVP schedule constants and the sign/verify shim are HAWK-256-only.
        print(
            f"error: end-to-end recovery supports HAWK-256 only "
            f"(input decodes as HAWK-{pk.params.n})",
            file=sys.stderr,
        )
        return 2

    # 2. cocycle lattice
    cached = ckpt.load("cocycle.sobj")
    if cached is not None:
        lat = CocycleLattice(
            params=pk.params, LB=cached["LB"], Gram_tau=cached["Gram_tau"]
        )
        lat.check()
        log("cocycle lattice: loaded from checkpoint (re-verified)")
    else:
        lat = cocycle_lattice(pk, verbose=verbose)
        ckpt.save("cocycle.sobj", {"LB": lat.LB, "Gram_tau": lat.Gram_tau})
        log("cocycle lattice: built")

    # 3. shortest vectors
    cached = ckpt.load("shell.sobj")
    if cached is not None:
        shell = ShellVectors(
            params=pk.params,
            W=list(cached["W"]),
            Gram_tau=lat.Gram_tau,
            LB=lat.LB,
            provenance=str(cached.get("provenance", "checkpoint")),
        )
        shell.check()
        log(f"shell vectors: loaded {len(shell.W)} from checkpoint (re-verified)")
    else:
        import os

        backend = select_backend(args.svp, threads=os.cpu_count() or 4)
        log(f"sieve backend: {backend.name}")
        shell = find_shell_vectors(
            lat, workdir=workdir / "sieve", backend=backend, verbose=verbose, ckpt=ckpt
        )
        ckpt.save("shell.sobj", {"W": shell.W, "provenance": shell.provenance})
        log(f"shell vectors: found {len(shell.W)}")

    # 4. intertwiner
    cached = ckpt.load("xtau.sobj")
    if cached is not None:
        itw = Intertwiner(
            params=pk.params,
            wt=cached["wt"],
            V=cached["V"],
            la_rank=int(cached["la_rank"]),
        )
        # A corrupted xtau.sobj must fail here, not as a wrong key three stages
        # later.
        from hawk_recover.ring import field_and_homs

        itw.check(field_and_homs(pk.params.ell)[0])
        log("intertwiner: loaded from checkpoint (re-verified)")
    else:
        itw = intertwiner(shell, pk, lat, verbose=verbose)
        ckpt.save("xtau.sobj", {"wt": itw.wt, "V": itw.V, "la_rank": itw.la_rank})
        log("intertwiner: a-side cocycle V_τ found")

    # 5. tower descent
    rb = recover_basis(itw, pk, verbose=verbose)
    log(f"tower descent: B_rec recovered (β-trace {rb.beta_trace})")

    # 6. encode + verify
    sk = encode(rb, pk)
    out = Path(args.out)
    if args.skip_verify:
        out.write_bytes(sk.raw)
        log(
            f"secret key written to {out} ({len(sk.raw)} bytes) — UNVERIFIED (--skip-verify)"
        )
        return 0

    # The output file only ever contains a key that passed the reference
    # implementation's sign/verify: write to a sidecar, verify, then rename.
    tmp = out.with_suffix(out.suffix + ".unverified")
    tmp.write_bytes(sk.raw)
    try:
        ok = shim.signverify(pk_bytes, sk.raw, b"hawk_recover e2e")
    except (RuntimeError, AssertionError, OSError) as e:
        # The shim itself could not run (not built, not executable, bad
        # inputs) — a setup defect, not a verdict on the key.
        print(f"error: {e} — recovered key left at {tmp}", file=sys.stderr)
        return 3
    log(f"reference-implementation sign/verify: {'PASS' if ok else 'FAIL'}")
    if ok:
        tmp.rename(out)
        log(f"secret key written to {out} ({len(sk.raw)} bytes)")
        return 0
    log(f"key FAILED verification — left at {tmp} for inspection")
    return 1


if __name__ == "__main__":
    sys.exit(main())
