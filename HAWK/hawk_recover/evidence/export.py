# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Trusted-side export: a run's .sobj checkpoints → the inert evidence format.

CLI: sage -python -m hawk_recover.evidence export <workdir> <out_dir> --pk <pk>
"""

from __future__ import annotations

import json
from pathlib import Path

from hawk_recover.evidence import EVIDENCE_VERSION, _mat_to_rows, _vec_to_row


def export_evidence(
    workdir: str | Path,
    out_dir: str | Path,
    *,
    pk_path: str | Path | None = None,
    sieve_log: str | Path | None = None,
    warm_path: str | Path | None = None,
) -> Path:
    """Convert a run's .sobj checkpoints (trusted side) into inert JSON/text
    that a third party can verify without executing code."""
    from sage.all import load

    workdir = Path(workdir)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ck = workdir / "checkpoints"

    coc = load(str(ck / "cocycle.sobj"))
    (out / "lattice.json").write_text(
        json.dumps(
            {
                "version": EVIDENCE_VERSION,
                "LB": _mat_to_rows(coc["LB"]),
                "Gram_tau": _mat_to_rows(coc["Gram_tau"]),
            }
        )
    )

    shl = load(str(ck / "shell.sobj"))
    (out / "shell.json").write_text(
        json.dumps(
            {
                "version": EVIDENCE_VERSION,
                "W": [_vec_to_row(w) for w in shl["W"]],
                "provenance": str(shl.get("provenance", "")),
            }
        )
    )

    wsrc = Path(warm_path) if warm_path else ck / "warmup.sobj"
    if wsrc.exists():
        wrm = load(str(wsrc))
        wdata = {
            "version": EVIDENCE_VERSION,
            "beta": int(wrm["beta"]),
            "Uim": _mat_to_rows(wrm["Uim"]),
        }
        if "Bint" in wrm:
            wdata["Bint"] = _mat_to_rows(wrm["Bint"])
            wdata["sb"] = int(wrm["sb"])
        (out / "warm.json").write_text(json.dumps(wdata))

    # The sieve log is optional evidence: only an explicitly named log is
    # exported (a HAWK-256 run's is at <workdir>/sieve_b<beta>/final.log).
    if sieve_log is not None:
        (out / "sieve.log").write_text(Path(sieve_log).read_text())

    if pk_path is not None:
        raw = Path(pk_path).read_bytes()
        try:
            hex_text = bytes.fromhex(raw.decode("ascii").strip()).hex()
        except (UnicodeDecodeError, ValueError):
            hex_text = raw.hex()
        (out / "pk.hex").write_text(hex_text)

    return out


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(
        prog="hawk_recover.evidence export",
        description="Export a run's checkpoints to the inert JSON evidence format.",
    )
    ap.add_argument("workdir", help="the run's --workdir")
    ap.add_argument("out_dir", help="evidence output directory")
    ap.add_argument("--pk", default=None, help="the run's input pk (.hex/.bin)")
    ap.add_argument(
        "--sieve-log",
        default=None,
        help="sieve log to include (a HAWK-256 run's is at "
        "<workdir>/sieve_b<beta>/final.log; omit for BKZ-only runs)",
    )
    ap.add_argument(
        "--warm",
        default=None,
        help="warm basis .sobj (default: <workdir>/checkpoints/warmup.sobj)",
    )
    args = ap.parse_args(argv)
    if args.sieve_log is not None and not Path(args.sieve_log).is_file():
        ap.error(f"--sieve-log {args.sieve_log} is not a file")
    out = export_evidence(
        args.workdir,
        args.out_dir,
        pk_path=args.pk,
        sieve_log=args.sieve_log,
        warm_path=args.warm,
    )
    print(f"evidence exported to {out}")
    for f in sorted(out.iterdir()):
        print(f"  {f.name}  ({f.stat().st_size} bytes)")
    return 0
