# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Inert evidence format for cross-boundary verification: JSON (decimal strings)
+ plain-text logs; loading executes no code. Sage pickles (.sobj) stay INTRA-run
only — pickle executes code on load, so they must never cross a trust boundary.

This module is the format itself (serializers + loaders); the CLI is
`sage -python -m hawk_recover.evidence {export,verify} …`.
"""

from __future__ import annotations

import json
from pathlib import Path

EVIDENCE_VERSION = 1


# serialization (Sage objects -> JSON-safe)


def _mat_to_rows(M) -> list[list[str]]:
    return [[str(int(x)) for x in row] for row in M.rows()]


def _vec_to_row(v) -> list[str]:
    return [str(int(x)) for x in v]


def _rows_to_mat(
    rows: list[list[str]], nrows: int | None = None, ncols: int | None = None
):
    from sage.all import ZZ, matrix

    data = [[int(x) for x in row] for row in rows]
    if nrows is not None and (len(data) != nrows or any(len(r) != ncols for r in data)):
        raise ValueError(
            f"matrix shape {len(data)}×{len(data[0]) if data else 0}, "
            f"expected {nrows}×{ncols}"
        )
    return matrix(ZZ, data)


def _row_to_vec(row: list[str]):
    from sage.all import ZZ, vector

    return vector(ZZ, [int(x) for x in row])


# load (untrusted side: JSON/text only, no code execution)


def load_lattice(evidence_dir: str | Path):
    """Returns (LB, Gram_tau) as Sage integer matrices from inert JSON."""
    d = json.loads((Path(evidence_dir) / "lattice.json").read_text())
    _check_version(d)
    return _rows_to_mat(d["LB"]), _rows_to_mat(d["Gram_tau"])


def load_shell(evidence_dir: str | Path):
    """Returns (W, provenance) from inert JSON."""
    d = json.loads((Path(evidence_dir) / "shell.json").read_text())
    _check_version(d)
    return [_row_to_vec(w) for w in d["W"]], d.get("provenance", "")


def load_warm(evidence_dir: str | Path):
    """Returns a dict with beta/Uim (+ Bint/sb when present) from inert JSON,
    or None when the run has no warm evidence."""
    p = Path(evidence_dir) / "warm.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    _check_version(d)
    out = {"beta": int(d["beta"]), "Uim": _rows_to_mat(d["Uim"])}
    if "Bint" in d:
        out["Bint"] = _rows_to_mat(d["Bint"])
        out["sb"] = int(d["sb"])
    return out


def load_sieve_log(evidence_dir: str | Path) -> str | None:
    p = Path(evidence_dir) / "sieve.log"
    return p.read_text() if p.exists() else None


def load_pk_bytes(evidence_dir: str | Path) -> bytes | None:
    p = Path(evidence_dir) / "pk.hex"
    return bytes.fromhex(p.read_text().strip()) if p.exists() else None


def _check_version(d: dict) -> None:
    v = d.get("version")
    if v != EVIDENCE_VERSION:
        raise ValueError(
            f"evidence format version {v}, this tool reads {EVIDENCE_VERSION}"
        )


