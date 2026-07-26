# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""SVP-stage unit tests: NTL serialization round-trip and ambient dedup.

The synthetic ℓ=5–6 pipeline (test_e2e_small) covers the BKZ-only path;
rank ≤ 128 never reaches the sieve, so the sieve backend and harvest_log
are exercised only by the CLI on a real HAWK-256 key (and the mocked
orchestration in test_dispatch).
"""

from __future__ import annotations

from pathlib import Path

from hawk_recover.infra.sieve_backends.bgj import write_ntl
from hawk_recover.svp import dedup_ambient


def test_write_ntl_roundtrip(tmp_path: Path) -> None:
    """NTL output parses back to the same integer matrix (spot-check corners)."""
    from sage.all import ZZ, matrix

    n = 8
    B = matrix(ZZ, n, n, lambda i, j: (37 * i + 101 * j - 250) * (2**35 + i - j))
    ntl = tmp_path / "basis.ntl"
    write_ntl(B, ntl)
    rows = [
        r
        for r in ntl.read_text().replace("[", " ").replace("]", " ").splitlines()
        if r.strip()
    ]
    assert len(rows) == n
    assert [int(x) for x in rows[0].split()] == [int(B[0, j]) for j in range(n)]
    assert [int(x) for x in rows[-1].split()] == [int(B[n - 1, j]) for j in range(n)]


def test_dedup_sign_invariance() -> None:
    """±w and exact duplicates collapse to one ambient representative."""
    from sage.all import ZZ, identity_matrix, vector

    LB = identity_matrix(ZZ, 4)
    v = vector(ZZ, [1, -2, 0, 3])
    w = vector(ZZ, [0, 1, 1, -1])
    uniq = dedup_ambient([v, -v, w, -w, v], LB)
    assert len(uniq) == 2
