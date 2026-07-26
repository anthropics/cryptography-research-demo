#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Time the four attack stages on a synthetic key at a small conductor exponent.

Reproduces the README's small-ℓ timing table: cocycle → svp → intertwiner →
tower on a freshly generated random key (the paper's three steps of §2.1, with
step 3 split into intertwiner + tower as in §6.2). decode and
encode need the real HAWK byte encoding, defined only for n ≥ 256, so they are
skipped; the correctness check here is B_rec*·B_rec = Q.

Usage: nix-shell nix/shell.nix --run 'sage -python scripts/synthetic_attack.py [ell] [seed]'
       (default ell=7, seed=1)
"""

from __future__ import annotations

import random
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from hawk_recover import ring  # noqa: E402
from hawk_recover.cocycle import cocycle_lattice  # noqa: E402
from hawk_recover.intertwiner import intertwiner  # noqa: E402
from hawk_recover.params import HawkParams  # noqa: E402
from hawk_recover.svp import find_shell_vectors  # noqa: E402
from hawk_recover.tower import recover_basis  # noqa: E402
from hawk_recover.types import PublicKey  # noqa: E402


def synthetic_key(ell: int, seed: int):
    """A random B ∈ GL₂(O_L) (product of elementary shears) and its Gram Q."""
    p = HawkParams.for_ell(ell)
    n = p.n
    L, z, c, _, _ = ring.field_and_homs(ell)
    rng = random.Random(seed)

    def small_elt():
        return sum(rng.randint(-1, 1) * z**k for k in range(n))

    B = ring.gidentity(L, 2)
    for _ in range(3):
        B = (
            B
            * ring.gmat(L, [[1, small_elt()], [0, 1]])
            * ring.gmat(L, [[1, 0], [small_elt(), 1]])
        )
    Q = ring.star(B, c) * B
    Q_list = tuple(
        [int(co) for co in list(Q[r, cc])] + [0] * (n - len(list(Q[r, cc])))
        for r in range(2)
        for cc in range(2)
    )
    return PublicKey(params=p, Q_list=Q_list, raw=b""), L, c, Q


def main(argv: list[str]) -> int:
    ell = int(argv[1]) if len(argv) > 1 else 7
    seed = int(argv[2]) if len(argv) > 2 else 1
    print(f"synthetic HAWK instance: ell={ell}, n={1 << (ell - 1)}, seed={seed}")

    pk, L, c, Q = synthetic_key(ell, seed)
    T = time.time()

    def stage(name: str, t0: float) -> float:
        print(f"  {name:<12} {time.time() - t0:8.1f}s")
        return time.time()

    t = time.time()
    lat = cocycle_lattice(pk)
    t = stage("cocycle", t)
    with tempfile.TemporaryDirectory() as d:
        shell = find_shell_vectors(lat, workdir=d)
    t = stage("svp", t)
    itw = intertwiner(shell, pk, lat)
    t = stage("intertwiner", t)
    rb = recover_basis(itw, pk)
    t = stage("tower", t)

    rb.check()
    Brec = rb.B(L)
    ok = ring.star(Brec, c) * Brec == Q
    print(f"B_rec*·B_rec = Q: {ok}")
    print(f"total: {time.time() - T:.1f}s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
