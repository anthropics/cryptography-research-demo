# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Shell vectors → the a-side shortest cocycle V_τ via the intertwiner
sublattice Λ_a (paper §6.2, Remark 6.2; the first half of step 3)."""

from __future__ import annotations

import time
from typing import Any

from hawk_recover import ring
from hawk_recover.ring import integer_kernel
from hawk_recover.svp import solve_svp_gram
from hawk_recover.types import (
    CocycleLattice,
    Intertwiner,
    PublicKey,
    SageMatrix,
    SageVector,
    ShellVectors,
    require,
)


def _log(verbose: bool, t0: float, msg: str) -> None:
    if verbose:
        print(f"[intw {time.time() - t0:6.1f}s] {msg}", flush=True)


def _lift_mat2(amb, n: int) -> list[list[list[int]]]:
    """Ambient 4n-coordinate vector → 2×2 matrix of coefficient lists."""
    return [
        [[int(amb[(2 * r + c) * n + k]) for k in range(n)] for c in range(2)]
        for r in range(2)
    ]


def _build_system(N, M_y, LB, n: int, *, verbose: bool, t0: float):
    """Constraint matrix A (4n × rank(Λ_τ)) for N·W − W·M_y = 0 over the Λ_τ basis."""
    from sage.all import block_matrix

    CN = [[ring.negacirc(N[r][k], n) for k in range(2)] for r in range(2)]
    CM = [[ring.negacirc(M_y[k][c], n) for c in range(2)] for k in range(2)]

    # LBᵀ slices: LBt[a][b] = (LB columns for ambient entry (a,b))ᵀ, n × rank.
    LBt = [
        [LB[:, (2 * a + b) * n : (2 * a + b + 1) * n].transpose() for b in range(2)]
        for a in range(2)
    ]

    blocks = []
    for r in range(2):
        for c in range(2):
            blk = (
                CN[r][0] * LBt[0][c]
                + CN[r][1] * LBt[1][c]
                - CM[0][c] * LBt[r][0]
                - CM[1][c] * LBt[r][1]
            )
            blocks.append([blk])
    A = block_matrix(blocks, subdivide=False)
    _log(verbose, t0, f"system A {A.nrows()}×{A.ncols()}, height {A.height().nbits()}b")
    return A


def _aside_cocycle(amb: SageVector, n: int, L: Any) -> SageMatrix | None:
    """V = W(amb) over L if amb is the a-side cocycle — W ≡ I (mod 2) and
    W·τ(W) = I, i.e. W = ±V_τ (paper App. D parity test) — else None."""
    from hawk_recover.tower import _amb_to_matrix, _is_aside

    if not _is_aside(amb, n):
        return None
    V = _amb_to_matrix(amb, L, L.gen(), n)
    tau = L.hom([-L.gen()])
    return V if V * ring.mat_apply(V, tau) == ring.gidentity(L, 2) else None


def intertwiner(
    shell: ShellVectors,
    pk: PublicKey,
    lat: CocycleLattice,
    *,
    verbose: bool = False,
    max_pairs: int = 6,
) -> Intertwiner:
    """Find the a-side cocycle V_τ from ≥2 shell vectors (paper Remark 6.2).

    A shell vector that is already the a-side cocycle *is* ±V_τ (paper's
    shortcut); otherwise Remark 6.2 needs a pair with neither member ±V_τ, so
    cocycle vectors are excluded from pairing.
    """
    from itertools import combinations

    from sage.all import ZZ, vector

    p = pk.params
    t0 = time.time()
    shell.check()

    L, _, _, _, _ = ring.field_and_homs(p.ell)
    rest = []
    for i, w in enumerate(shell.W):
        V = _aside_cocycle(vector(ZZ, w * lat.LB), p.n, L)
        if V is not None:
            _log(verbose, t0, f"shell vector {i} is the a-side cocycle ✓")
            result = Intertwiner(params=p, wt=vector(ZZ, w), V=V, la_rank=p.n // 2)
            result.check(L)
            return result
        rest.append(i)

    pairs = list(combinations(rest, 2))[:max_pairs]
    last: AssertionError | None = None
    for i, j in pairs:
        try:
            return _intertwiner_pair(shell, pk, lat, i, j, verbose=verbose, t0=t0)
        except AssertionError as e:
            last = e
            _log(verbose, t0, f"pair ({i},{j}): {e} — trying next pair")
    raise AssertionError(
        f"no shell-vector pair yielded the a-side cocycle ({len(pairs)} pairs tried); "
        f"last failure: {last}"
    )


def _intertwiner_pair(
    shell: ShellVectors,
    pk: PublicKey,
    lat: CocycleLattice,
    i: int,
    j: int,
    *,
    verbose: bool,
    t0: float,
) -> Intertwiner:
    p = pk.params
    n = p.n

    # N and M_y (exact Python-int 2×2 polynomial algebra)
    V0 = _lift_mat2(shell.W[i] * shell.LB, n)
    V1 = _lift_mat2(shell.W[j] * shell.LB, n)
    Q = [[list(pk.Q_list[2 * r + c]) for c in range(2)] for r in range(2)]

    tauV1 = ring.mat2_apply(ring.nc_tau, V1, n)
    N = ring.mat2_mul(V0, tauV1, n)

    tauQ = ring.mat2_apply(ring.nc_tau, Q, n)
    adj_tauQ = ring.mat2_adj(tauQ, n)
    V0_star = ring.mat2_transpose_conj(V0, n)
    M_y = ring.mat2_mul(ring.mat2_mul(ring.mat2_mul(adj_tauQ, V0_star, n), Q, n), V1, n)
    _log(verbose, t0, "N, M_y computed")

    # Λ_a = kernel of the intertwining system, inside Λ_τ
    A = _build_system(N, M_y, lat.LB, n, verbose=verbose, t0=t0)
    La = integer_kernel(A, verbose=verbose)
    rank = La.nrows()
    _log(verbose, t0, f"Λ_a rank = {rank} (expect {n // 2})")
    require(rank == n // 2, f"Λ_a rank {rank} ≠ n/2 = {n // 2}")

    # SVP on Λ_a (gap √2 ⇒ moderate β suffices)
    Ga = La * lat.Gram_tau * La.transpose()
    hits = solve_svp_gram(Ga, p.lam1sq, verbose=verbose, tag="Λ_a")
    require(bool(hits), "no vector of norm λ₁² found in Λ_a")

    from sage.all import ZZ, vector

    L, _, _, _, _ = ring.field_and_homs(p.ell)
    for u in hits:
        wt = vector(ZZ, u * La)
        V = _aside_cocycle(wt * lat.LB, n, L)
        if V is None:
            continue
        result = Intertwiner(params=p, wt=wt, V=V, la_rank=rank)
        result.check(L)
        _log(verbose, t0, "a-side cocycle V_τ found ✓")
        return result

    raise AssertionError(
        f"none of {len(hits)} norm-λ₁² vectors in Λ_a is an a-side cocycle"
    )
