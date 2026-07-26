# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Q → the τ-cocycle lattice Λ_τ(Q) with its Gram matrix (paper §4, step 1)."""

from __future__ import annotations

import time

import numpy as np

from hawk_recover import ring
from hawk_recover.params import HawkParams
from hawk_recover.ring import integer_kernel
from hawk_recover.types import CocycleLattice, PublicKey, require


def _log(verbose: bool, t0: float, msg: str) -> None:
    if verbose:
        print(f"[cocycle {time.time() - t0:6.1f}s] {msg}", flush=True)


def _Q_to_array(pk: PublicKey) -> np.ndarray:
    n = pk.params.n
    Qc = np.zeros((2, 2, n), dtype=np.int64)
    for r in range(2):
        for c in range(2):
            v = pk.Q_list[2 * r + c]
            Qc[r, c, : len(v)] = v
    return Qc


def _build_cBr(Qc: np.ndarray, dL: int, dK: int) -> np.ndarray:
    r"""The Q-constraint matrix on the explicit (C1)-basis P.

    Column j encodes the 4dL coefficients of Q·P_j − σ(P_j)^T·τ(Q)."""
    tauQc = ring.tau(Qc)
    cBr = np.zeros((4 * dL, 2 * dL), dtype=np.int64)
    sgn = np.array([1 - 2 * (k & 1) for k in range(dL)], dtype=np.int64)
    # P_k = diag(ζ^k, (−1)^k ζ^k), k = 0..dL-1.
    for k in range(dL):
        sk = sgn[k]
        cBr[0 * dL : 1 * dL, k] = ring.nshift(Qc[0, 0], k) - sk * ring.nshift(
            tauQc[0, 0], -k
        )
        cBr[1 * dL : 2 * dL, k] = sk * ring.nshift(Qc[0, 1], k) - sk * ring.nshift(
            tauQc[0, 1], -k
        )
        cBr[2 * dL : 3 * dL, k] = ring.nshift(Qc[1, 0], k) - ring.nshift(
            tauQc[1, 0], -k
        )
        cBr[3 * dL : 4 * dL, k] = sk * ring.nshift(Qc[1, 1], k) - ring.nshift(
            tauQc[1, 1], -k
        )
    # P_{dL+j} = [[0, ζ^m],[0,0]], m = 2j+1.
    for j in range(dK):
        m, col = 2 * j + 1, dL + j
        cBr[1 * dL : 2 * dL, col] = ring.nshift(Qc[0, 0], m)
        cBr[2 * dL : 3 * dL, col] = ring.nshift(tauQc[0, 0], -m)
        cBr[3 * dL : 4 * dL, col] = ring.nshift(Qc[1, 0], m) + ring.nshift(
            tauQc[0, 1], -m
        )
    # P_{dL+dK+j} = [[0,0],[ζ^m,0]], m = 2j+1.
    for j in range(dK):
        m, col = 2 * j + 1, dL + dK + j
        cBr[0 * dL : 1 * dL, col] = ring.nshift(Qc[0, 1], m) + ring.nshift(
            tauQc[1, 0], -m
        )
        cBr[1 * dL : 2 * dL, col] = ring.nshift(tauQc[1, 1], -m)
        cBr[2 * dL : 3 * dL, col] = ring.nshift(Qc[1, 1], m)
    return cBr


def _saturated_kernel(cBr: np.ndarray, dL: int, *, verbose: bool, t0: float):
    """Saturated Z-kernel of cBr (rank dL): float64 QR picks dL independent rows,
    then integer_kernel reads the exact kernel off the HNF transform. Its LLL
    pass shrinks entries ~6000 → ~10 bits so the downstream Gram stays small.
    """
    from sage.all import ZZ, matrix
    from scipy.linalg import qr as scipy_qr

    piv = scipy_qr(cBr.T.astype(np.float64), pivoting=True, mode="economic")[-1]
    rows = sorted(int(i) for i in piv[:dL])
    M = cBr[rows, :]
    rk = int(np.linalg.matrix_rank(M.astype(np.float64)))
    require(rk == dL, f"selected-row rank {rk} ≠ dL={dL}")

    Kr = integer_kernel(matrix(ZZ, M.tolist()), verbose=verbose)
    require(Kr.nrows() == dL, f"kernel rank {Kr.nrows()} ≠ dL={dL}")
    return Kr


def _embed_LB(Kr, dL: int, dK: int):
    """Rewrite Kr (rows in the 2dL-dim P-basis) as ambient M₂(O_L)-coords
    (rows of length 4dL). LB[·, (2r+c)·dL + k] = coefficient of ζ^k in W[r,c]."""
    from sage.all import ZZ, matrix

    KrN = np.asarray(Kr, dtype=object)
    sgn = np.array([1 - 2 * (k & 1) for k in range(dL)], dtype=object)
    LBn = np.zeros((dL, 4 * dL), dtype=object)
    LBn[:, 0 * dL : 1 * dL] = KrN[:, 0:dL]
    LBn[:, 3 * dL : 4 * dL] = KrN[:, 0:dL] * sgn
    LBn[:, 1 * dL + 1 : 2 * dL : 2] = KrN[:, dL : dL + dK]
    LBn[:, 2 * dL + 1 : 3 * dL : 2] = KrN[:, dL + dK : 2 * dL]
    return matrix(ZZ, LBn.tolist())


def _gram_tau(LB, dL: int):
    r"""Gram_τ[i,j] = Tr_{L/Q}(det-bilinear(W_i,W_j)).  Since Tr_{L/Q}(Σ a_k ζ^k)
    = dL·a_0, and (f·g)[0] = f J g, this is dL·(A J Dᵀ + D J Aᵀ − B J Cᵀ − C J Bᵀ)."""
    from sage.all import ZZ, matrix

    J = matrix(ZZ, ring.tr_mul_matrix(dL).tolist())
    A = LB[:, 0 * dL : 1 * dL]
    B = LB[:, 1 * dL : 2 * dL]
    C = LB[:, 2 * dL : 3 * dL]
    D = LB[:, 3 * dL : 4 * dL]
    AD = A * J * D.transpose()
    BC = B * J * C.transpose()
    return dL * (AD + AD.transpose() - BC - BC.transpose())


def cocycle_lattice(pk: PublicKey, *, verbose: bool = False) -> CocycleLattice:
    p = pk.params
    dL, dK = p.dL, p.dK
    t0 = time.time()

    Qc = _Q_to_array(pk)
    cBr = _build_cBr(Qc, dL, dK)
    _log(verbose, t0, f"cBr built, max|·|={int(np.abs(cBr).max())}")

    Kr = _saturated_kernel(cBr, dL, verbose=verbose, t0=t0)
    # Kernel check over Z (cBr has small entries; Kr may be large, so use Sage ints).
    from sage.all import ZZ, matrix

    cBrZ = matrix(ZZ, cBr.tolist())
    require((cBrZ * Kr.transpose()).is_zero(), "cBr · Krᵀ ≠ 0")

    LB = _embed_LB(Kr, dL, dK)
    Gt = _gram_tau(LB, dL)
    _log(
        verbose, t0, f"Gram_τ {Gt.nrows()}×{Gt.ncols()}, height {Gt.height().nbits()}b"
    )

    lat = CocycleLattice(params=p, LB=LB, Gram_tau=Gt)
    lat.check()
    # det check against the paper's Prop 4.5 value (basis-independent).
    from sage.all import RR

    dl2 = float(RR(Gt.det().abs()).log(2))
    require(
        abs(dl2 - p.logdet_gram_tau) < 0.5,
        f"log₂ det(Gram_τ) = {dl2:.2f}, expected {p.logdet_gram_tau:.2f}",
    )
    return lat


def cocycle_lattice_small(
    Q_list: tuple[list[int], list[int], list[int], list[int]],
    params: HawkParams,
    *,
    verbose: bool = False,
) -> CocycleLattice:
    """Same as cocycle_lattice but takes bare Q_list (for sub-level Q' in the
    tower, where no encoded pk exists)."""
    pk = PublicKey(params=params, Q_list=Q_list, raw=b"")
    return cocycle_lattice(pk, verbose=verbose)
