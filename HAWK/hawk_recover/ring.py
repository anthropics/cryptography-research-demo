# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Coefficient-vector arithmetic in O_L = Z[ζ_m] ≅ Z[x]/(x^n+1) via numpy int64.

Elements are length-n int64 arrays in the power basis 1,ζ,…,ζ^{n-1}. The
Galois involutions τ:ζ↦-ζ, c:ζ↦ζ⁻¹, σ=cτ:ζ↦-ζ⁻¹ act coefficientwise.
"""

from __future__ import annotations

import numpy as np

from hawk_recover.types import require


def nshift(f: np.ndarray, m: int) -> np.ndarray:
    """ζ^m · f in Z[x]/(x^n+1), f shape (..., n)."""
    n = f.shape[-1]
    m = m % (2 * n)
    s = 1
    if m >= n:
        m -= n
        s = -1
    if m == 0:
        return s * f
    out = np.empty_like(f)
    out[..., :m] = -f[..., n - m :]
    out[..., m:] = f[..., : n - m]
    return s * out


def tau(f: np.ndarray) -> np.ndarray:
    """τ: ζ ↦ -ζ, i.e. f_k ↦ (-1)^k f_k."""
    n = f.shape[-1]
    sgn = np.array([1 - 2 * (k & 1) for k in range(n)], dtype=f.dtype)
    return f * sgn


def tr_mul_matrix(n: int) -> np.ndarray:
    """The n×n matrix J such that (f·g)[0] = f @ J @ g, i.e. the constant
    coefficient of the negacyclic product. J[0,0]=1, J[k,n-k]=-1 (k≥1)."""
    J = np.zeros((n, n), dtype=np.int64)
    J[0, 0] = 1
    for k in range(1, n):
        J[k, n - k] = -1
    return J


# exact Python-int negacyclic arithmetic (no overflow concerns)


def nc_mul(a: list[int], b: list[int], n: int) -> list[int]:
    """Product in Z[x]/(x^n + 1), schoolbook, arbitrary-precision ints."""
    c = [0] * n
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            k = i + j
            if k < n:
                c[k] += ai * bj
            else:
                c[k - n] -= ai * bj
    return c


def nc_conj(a: list[int], n: int) -> list[int]:
    """c: ζ ↦ ζ⁻¹ on coefficient lists."""
    return [a[0]] + [-a[n - k] for k in range(1, n)]


def nc_tau(a: list[int], n: int) -> list[int]:
    """τ: ζ ↦ −ζ on coefficient lists."""
    return [(-x if (k & 1) else x) for k, x in enumerate(a)]


def mat2_mul(
    A: list[list[list[int]]], B: list[list[list[int]]], n: int
) -> list[list[list[int]]]:
    """Product of 2×2 matrices whose entries are negacyclic coefficient lists."""
    out = [[[0] * n for _ in range(2)] for _ in range(2)]
    for r in range(2):
        for c in range(2):
            acc = [0] * n
            for k in range(2):
                p = nc_mul(A[r][k], B[k][c], n)
                for i in range(n):
                    acc[i] += p[i]
            out[r][c] = acc
    return out


def mat2_apply(f, A: list[list[list[int]]], n: int) -> list[list[list[int]]]:
    """Apply a coefficient-list map f (e.g. nc_tau, nc_conj) entrywise."""
    return [[f(A[r][c], n) for c in range(2)] for r in range(2)]


def mat2_adj(A: list[list[list[int]]], n: int) -> list[list[list[int]]]:
    """Adjugate of a 2×2 matrix: [[d,−b],[−c,a]]."""
    neg = lambda v: [-x for x in v]
    return [[A[1][1], neg(A[0][1])], [neg(A[1][0]), A[0][0]]]


def mat2_transpose_conj(A: list[list[list[int]]], n: int) -> list[list[list[int]]]:
    """A* = conjugate-transpose (c entrywise, then transpose)."""
    Ac = mat2_apply(nc_conj, A, n)
    return [[Ac[0][0], Ac[1][0]], [Ac[0][1], Ac[1][1]]]


def negacirc(f: list[int], n: int):
    """Sage ZZ matrix of multiplication by f in Z[x]/(x^n+1):
    column k is the coefficient vector of ζ^k·f."""
    from sage.all import ZZ, matrix

    cols: list[list[int]] = []
    cur = list(f)
    for _ in range(n):
        cols.append(list(cur))
        # multiply by ζ: shift up one degree, wrap with sign flip.
        cur = [-cur[-1]] + cur[:-1]
    return matrix(ZZ, cols).transpose()


def field_and_homs(ell: int):
    """Sage CyclotomicField L with the three involutions. For verification
    and the slow-path / small-ℓ stages only; hot loops use the numpy ops."""
    from sage.all import CyclotomicField

    m = 1 << ell
    L = CyclotomicField(m, f"z{ell}")
    z = L.gen()
    return L, z, L.hom([z**-1]), L.hom([-z]), L.hom([-(z**-1)])


# safe matrix algebra over cyclotomic fields
#
# Sage's specialized Matrix_cyclo_dense class computes matrix products with a
# multimodular algorithm whose height bound is WRONG for moderately large
# entries (observed: conductor 64, ~12-bit coefficients give silently incorrect
# products). Every matrix over a cyclotomic field in this package must
# therefore use the generic (naive, exact) implementation:
#   * construct with gmat / gidentity,
#   * apply Galois homs with mat_apply (Sage's .apply_map() reconstructs the
#     result in the *default* MatrixSpace, silently re-entering the buggy class),
#   * convert externally produced matrices with to_generic before use.
# Matrices over ZZ/QQ are unaffected (different classes).


def mat_space(L, nrows: int = 2, ncols: int | None = None):
    """Generic-implementation MatrixSpace over L (Sage caches these)."""
    from sage.all import MatrixSpace

    return MatrixSpace(
        L, nrows, ncols if ncols is not None else nrows, implementation="generic"
    )


def gmat(L, entries, nrows: int = 2, ncols: int | None = None):
    """Construct a matrix over L in the generic implementation."""
    return mat_space(L, nrows, ncols)(entries)


def gidentity(L, n: int = 2):
    return mat_space(L, n).one()


def mat_apply(M, hom):
    """Entrywise ring-hom application, staying in the generic implementation."""
    MS = mat_space(M.base_ring(), M.nrows(), M.ncols())
    return MS([hom(e) for e in M.list()])


def to_generic(M):
    """Re-create any Sage matrix over a cyclotomic field in the generic class."""
    return gmat(M.base_ring(), M.list(), M.nrows(), M.ncols())


def star(M, c):
    """Conjugate-transpose M* = c(M)ᵀ in the generic implementation."""
    return mat_apply(M, c).transpose()


def integer_kernel(M, *, lll: bool = True, verbose: bool = False):
    """Saturated Z-basis of the right kernel {w : M·w = 0}, as rows, via
    HNF-with-transform of Mᵀ. The optional LLL pass keeps entry sizes small."""
    H, U = M.transpose().hermite_form(transformation=True)
    nz = sum(1 for i in range(H.nrows()) if not H.row(i).is_zero())
    K = U[nz:, :]
    require((M * K.transpose()).is_zero(), "kernel check failed")
    if lll and K.nrows() > 0:
        K = K.LLL()
    return K
