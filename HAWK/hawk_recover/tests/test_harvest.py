# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""harvest_log: the only bridge from svp_tool's text output back to Λ_τ.

A vector w with exact norm λ₁² is written into a synthetic sieve log in the
scaled-Cholesky coordinates svp_tool prints; harvest must recover exactly w
(and drop garbage lines). Also pins the Bint scaling invariant harvest relies on."""

from __future__ import annotations

from hawk_recover.svp import _make_warm_basis, harvest_log
from hawk_recover.types import CocycleLattice


def _identity_lattice():
    """Gram_tau = lam1sq·I: every unit vector has exact norm lam1sq."""
    from sage.all import ZZ, identity_matrix, matrix

    from hawk_recover.params import HAWK256

    n = HAWK256.rank_tau
    Gt = HAWK256.lam1sq * identity_matrix(ZZ, n)
    LB = matrix(ZZ, n, 4 * HAWK256.dL, {(i, i): 1 for i in range(n)})
    return CocycleLattice(params=HAWK256, LB=LB, Gram_tau=Gt)


def test_bint_scaling_invariant() -> None:
    """Bint·Bintᵀ ≈ 4^sb·(Uim·Gram_tau·Uim^T): the rounding must be tighter
    than half a λ₁² step so a target vector can't cross the harvest threshold."""
    from sage.all import ZZ, identity_matrix

    lat = _identity_lattice()
    n = lat.params.rank_tau
    warm = _make_warm_basis(lat, identity_matrix(ZZ, n), beta=44)
    err = warm.Bint * warm.Bint.transpose() - (4**warm.sb) * lat.Gram_tau
    assert max(abs(x) for x in err.list()) * 8 < 4**warm.sb * warm.params.lam1sq


def test_harvest_recovers_vector_and_ignores_garbage() -> None:
    from sage.all import ZZ, RealField, identity_matrix, vector

    lat = _identity_lattice()
    n = lat.params.rank_tau
    warm = _make_warm_basis(lat, identity_matrix(ZZ, n), beta=44)

    # A λ₁²-vector w (a unit vector, since Gram = lam1sq·I), in the coordinates
    # svp_tool prints: row = tᵀ·Bint with t = w (Uim = I).
    w = vector(ZZ, [1 if j == 3 else 0 for j in range(n)])
    RF = RealField(100)
    row = [RF(sum(w[i] * warm.Bint[i, j] for i in range(n))) for j in range(n)]
    length = (sum(x * x for x in row)) ** 0.5
    good = f"length = {float(length):.6f}(0.0), vec = [{' '.join(str(float(x)) for x in row)}]"
    garbage = "length = 999999.0(0.0), vec = [" + " ".join("7.5" for _ in range(n)) + "]"
    shortline = "length = 1.0(0.0), vec = [1.0 2.0]"  # wrong dimension → skipped

    hits = harvest_log("\n".join([garbage, good, shortline]), warm)
    assert len(hits) == 1
    assert list(hits[0]) == list(w)
    assert hits[0] * warm.Gram_tau * hits[0] == warm.params.lam1sq


def test_harvest_empty_log() -> None:
    from sage.all import ZZ, identity_matrix

    lat = _identity_lattice()
    warm = _make_warm_basis(lat, identity_matrix(ZZ, lat.params.rank_tau), beta=44)
    assert harvest_log("no vectors here\n", warm) == []
