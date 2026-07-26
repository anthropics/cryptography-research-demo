# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""SVP dispatch tests for the HAWK-256-scale path (rank > 128).

Every other test uses small conductors, which take the BKZ-only early return
in find_shell_vectors — so nothing else executes the adaptive warm-up/sieve
loop. These tests exercise that loop with mocked BKZ and a mock backend, so
bugs in the orchestration code (imports, path handling, deepening logic) are
caught without running a real multi-hour sieve.
"""

from __future__ import annotations

import pytest

from hawk_recover.params import HAWK256
from hawk_recover.types import CocycleLattice, WarmBasis


@pytest.fixture()
def synthetic_256_lattice():
    """A rank-256 'lattice' whose minimal vectors are trivially known.

    Gram_tau = lam1sq·I means every standard basis vector has exact norm
    lam1sq — the orchestration code can't tell it from a real Λ_τ, and no
    BKZ/sieve math is needed to find its shortest vectors.
    """
    from sage.all import ZZ, identity_matrix, matrix

    n = HAWK256.rank_tau
    Gt = HAWK256.lam1sq * identity_matrix(ZZ, n)
    LB = matrix(
        ZZ, n, 4 * HAWK256.dL, {(i, i): 1 for i in range(n)}
    )  # any full-rank embedding works
    return CocycleLattice(params=HAWK256, LB=LB, Gram_tau=Gt)


class MockBackend:
    """Returns canned vectors; records how it was called."""

    name = "mock"

    def __init__(self, results_per_call: list[list]):
        self.results_per_call = list(results_per_call)
        self.calls: list[dict] = []

    @classmethod
    def availability(cls):
        return True, ""

    def sieve_to_target(self, warm: WarmBasis, *, workdir, verbose=False) -> list:
        self.calls.append({"beta": warm.beta, "workdir": str(workdir)})
        return self.results_per_call.pop(0) if self.results_per_call else []


def _fake_warmup_factory(hits_per_beta=None):
    """A bkz_warmup stand-in: returns a structurally valid WarmBasis at the
    requested launch β without running BKZ."""
    from sage.all import ZZ, identity_matrix

    def fake_warmup(lat, *, launch_beta=None, verbose=False, ckpt=None, **kw):
        n = lat.params.rank_tau
        beta = launch_beta if launch_beta is not None else lat.params.beta_warm_launch
        warm = WarmBasis(
            params=lat.params,
            beta=beta,
            Uim=identity_matrix(ZZ, n),
            Bint=identity_matrix(ZZ, n),
            sb=20,
            Gram_tau=lat.Gram_tau,
            LB=lat.LB,
        )
        hits = (hits_per_beta or {}).get(beta, [])
        return warm, hits

    return fake_warmup


def _unit_vectors(lat, count):
    from sage.all import ZZ, vector

    n = lat.params.rank_tau
    return [vector(ZZ, [1 if j == i else 0 for j in range(n)]) for i in range(count)]


def test_dispatch_sieve_path_runs(synthetic_256_lattice, monkeypatch, tmp_path):
    """rank>128 dispatches to the warm-up + backend path and returns the
    backend's vectors. Catches import/path/orchestration errors in the loop."""
    import hawk_recover.svp as svp

    lat = synthetic_256_lattice
    monkeypatch.setattr(svp, "bkz_warmup", _fake_warmup_factory())

    backend = MockBackend([_unit_vectors(lat, 3)])
    sv = svp.find_shell_vectors(
        lat, workdir=tmp_path / "sieve", backend=backend, verbose=False
    )
    assert len(sv.W) == 3
    assert backend.calls[0]["beta"] == HAWK256.beta_warm_launch
    assert "mock" in sv.provenance


def test_dispatch_deepening_on_insufficient_vectors(
    synthetic_256_lattice, monkeypatch, tmp_path
):
    """<2 vectors from the first sieve attempt deepens the warm-up by one
    step and re-sieves; the second attempt's vectors are returned."""
    import hawk_recover.svp as svp

    lat = synthetic_256_lattice
    monkeypatch.setattr(svp, "bkz_warmup", _fake_warmup_factory())

    backend = MockBackend(
        [
            _unit_vectors(lat, 1),  # attempt 1: insufficient
            _unit_vectors(lat, 2),
        ]
    )  # attempt 2: enough
    sv = svp.find_shell_vectors(
        lat, workdir=tmp_path / "sieve", backend=backend, verbose=False
    )
    assert len(sv.W) == 2
    assert len(backend.calls) == 2
    assert backend.calls[1]["beta"] == HAWK256.beta_warm_launch + 4
    # Attempts use separate workdirs (the basis is warm-up-specific).
    assert backend.calls[0]["workdir"] != backend.calls[1]["workdir"]


def test_dispatch_exhaustion_raises(synthetic_256_lattice, monkeypatch, tmp_path):
    """Reaching beta_warm_max with <2 vectors raises rather than looping."""
    import hawk_recover.svp as svp

    lat = synthetic_256_lattice
    monkeypatch.setattr(svp, "bkz_warmup", _fake_warmup_factory())

    n_attempts = (HAWK256.beta_warm_max - HAWK256.beta_warm_launch) // 4 + 1
    backend = MockBackend([[] for _ in range(n_attempts)])
    with pytest.raises(RuntimeError, match="deepest warm-up"):
        svp.find_shell_vectors(
            lat, workdir=tmp_path / "sieve", backend=backend, verbose=False
        )


def test_dispatch_direct_hits_skip_sieve(synthetic_256_lattice, monkeypatch, tmp_path):
    """≥2 direct hits from the warm-up skip the backend entirely."""
    import hawk_recover.svp as svp

    lat = synthetic_256_lattice
    direct = {HAWK256.beta_warm_launch: _unit_vectors(lat, 2)}
    monkeypatch.setattr(svp, "bkz_warmup", _fake_warmup_factory(direct))

    backend = MockBackend([])
    sv = svp.find_shell_vectors(
        lat, workdir=tmp_path / "sieve", backend=backend, verbose=False
    )
    assert len(sv.W) == 2
    assert backend.calls == []
    assert "mock" not in sv.provenance


def test_bkz_warmup_stops_at_launch_beta(monkeypatch):
    """The warm-up must stop at launch_beta, not run its whole schedule.

    Without the break, the warm-up silently runs to the end of the schedule,
    re-creating the conservative recipe and making find_shell_vectors's
    deepening logic unreachable (regression found in review). lam1sq is patched to an
    unreachable value so direct hits cannot provide an alternative early exit:
    the ONLY way this test passes is the launch_beta break itself.
    """
    import random
    from unittest.mock import PropertyMock

    from hawk_recover import ring
    from hawk_recover.cocycle import cocycle_lattice
    from hawk_recover.params import HawkParams
    from hawk_recover.svp import bkz_warmup
    from hawk_recover.types import PublicKey

    # Small synthetic lattice (ℓ=6, dim 32) so real BKZ runs in seconds.
    ell = 6
    p = HawkParams.for_ell(ell)
    L, z, c, _, _ = ring.field_and_homs(ell)
    rng = random.Random(7)

    def ro():
        return sum(rng.randint(-1, 1) * z**k for k in range(p.n))

    B = ring.gidentity(L, 2)
    for _ in range(3):
        B = B * ring.gmat(L, [[1, ro()], [0, 1]]) * ring.gmat(L, [[1, 0], [ro(), 1]])
    Q = ring.star(B, c) * B
    Q_list = tuple(
        [int(co) for co in list(Q[r, cc])] + [0] * (p.n - len(list(Q[r, cc])))
        for r in range(2)
        for cc in range(2)
    )
    pk = PublicKey(params=p, Q_list=Q_list, raw=b"")
    lat = cocycle_lattice(pk)

    # No vector can have this norm -> direct_hits is always empty -> the
    # launch_beta break is the only early exit.
    monkeypatch.setattr(
        HawkParams, "lam1sq", PropertyMock(return_value=10**9), raising=True
    )
    # PropertyMock on the class: accessing p.lam1sq calls it.
    monkeypatch.setattr(type(lat.params), "lam1sq", property(lambda self: 10**9))

    warm, hits = bkz_warmup(lat, schedule=(4, 8, 12, 16, 20, 24, 28), launch_beta=12)
    assert hits == [], "direct hits should be impossible with patched lam1sq"
    assert warm.beta == 12, (
        f"warm-up ran to β={warm.beta}; the launch_beta=12 break is missing or broken"
    )


def test_bgj_backend_clamps_threads_to_binary_limit():
    """svp_tool has a compiled-in MAX_NTHREADS; more threads crash it with a
    stack overflow and no output. The Nix build raises the upstream limit of
    112 to 192 (postPatch in nix/sieve-src.nix); the clamp must stay equal
    to that build-time value.
    """
    from hawk_recover.infra.sieve_backends.bgj import (
        _BGJ_MAX_THREADS,
        BGJBackend,
    )

    assert _BGJ_MAX_THREADS == 192
    # Requests under the limit pass through unchanged...
    assert BGJBackend(threads=4)._effective_threads() == 4
    assert BGJBackend(threads=100)._effective_threads() == 100
    # ...and requests beyond it clamp to the compiled-in limit.
    assert BGJBackend(threads=256)._effective_threads() == _BGJ_MAX_THREADS


def test_bkz_warmup_resume_at_launch_runs_no_bkz(synthetic_256_lattice, tmp_path):
    """Resuming from a checkpoint at (or past) launch_beta must rebuild the
    warm basis from the checkpoint and run ZERO further BKZ steps — the
    re-validation run crashed because resume ran an unnecessary extra step
    (whose numerical failure modes the warm-up's float ladder doesn't cover).
    """
    from sage.all import ZZ, identity_matrix

    from hawk_recover.infra.checkpoint import CheckpointDir
    from hawk_recover.params import HAWK256
    from hawk_recover.svp import bkz_warmup

    lat = synthetic_256_lattice
    n = HAWK256.rank_tau
    ck = CheckpointDir(tmp_path / "checkpoints")
    # A checkpoint exactly at the launch point (the re-validation's state).
    ck.save(
        "warmup.sobj", {"Uim": identity_matrix(ZZ, n), "beta": HAWK256.beta_warm_launch}
    )

    warm, hits = bkz_warmup(lat, ckpt=ck)
    # If any BKZ step had run, this would take minutes on a dim-256 lattice;
    # the identity transform surviving proves no reduction happened.
    assert warm.beta == HAWK256.beta_warm_launch
    assert warm.Uim == identity_matrix(ZZ, n)
    # Gram = lam1sq*I means every basis row is a direct hit.
    assert len(hits) == n
