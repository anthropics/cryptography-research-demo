# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Shortest vectors of Λ_τ (paper step 2, heuristic version of §6.2): BKZ
warm-up, then sieve — in place of the provable block reduction of §5."""

from __future__ import annotations

import re
import time
from pathlib import Path

from hawk_recover.types import CocycleLattice, ShellVectors, WarmBasis, require


def _log(verbose: bool, t0: float, msg: str) -> None:
    if verbose:
        print(f"[svp {time.time() - t0:7.1f}s] {msg}", flush=True)


# Phase 1: BKZ warm-up (fpylll, Gram mode)


def bkz_warmup(
    lat: CocycleLattice,
    *,
    schedule: tuple[int, ...] | None = None,
    launch_beta: int | None = None,
    max_loops: int = 6,
    verbose: bool = False,
    ckpt=None,
) -> tuple[WarmBasis, list]:
    """Progressive Gram-mode BKZ on Gram_τ; returns (warm_basis, direct_shell_vectors).

    If ≥2 direct shell vectors are found the caller can skip the sieve.
    `ckpt` saves/resumes the accumulated transform after each β step.
    """
    from fpylll import BKZ, FPLLL, GSO, LLL, IntegerMatrix
    from fpylll.algorithms.bkz2 import BKZReduction
    from sage.all import ZZ, matrix, vector

    p = lat.params
    n = p.rank_tau
    Gt = lat.Gram_tau
    schedule = schedule if schedule is not None else p.beta_warm_schedule
    launch_beta = launch_beta if launch_beta is not None else p.beta_warm_launch
    t0 = time.time()

    U_prev = None
    done_beta = -1
    saved = ckpt.load("warmup.sobj") if ckpt is not None else None
    if saved is not None:
        U_prev = saved["Uim"]
        done_beta = int(saved["beta"])
        if done_beta >= launch_beta:
            # Running further BKZ here would pre-empt find_shell_vectors's
            # own deepening decision.
            from sage.all import ZZ, vector

            _log(
                verbose,
                t0,
                f"resume: checkpoint β={done_beta} ≥ launch β={launch_beta}, no BKZ needed",
            )
            hits = []
            for ri in range(n):
                u = vector(ZZ, U_prev.row(ri))
                if u * Gt * u == p.lam1sq:
                    hits.append(u)
            return _make_warm_basis(lat, U_prev, done_beta), hits
        G_start = U_prev * Gt * U_prev.transpose()
        _log(verbose, t0, f"resuming warm-up from checkpoint at β={done_beta}")
    else:
        G_start = Gt

    A = IntegerMatrix.from_matrix(
        [[int(G_start[i, j]) for j in range(n)] for i in range(n)]
    )
    Ui = IntegerMatrix.identity(n)
    # Float-type ladder: 'd' is 5-10x faster than mpfr and exact enough for
    # LLL-reduced Gram entries; escalate precision on fpylll's reduction/
    # precision errors only (a real bug must propagate, not become "ladder
    # exhausted").
    from fpylll.util import ReductionError

    Mg = None
    last_err: BaseException | None = None
    for ft in ("d", "ld", "mpfr"):
        try:
            if ft == "mpfr":
                FPLLL.set_precision(max(100, Gt.height().nbits() + 50))
            Mg = GSO.Mat(A, U=Ui, gram=True, float_type=ft, flags=GSO.INT_GRAM)
            Mg.update_gso()
            LLL.Reduction(Mg)()
            Mg.update_gso()
            _log(verbose, t0, f"gram-LLL({ft}) r0={Mg.get_r(0, 0):.3e}")
            break
        except (ReductionError, RuntimeError, ValueError) as e:
            _log(verbose, t0, f"float_type={ft} failed: {e}; escalating precision")
            last_err = e
            Mg = None
    if Mg is None:
        raise RuntimeError("gram-LLL failed at every fpylll float type") from last_err

    bkz = BKZReduction(Mg)

    def U_sage():
        """Total transform: (this session's U) · (any resumed prior transform)."""
        U_now = matrix(ZZ, n, n, [Ui[i][j] for i in range(n) for j in range(n)])
        return U_now * U_prev if U_prev is not None else U_now

    def direct_hits(Um) -> list:
        out = []
        for ri in range(n):
            u = vector(ZZ, Um.row(ri))
            if u * Gt * u == p.lam1sq:
                out.append(u)
        return out

    warm = None
    for b in schedule:
        if b <= done_beta:
            continue  # already completed before the resume point
        bkz(
            BKZ.Param(
                b,
                strategies=BKZ.DEFAULT_STRATEGY,
                max_loops=max_loops,
                flags=BKZ.GH_BND | BKZ.AUTO_ABORT,
            )
        )
        Mg.update_gso()
        Um = U_sage()
        if ckpt is not None:
            ckpt.save("warmup.sobj", {"Uim": Um, "beta": b})
        hits = direct_hits(Um)
        if b >= launch_beta or hits:
            warm = _make_warm_basis(lat, Um, b)
        if len(hits) >= 2:
            _log(verbose, t0, f"*** {len(hits)} direct shell vectors at β={b} ***")
            return warm, hits
        if b >= launch_beta:
            # Deepening past the launch point is find_shell_vectors's decision,
            # made only if the sieve comes back short.
            break
    if warm is None and done_beta >= launch_beta:
        # Resumed past the launch point with no further schedule steps to run.
        warm = _make_warm_basis(lat, U_sage(), done_beta)
    require(
        warm is not None, f"schedule {schedule} never reached launch_beta={launch_beta}"
    )
    return warm, direct_hits(U_sage())


def _make_warm_basis(lat: CocycleLattice, Um, beta: int) -> WarmBasis:
    """Scaled integer Cholesky of the reduced Gram, for the sieve backend."""
    from sage.all import ZZ, Integer, RealField, ceil, log, matrix

    p = lat.params
    n = p.rank_tau
    Gred = Um * lat.Gram_tau * Um.transpose()
    RF = RealField(200)
    Lch = matrix(RF, Gred).cholesky()
    mx = max(abs(Lch[i, j]) for i in range(n) for j in range(i + 1))
    # Scale so entries use ~40 bits: enough precision that rounding error
    # cannot shift a λ₁²-vector past the harvest threshold, small enough that
    # svp_tool's doubles don't overflow.
    sb = max(14, 40 - int(ceil(log(mx, 2))))
    Bint = matrix(
        ZZ,
        n,
        n,
        [
            [
                Integer((Lch[i, j] * RF(2) ** sb).round()) if j <= i else 0
                for j in range(n)
            ]
            for i in range(n)
        ],
    )
    return WarmBasis(
        params=p, beta=beta, Uim=Um, Bint=Bint, sb=sb, Gram_tau=lat.Gram_tau, LB=lat.LB
    )


def solve_svp_gram(
    Ga,
    target: int,
    *,
    schedule: range | None = None,
    max_loops: int = 8,
    extra_rounds: int = 3,
    verbose: bool = False,
    tag: str = "svp",
    return_beta: bool = False,
):
    """Progressive Gram-mode BKZ on Ga; returns transform rows of exact squared
    norm `target` (or (rows, beta_first_hit) when return_beta=True). Heuristic
    only — the paper's provable route is Ducas's block reduction (§5); the
    target's gap in these lattices makes moderate β suffice in practice
    (§6.2). `extra_rounds` keeps going after the first hit to collect more
    candidates.
    """
    from fpylll import BKZ, FPLLL, GSO, LLL, IntegerMatrix
    from fpylll.algorithms.bkz2 import BKZReduction
    from sage.all import ZZ, matrix, vector

    t0 = time.time()
    d = Ga.nrows()

    def hits_of(Um) -> list:
        out = []
        for ri in range(d):
            u = vector(ZZ, Um.row(ri))
            if u * Ga * u == target:
                out.append(u)
        return out

    # Tiny dimensions: exact enumeration is instant and complete.
    if d <= 16:
        from sage.all import QuadraticForm

        QF = QuadraticForm(ZZ, 2 * Ga)
        sv = QF.short_vector_list_up_to_length(target + 1)
        out = [vector(ZZ, v) for v in sv[target]] if target < len(sv) else []
        return (out, 0) if return_beta else out

    A = IntegerMatrix.from_matrix([[int(Ga[i, j]) for j in range(d)] for i in range(d)])
    Ui = IntegerMatrix.identity(d)
    FPLLL.set_precision(max(100, Ga.height().nbits() + 50))
    Mg = GSO.Mat(A, U=Ui, gram=True, float_type="mpfr", flags=GSO.INT_GRAM)
    Mg.update_gso()
    LLL.Reduction(Mg)()
    Mg.update_gso()
    _log(verbose, t0, f"[{tag}] dim {d} LLL r0={Mg.get_r(0, 0):.4g} target={target}")

    def U_sage():
        return matrix(ZZ, d, d, [Ui[i][j] for i in range(d) for j in range(d)])

    bkz = BKZReduction(Mg)
    schedule = schedule if schedule is not None else range(4, 91, 4)
    rounds_after_hit = 0
    beta_first_hit = 0  # 0 = already at target after plain LLL
    for b in schedule:
        if int(round(Mg.get_r(0, 0))) <= target:
            rounds_after_hit += 1
            if rounds_after_hit > extra_rounds:
                break
        else:
            beta_first_hit = b  # overwritten until a step reaches the target
        bkz(
            BKZ.Param(
                b,
                strategies=BKZ.DEFAULT_STRATEGY,
                max_loops=max_loops,
                flags=BKZ.GH_BND | BKZ.AUTO_ABORT,
            )
        )
        Mg.update_gso()

    hits = hits_of(U_sage())
    return (hits, beta_first_hit) if return_beta else hits


# Phase 2: harvest — map sieve output back to Λ_τ coordinates

# svp_tool's output line format: "length = <float>(<float>), vec = [x0 x1 ...]"
_VEC_RE = re.compile(r"length = ([\d.]+)\([^)]*\), vec = \[([^\]]+)\]")


def harvest_log(log_text: str, warm: WarmBasis, *, slack: float = 3.0) -> list:
    """Extract integer shell vectors (exact Gram_τ-norm λ₁²) from svp_tool output."""
    from sage.all import ZZ, RealField, matrix, vector

    p = warm.params
    n = p.rank_tau
    RF = RealField(100)
    BiR = matrix(RF, warm.Bint)
    tgt_scaled = RF(warm.target_sq_scaled)
    out = []
    for lstr, vstr in _VEC_RE.findall(log_text):
        # Cheap float filter before the exact check: scaled length² ≤ slack·target.
        if RF(lstr) ** 2 / tgt_scaled > slack:
            continue
        row = [RF(x) for x in vstr.split()]
        if len(row) != n:
            continue
        t = [RF(0)] * n
        for j in range(n - 1, -1, -1):
            s = row[j] - sum(t[i] * BiR[i, j] for i in range(j + 1, n))
            t[j] = s / BiR[j, j]
        ti = vector(ZZ, [ZZ(x.round()) for x in t])
        w = ti * warm.Uim
        if w * warm.Gram_tau * w == p.lam1sq:
            out.append(w)
    return out


def dedup_ambient(ws: list, LB) -> list:
    """Dedup by sign-canonicalized ambient coordinates V = w·LB (basis-independent)."""
    seen: dict[tuple, object] = {}
    for w in ws:
        V = w * LB
        nz = next((i for i, x in enumerate(V) if x != 0), 0)
        key = tuple(V) if V[nz] > 0 else tuple(-V)
        seen.setdefault(key, w)
    return list(seen.values())


def find_shell_vectors(
    lat: CocycleLattice,
    *,
    workdir: str | Path,
    backend=None,
    verbose: bool = False,
    ckpt=None,
) -> ShellVectors:
    """CocycleLattice → ≥2 independent shortest vectors (ShellVectors).

    Dim ≤ 128: progressive BKZ alone reaches λ₁ in practice (paper App. A);
    above that (HAWK-256, dim 256), BKZ warm-up is followed by the sieve
    backend (§6.2). Heuristic throughout: no provable claim is inherited from
    the paper's block reduction (§5).
    """
    p = lat.params
    t0 = time.time()

    if p.rank_tau <= 128:
        hits = dedup_ambient(
            solve_svp_gram(lat.Gram_tau, p.lam1sq, verbose=verbose, tag="Λ_τ"), lat.LB
        )
        _log(verbose, t0, f"BKZ-only: {len(hits)} unique shell vectors")
        sv = ShellVectors(
            params=p, W=hits, Gram_tau=lat.Gram_tau, LB=lat.LB, provenance="bkz"
        )
        sv.check()
        return sv

    if backend is None:
        from hawk_recover.infra.sieve_backends import select_backend

        backend = select_backend("auto")
    _log(
        verbose,
        t0,
        f"sieve backend: {getattr(backend, 'name', type(backend).__name__)}",
    )

    # Adaptive schedule (paper §6.2): sieve from a light warm-up; if <2 shell
    # vectors, deepen by one β step and re-sieve, up to beta_warm_max. Each
    # attempt sieves in a directory named by its warm-up β: the sieve basis is
    # warm-up-specific, and a resume that lands at a deeper β (crash after a
    # deepening step) must not reuse the shallower attempt's basis file.
    workdir = Path(workdir)
    launch = p.beta_warm_launch
    while True:
        warm, hits = bkz_warmup(lat, launch_beta=launch, verbose=verbose, ckpt=ckpt)
        # A resumed checkpoint can be deeper than requested (crash mid-deepening);
        # sync the loop index so deepening always means *new* BKZ work.
        launch = max(launch, warm.beta)
        hits = dedup_ambient(hits, lat.LB)
        _log(verbose, t0, f"warm-up at β={warm.beta}: {len(hits)} direct hits")
        if len(hits) >= 2:
            provenance = f"bkz{warm.beta}"
            break

        attempt_dir = workdir.parent / f"{workdir.name}_b{warm.beta}"
        hits = dedup_ambient(
            hits + backend.sieve_to_target(warm, workdir=attempt_dir, verbose=verbose),
            lat.LB,
        )
        _log(verbose, t0, f"sieve from β={warm.beta}: {len(hits)} unique shell vectors")
        if len(hits) >= 2:
            provenance = f"bkz{warm.beta}+{backend.name}-sieve"
            break

        if launch >= p.beta_warm_max:
            raise RuntimeError(
                f"<2 shell vectors even at the deepest warm-up (β={p.beta_warm_max}); "
                f"this exceeds every observed HAWK-{p.n} instance — check the input key."
            )
        launch += 4
        _log(
            verbose,
            t0,
            f"insufficient — deepening warm-up to β={launch} and re-sieving",
        )

    sv = ShellVectors(
        params=p, W=hits, Gram_tau=lat.Gram_tau, LB=lat.LB, provenance=provenance
    )
    sv.check()
    return sv
