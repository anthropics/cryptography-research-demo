# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Stage boundary dataclasses. All matrices/vectors are Sage objects; lists of
ints are plain Python. Each carries the invariants the next stage relies on.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from hawk_recover.params import HawkParams

# Sage types are not importable for static annotation without pulling sage.all;
# we alias them as Any and rely on runtime invariant checks.
SageMatrix = Any
SageVector = Any
SageElement = Any


class VerificationError(AssertionError):
    """A verification invariant failed.

    Subclasses AssertionError but is raised explicitly, so `python -O` cannot
    strip it; paper-correctness checks must raise this, never bare `assert`.
    """


def require(cond: bool, msg: str) -> None:
    """`-O`-proof assertion for verification-path invariants."""
    if not cond:
        raise VerificationError(msg)


@dataclass
class PublicKey:
    """Decoded HAWK public Gram matrix Q = B*·B over O_L (paper §3.1)."""

    params: HawkParams
    # Coefficient lists (length n each) for q00, q01, q10, q11 in the power
    # basis of L. q10 = c(q01), q11 = (1 + q01·q10)/q00; all integral.
    Q_list: tuple[list[int], list[int], list[int], list[int]]
    raw: bytes

    @property
    def ell(self) -> int:
        return self.params.ell

    def Q(self, L: Any) -> SageMatrix:
        from sage.all import ZZ

        from hawk_recover import ring

        z = L.gen()
        n = self.params.n
        return ring.gmat(
            L,
            [
                [
                    sum(ZZ(self.Q_list[2 * r + c][k]) * z**k for k in range(n))
                    for c in range(2)
                ]
                for r in range(2)
            ],
        )


@dataclass
class CocycleLattice:
    """The rank-n τ-cocycle lattice Λ_τ(Q) (paper §4), in a concrete Z-basis.

    LB: n×4n Z-basis in ambient (2×2)·O_L coords. Gram_tau = 8·Q^{(τ)} of §4.1
    ([L:F]=4 from Tr_L vs Tr_F, ×2 from omitting the paper's ½ in the associated
    bilinear form); λ₁² = 2·dL (Lemma 4.4).
    """

    params: HawkParams
    LB: SageMatrix
    Gram_tau: SageMatrix

    @property
    def lam1sq(self) -> int:
        return self.params.lam1sq

    def check(self) -> None:
        p = self.params
        require(
            self.LB.nrows() == p.rank_tau and self.LB.ncols() == 4 * p.dL,
            f"LB shape {self.LB.nrows()}×{self.LB.ncols()}, expected {p.rank_tau}×{4 * p.dL}",
        )
        require(self.Gram_tau.is_symmetric(), "Gram_tau not symmetric")
        require(
            self.Gram_tau.nrows() == p.rank_tau,
            f"Gram_tau has {self.Gram_tau.nrows()} rows, expected {p.rank_tau}",
        )


@dataclass
class WarmBasis:
    """BKZ-reduced basis of Λ_τ at a given β, plus scaled Cholesky for svp_tool."""

    params: HawkParams
    beta: int
    Uim: SageMatrix  # n×n unimodular, Uim·LB is the reduced basis
    Bint: SageMatrix  # n×n lower-tri int Cholesky of Uim·Gram_tau·Uim^T, scaled 2^sb
    sb: int
    Gram_tau: SageMatrix
    LB: SageMatrix

    @property
    def target_sq_scaled(self) -> int:
        return self.params.lam1sq * (4**self.sb)


@dataclass
class ShellVectors:
    """≥2 independent shortest vectors of Λ_τ (paper Lemma 4.4), in the LB basis."""

    params: HawkParams
    W: list[SageVector]  # each w ∈ Z^n with w·Gram_tau·w = lam1sq
    Gram_tau: SageMatrix
    LB: SageMatrix
    provenance: str = ""

    def check(self) -> None:
        require(len(self.W) >= 2, f"need ≥2 shell vectors, have {len(self.W)}")
        t = self.params.lam1sq
        for i, w in enumerate(self.W):
            nw = w * self.Gram_tau * w
            require(nw == t, f"W[{i}]·Gt·W[{i}] = {nw} ≠ {t}")


@dataclass
class Intertwiner:
    """The a-side shortest cocycle V_τ with V_τ·τ(V_τ)=I (paper Lemma 4.1),
    V_τ ≡ I (mod 2) — the mod-2 congruence is what singles out ±V_τ among the
    minimal vectors (τ ≡ id mod 2, so V_τ = B⁻¹τ(B) ≡ I; see intertwiner.py)."""

    params: HawkParams
    wt: SageVector | None  # Z^n coords in the LB basis (None if V was given directly)
    V: SageMatrix  # 2×2 over L, = lift_LB(wt)
    la_rank: int  # rank of Λ_a used to isolate it (should be n/2)

    def check(self, L: Any) -> None:
        from hawk_recover import ring

        tau = L.hom([-L.gen()])
        V = ring.to_generic(self.V)
        require(
            V * ring.mat_apply(V, tau) == ring.gidentity(L, 2),
            "V·τ(V) ≠ I (not a cocycle)",
        )


@dataclass
class RecoveredBasis:
    """B_rec ∈ GL₂(O_L) with B_rec*·B_rec = Q — a solution to HAWK-KR(n)
    (paper §3.1, Thm 6.1)."""

    params: HawkParams
    # Coefficient lists for the 4 entries of B_rec.
    B_rec_list: tuple[list[int], list[int], list[int], list[int]]
    beta_trace: dict[int, int | str]
    eq_G: bool
    in_OL: bool
    det_unit: bool

    def B(self, L: Any) -> SageMatrix:
        from sage.all import QQ

        from hawk_recover import ring

        z = L.gen()
        return ring.gmat(
            L,
            [
                [
                    sum(
                        QQ(self.B_rec_list[2 * r + c][k]) * z**k
                        for k in range(len(self.B_rec_list[2 * r + c]))
                    )
                    for c in range(2)
                ]
                for r in range(2)
            ],
        )

    def check(self) -> None:
        require(self.eq_G, "B_rec*·B_rec ≠ Q")
        require(self.in_OL, "B_rec ∉ M₂(O_L)")
        require(self.det_unit, "|N(det B_rec)| ≠ 1")


@dataclass
class SecretKey:
    """Recovered private key in the reference implementation's *decoded*
    layout — f,g (int8) ‖ F,G (mod-2 packed) [‖ hpub], as consumed by
    hawk_sign_finish_alt() (592 bytes for HAWK-256). This is not the NIST
    encoded sk, whose (f,g) are regenerated from a seed that no attack can
    reconstruct; the two are functionally equivalent for signing."""

    params: HawkParams
    f: list[int]
    g: list[int]
    F: list[int]
    G: list[int]
    core: bytes  # f‖g‖F2‖G2, len = sk_core_bytes
    hpub: bytes = b""

    @property
    def raw(self) -> bytes:
        return self.core + self.hpub

    def check(self) -> None:
        p = self.params
        for name, v in (("f", self.f), ("g", self.g)):
            require(
                len(v) == p.n and all(-128 <= x < 128 for x in v),
                f"{name} out of int8 range",
            )
        require(
            len(self.core) == p.sk_core_bytes,
            f"core is {len(self.core)} bytes, expected {p.sk_core_bytes}",
        )
