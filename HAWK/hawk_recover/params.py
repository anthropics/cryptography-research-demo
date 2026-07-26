# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Parameter sets. Every magic constant is derived from ell or cited from the
HAWK specification / the attack paper.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class HawkParams:
    """All dimensions follow paper §3.1 notation: L = Q(ζ_{2^ell}), n = [L:Q]."""

    ell: int

    # HAWK encoding constants (ref impl hawk_vrfy.c / hawk_kgen.c)
    # Segregated Golomb–Rice: low-bit width and |coeff| bit-limit for q00, q01.
    k_gr_q00: int
    lim_bits_q00: int
    k_gr_q01: int
    lim_bits_q01: int
    # Encoded pk byte length (hawk.h:HAWK_PUBKEY_SIZE).
    pk_bytes: int
    # hpub = SHAKE256(pk)[:hpub_bytes], appended to sk. hpub_bytes = 1<<(logn-4).
    hpub_bytes: int

    # Derived field dimensions (paper §3.1)
    @property
    def m(self) -> int:
        """Cyclotomic conductor 2^ell."""
        return 1 << self.ell

    @property
    def n(self) -> int:
        """Ring degree [L:Q] = φ(2^ell) = 2^{ell-1}. HAWK-n's n."""
        return 1 << (self.ell - 1)

    @property
    def dL(self) -> int:
        return self.n

    @property
    def dK(self) -> int:
        """[K:Q] where K = L^c is the maximal totally real subfield (§3.1)."""
        return self.n >> 1

    # τ-cocycle lattice constants (paper §4)
    @property
    def rank_tau(self) -> int:
        """rank Λ_τ = n (Lemma 4.3)."""
        return self.n

    @property
    def lam1sq(self) -> int:
        """λ₁(Λ_τ)² under Tr_L-form = 2·dL (Lemma 4.4, V_τ has det 1)."""
        return 2 * self.dL

    @property
    def logdet_gram_tau(self) -> float:
        """log₂ det(Gram_τ) under the Tr_L∘det-bilinear form. By Prop 4.5 (in
        the ×8 Gram convention of types.CocycleLattice) Λ_τ ≅ √(2dL)·(Z^{dK+1}
        ⊕ √2·Z^{dK−1}), so det = (2dL)^{dK+1}·(4dL)^{dK−1};
        e.g. 2431 at ell=9. Used only as a basis-independent self-check."""
        from math import log2

        return (self.dK + 1) * log2(2 * self.dL) + (self.dK - 1) * log2(4 * self.dL)

    # SVP stage schedule (paper §6.2; research-run values).
    # β=44 launch + the extended slide schedule below measured ~2× faster
    # end-to-end than β=56 + (80,95,108); automatic warm-up deepening means the
    # fast launch point can never produce a worse outcome, only an earlier one.
    beta_warm_schedule: tuple[int, ...] = field(
        default_factory=lambda: tuple(range(20, 60, 4))
    )
    beta_warm_launch: int = 44
    beta_warm_max: int = 56
    # The two low ease-in levels (60, 70) do the reduction work BKZ β=48..56
    # would have done, at ~4× better wall-clock on multithreaded sieve kernels.
    sieve_msd_slide: tuple[int, ...] = (60, 70, 80, 95, 108)
    sieve_msd_final: int = 118
    sieve_esd_final: int = 8

    # tower max-β per level (research-run values; cf. paper §6.2)
    tower_max_beta: dict[int, int] = field(
        default_factory=lambda: {9: 120, 8: 60, 7: 40, 6: 20, 5: 10, 4: 5}
    )

    @property
    def sk_core_bytes(self) -> int:
        """|f|+|g| as int8 + |F mod 2|+|G mod 2| as bit-packed (spec §3.4)."""
        return 2 * self.n + 2 * (self.n // 8)

    @classmethod
    def for_ell(cls, ell: int) -> HawkParams:
        # Per hawk_vrfy.c decode_q00/decode_q01: (low, lim_bits).
        match ell:
            case 9:
                return cls(
                    ell=9,
                    k_gr_q00=5,
                    lim_bits_q00=9,
                    k_gr_q01=8,
                    lim_bits_q01=11,
                    pk_bytes=450,
                    hpub_bytes=16,
                )
            case 10:
                return cls(
                    ell=10,
                    k_gr_q00=5,
                    lim_bits_q00=9,
                    k_gr_q01=9,
                    lim_bits_q01=12,
                    pk_bytes=1024,
                    hpub_bytes=32,
                )
            case 11:
                return cls(
                    ell=11,
                    k_gr_q00=6,
                    lim_bits_q00=10,
                    k_gr_q01=10,
                    lim_bits_q01=14,
                    pk_bytes=2440,
                    hpub_bytes=64,
                )
            case _:
                # Synthetic small-ℓ test instances; encoding widths only matter
                # for real encode/decode, which small-ℓ tests skip.
                return cls(
                    ell=ell,
                    k_gr_q00=5,
                    lim_bits_q00=9,
                    k_gr_q01=8,
                    lim_bits_q01=11,
                    pk_bytes=0,
                    hpub_bytes=16,
                )

    @classmethod
    def from_pk_len(cls, nbytes: int) -> HawkParams:
        for ell in (9, 10, 11):
            p = cls.for_ell(ell)
            if p.pk_bytes == nbytes:
                return p
        raise ValueError(f"no HAWK parameter set has pk_bytes={nbytes}")


HAWK256 = HawkParams.for_ell(9)
HAWK512 = HawkParams.for_ell(10)
HAWK1024 = HawkParams.for_ell(11)
