# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Pluggable sieve backends for the SVP stage (paper §6.2, step 2).

A backend turns a warm-reduced basis (WarmBasis) into exact-λ₁² vectors in
Λ_τ coordinates. One implementation ships: `bgj` (the BGJ-Sieve-AMX binary,
needs a VNNI CPU); it reports why when it cannot run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from hawk_recover.types import WarmBasis


@runtime_checkable
class SieveBackend(Protocol):
    """What the SVP stage needs from a sieve implementation."""

    #: registry key and display name ("bgj", …)
    name: str

    @classmethod
    def availability(cls) -> tuple[bool, str]:
        """Can this backend run on this machine?

        Returns (True, "") when usable, or (False, reason) naming the missing
        requirement (CPU feature, missing binary/package, …). Must be cheap
        and must not raise."""
        ...

    def sieve_to_target(
        self, warm: WarmBasis, *, workdir: str | Path, verbose: bool = False
    ) -> list:
        """Run the sieve and return all found vectors w (Λ_τ-basis integer
        coordinate vectors) with w·Gram_τ·w = λ₁² exactly.

        Implementations should be restart-friendly: re-invocation with the
        same workdir should resume or at least not redo completed work."""
        ...


def available_backends() -> dict[str, tuple[bool, str]]:
    """Availability of every registered backend on this machine."""
    from hawk_recover.infra.sieve_backends.bgj import BGJBackend

    return {cls.name: cls.availability() for cls in (BGJBackend,)}


def select_backend(name: str = "auto", *, threads: int = 50) -> SieveBackend:
    """Instantiate a sieve backend ("auto" or a registry name); raise
    RuntimeError naming the blocker when it cannot run here."""
    from hawk_recover.infra.sieve_backends.bgj import BGJBackend

    registry: dict[str, type] = {BGJBackend.name: BGJBackend}

    if name == "auto":
        reasons = []
        for cls in registry.values():
            ok, why = cls.availability()
            if ok:
                return cls(threads=threads)
            reasons.append(f"  {cls.name}: {why}")
        raise RuntimeError(
            "no sieve backend can run on this machine:\n"
            + "\n".join(reasons)
            + "\nSee README §Requirements & hardware."
        )

    if name not in registry:
        raise ValueError(
            f"unknown sieve backend {name!r}; choose from {sorted(registry)} or 'auto'"
        )
    cls = registry[name]
    ok, why = cls.availability()
    if not ok:
        raise RuntimeError(f"sieve backend {name!r} is not usable here: {why}")
    return cls(threads=threads)
