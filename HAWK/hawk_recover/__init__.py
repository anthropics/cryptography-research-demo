# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""HAWK key recovery via the τ-cocycle lattice attack (paper §§2–6).

Pipeline: decode → cocycle → svp → intertwiner → tower → encode.
Each stage is a pure function over the dataclasses in `types`; `checkpoint`
provides opt-in persistence for long-running stages.
"""

from hawk_recover.params import HawkParams
from hawk_recover.types import (
    CocycleLattice,
    Intertwiner,
    PublicKey,
    RecoveredBasis,
    SecretKey,
    ShellVectors,
    WarmBasis,
)

__all__ = [
    "HawkParams",
    "PublicKey",
    "CocycleLattice",
    "WarmBasis",
    "ShellVectors",
    "Intertwiner",
    "RecoveredBasis",
    "SecretKey",
]
