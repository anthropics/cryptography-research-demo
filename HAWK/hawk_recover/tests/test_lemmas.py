# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Smoke test for scripts/verify_lemmas.py: the paper's §4 statements
(Lemmas 4.1–4.4, Prop 4.5) hold on a random key at ℓ=4."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))


def test_lemmas_ell4() -> None:
    import verify_lemmas

    verify_lemmas.verify(4)
