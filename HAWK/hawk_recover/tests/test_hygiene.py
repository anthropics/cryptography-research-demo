# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Code-hygiene guards.

These tests read the package source and fail when a forbidden pattern is
reintroduced. They exist because the failure modes they guard against are
silent (wrong results, not errors).
"""

from __future__ import annotations

import re
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent

# Modules whose matrix arithmetic feeds the paper's correctness claims.
VERIFICATION_MODULES = [
    "types.py",
    "encode.py",
    "cocycle.py",
    "intertwiner.py",
    "tower.py",
    "ring.py",
]

# The Sage Matrix_cyclo_dense multimodular bug (see README §Known Sage issue)
# produces silently wrong products. Every cyclotomic matrix must be built
# through the ring.py helpers (gmat / mat_space / gidentity), and entrywise
# homomorphisms applied via ring.mat_apply — never the patterns below.
FORBIDDEN = [
    (
        re.compile(r"(?<!ring\.mat_space)\(L, *\d+, *implementation=(?!\"generic\")"),
        'MatrixSpace over L without implementation="generic"',
    ),
    (
        re.compile(r"\.apply_map\("),
        ".apply_map() silently reverts to the buggy Matrix_cyclo_dense class; "
        "use ring.mat_apply",
    ),
    (
        re.compile(r"(?<!\w)matrix\(L[,)]"),
        "matrix(L, ...) constructs Matrix_cyclo_dense; use ring.gmat",
    ),
    (
        re.compile(r"(?<!\w)matrix\(Ls[,)]"),
        "matrix(Ls, ...) constructs Matrix_cyclo_dense; use ring.gmat",
    ),
]


def test_no_raw_cyclotomic_matrix_construction() -> None:
    violations = []
    for mod in VERIFICATION_MODULES:
        src = (PKG / mod).read_text()
        # ring.py legitimately contains the MatrixSpace(..., implementation=
        # "generic") calls itself; the patterns are written to allow those.
        for lineno, line in enumerate(src.splitlines(), 1):
            if line.strip().startswith("#"):
                continue
            for pat, why in FORBIDDEN:
                if pat.search(line):
                    violations.append(f"{mod}:{lineno}: {line.strip()}\n    -> {why}")
    assert not violations, (
        "forbidden cyclotomic-matrix patterns found (Sage bug guard):\n"
        + "\n".join(violations)
    )


def test_no_bare_assert_on_verification_path() -> None:
    """The paper's stated guarantees must survive `python -O`: verification
    modules use types.require(), never bare assert statements."""
    offenders = []
    all_modules = [
        str(p.relative_to(PKG)) for p in PKG.rglob("*.py") if "tests" not in p.parts
    ]
    for mod in all_modules:
        src = (PKG / mod).read_text()
        for lineno, line in enumerate(src.splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("assert ") and not stripped.startswith(
                "assert not violations"
            ):
                offenders.append(f"{mod}:{lineno}: {stripped}")
    assert not offenders, (
        "bare `assert` on a verification path (stripped under python -O); "
        "use types.require():\n" + "\n".join(offenders)
    )
