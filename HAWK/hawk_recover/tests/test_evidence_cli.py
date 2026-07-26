# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Guards for the evidence CLI — the third-party verification boundary."""

import pytest

from hawk_recover.evidence import verify
from hawk_recover.evidence.__main__ import main as cli_main


def test_verify_fails_at_evidence_loading_not_import() -> None:
    """Verifying a nonexistent evidence directory must fail at evidence
    loading (FileNotFoundError) — NOT at import time. An ImportError here
    means a package reorganization broke the tool's deferred imports."""
    with pytest.raises(FileNotFoundError):
        verify.main(["/nonexistent/evidence"])


def test_dispatcher_routes_and_rejects() -> None:
    assert cli_main([]) == 2
    assert cli_main(["bogus"]) == 2
    with pytest.raises(FileNotFoundError):
        cli_main(["verify", "/nonexistent/evidence"])
