# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

"""Wrappers around the HAWK reference implementation for hpub and sign/verify.

`hpub` is computed natively (SHAKE256 of pk bytes, spec §2.2); `signverify`
shells out to the bundled `hawk_shim` binary built from the NIST reference
implementation.
"""

from __future__ import annotations

import hashlib
import subprocess
import tempfile
from pathlib import Path

from hawk_recover.params import HawkParams
from hawk_recover.types import require

_SHIM = Path(__file__).parent.parent / "ref_impl" / "hawk_shim"


def hpub(pk: bytes, params: HawkParams) -> bytes:
    return hashlib.shake_256(pk).digest(params.hpub_bytes)


def signverify(pk: bytes, sk_full: bytes, msg: bytes = b"test") -> bool:
    """Sign `msg` with sk_full (core‖hpub) and verify against pk via hawk_shim.

    Returns False only for a genuine key failure (shim rc 2 SIGN FAILED /
    rc 3 VERIFY FAILED); a shim that could not run (rc 1: wrong input sizes,
    usage) is a defect in our inputs and raises rather than reading as "wrong
    key"."""
    require(_SHIM.exists(), f"hawk_shim missing at {_SHIM}")
    with tempfile.TemporaryDirectory() as d:
        pkp, skp = Path(d) / "pk.bin", Path(d) / "sk.bin"
        pkp.write_bytes(pk)
        skp.write_bytes(sk_full)
        r = subprocess.run(
            [str(_SHIM), "signverify", str(pkp), str(skp), msg.decode("latin-1")],
            capture_output=True,
            text=True,
        )
    out = (r.stdout + r.stderr).strip()
    if r.returncode == 0 and "PASS" in out.upper():
        return True
    if r.returncode in (2, 3):
        return False
    raise RuntimeError(f"hawk_shim signverify failed to run (rc={r.returncode}): {out}")


def shim_hpub(pk: bytes) -> bytes:
    """Cross-check: compute hpub via the shim (used in tests)."""
    require(_SHIM.exists(), f"hawk_shim missing at {_SHIM}")
    with tempfile.TemporaryDirectory() as d:
        pkp = Path(d) / "pk.bin"
        pkp.write_bytes(pk)
        # The shim's argc check requires a third positional even for hpub.
        r = subprocess.run(
            [str(_SHIM), "hpub", str(pkp), "/dev/null"],
            capture_output=True,
        )
    require(r.returncode == 0, f"shim hpub failed: {r.stderr!r}")
    return r.stdout
