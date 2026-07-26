# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""χ-CANON: the AGL(1,GF(256))-invariant fingerprint (Layer 5).
This is the "scale-and-shift-proof hash" that lets online match offline
without knowing the bridging byte a."""

import hashlib
from gf256 import mul, inv, gpow

def _permute_chi(chi, a, b):
    """Permute a 256-bit int: bit d → bit (a·d ⊕ b)."""
    out = 0
    for d in range(256):
        if (chi >> d) & 1:
            out |= 1 << (mul(a, d) ^ b)
    return out

def chi_canon(seq):
    """seq[0..255] with seq[0]=0.  Return 12-byte AGL-invariant fingerprint.

    NAIVE-CORRECT version: canonical form = min over all 255×256 affine maps.
    This is ~65K ops/call — obviously right, obviously slow.  The production
    version (C_post=8–15, ../A1-e2e/chi256_fast.h) computes (α*,β*) directly
    from power-sums S₁,S₃,S₅,S₇ via a moving-frame formula, skipping the
    255× loop.  Same output; this version is for verification."""
    # Step 1: build χ = 256-bit parity indicator of the multiset {seq[1..255]}.
    chi = 0
    for w in range(1, 256):
        chi ^= 1 << seq[w]
    # Steps 2–4: canonical = lexicographic min over the full AGL(1,256) orbit.
    canon = min(_permute_chi(chi, a, b)
                for a in range(1, 256) for b in range(256))
    return hashlib.blake2b(canon.to_bytes(32, 'little'), digest_size=12).digest()

def chi_canon_fast_stub(seq):
    """Sketch of the C_post=15 version (see chi256_fast.h for working code).
    Shown here so the HTML's §2F pseudocode has a runnable counterpart."""
    chi = 0; S1 = S3 = S5 = S7 = 0
    for w in range(1, 256):
        d = seq[w]
        chi ^= 1 << d
        S1 ^= d; S3 ^= gpow(d,3); S5 ^= gpow(d,5); S7 ^= gpow(d,7)
    # α*,β* derivation (exact polynomial in chi256_fast.h; placeholder here):
    # alpha = frame_alpha(S1,S3,S5,S7); beta = frame_beta(S1,S3,alpha)
    # return hash(min(_permute_chi(chi,alpha,beta), _permute_chi(chi,alpha,beta^1)))
    raise NotImplementedError("use chi_canon() for correctness; C code for speed")

def verify_agl_invariance(trials=5):
    """The bridge property: offline has d=1/e; online has h=a⊕a²·d (unknown a).
    χ-CANON must give the same output on {d_w} and {h_w}."""
    import random
    for _ in range(trials):
        e = [0] + [random.randrange(1,256) for _ in range(255)]
        d = [0] + [inv(e[w]) for w in range(1,256)]          # offline side
        a = random.randrange(1,256); a2 = mul(a,a)
        h = [0] + [a ^ mul(a2, d[w]) for w in range(1,256)]  # online side
        assert chi_canon(d) == chi_canon(h), f"FAIL at a={a}"
    return True

if __name__ == "__main__":
    assert verify_agl_invariance(5)
    print("chi_canon AGL-invariance: 5/5 OK (naive canon ~13s/call; raise trials if patient)")
