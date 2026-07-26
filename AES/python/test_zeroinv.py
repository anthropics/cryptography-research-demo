#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Empirical check of the 0^{-1} parity patch (report.tex lines 651-672).

For random 256-multisets A over GF(256) with a_0 != 0:
  online  G  = {(s^-1 xor a_w^-1)^-1 : w=1..255}   (0^{-1}:=0 convention)
  offline D  = {(s xor a_w)^-1       : w=1..255}
Claim: with G' = G ∪ {0}, D' = D ∪ {0},
  exactly one of  I(G)==I(D)  or  I(G')==I(D')  holds,
and which one is determined by the parity of b = #{w : a_w in {0,s}}.
"""

import numpy as np
from test_invariants import (gf_inv, pairwise_power_sums, ALL_EXP,
                             gf_pow, gf_div, N_SET)

def fingerprint_fb_single(g):
    """g: 1-D uint8 array. Returns (12,) uint8 fingerprint with fallback."""
    P = pairwise_power_sums(g[None, :])[0]         # (15,)
    nz = np.where(P != 0)[0]
    I = np.zeros(12, dtype=np.uint8)
    if len(nz) < 13:
        return I                                    # sentinel
    m_pos, n_pos = nz[0], nz[1:13]
    m = ALL_EXP[m_pos]; Pm = int(P[m_pos])
    for j, npos in enumerate(n_pos):
        n = ALL_EXP[npos]; Pn = int(P[npos])
        I[j] = gf_div(gf_pow(Pm, n), gf_pow(Pn, m))
    return I

def run(N=3000, seed=0):
    rng = np.random.default_rng(seed)
    n_match_raw = n_match_app = n_match_both = n_match_neither = 0
    parity_consistent = 0
    bad_hist = np.zeros(32, dtype=int)

    for trial in range(N):
        A = rng.integers(0, 256, size=256, dtype=np.uint8)
        if A[0] == 0:
            A[0] = rng.integers(1, 256)            # enforce s != 0
        s = int(A[0])
        a = A[1:]                                   # (255,)

        # bad-index count
        b = int(((a == 0) | (a == s)).sum())
        bad_hist[min(b, 31)] += 1

        # online: g_w = (s^-1 xor a_w^-1)^-1  with 0^{-1}=0
        s_inv = int(gf_inv(s))
        g = gf_inv(np.uint8(s_inv) ^ gf_inv(a))     # (255,)
        # offline: d_w^-1 = (s xor a_w)^-1
        dinv = gf_inv(np.uint8(s) ^ a)              # (255,)

        IG  = fingerprint_fb_single(g)
        ID  = fingerprint_fb_single(dinv)
        IGp = fingerprint_fb_single(np.concatenate([g,    [np.uint8(0)]]))
        IDp = fingerprint_fb_single(np.concatenate([dinv, [np.uint8(0)]]))

        raw = np.array_equal(IG,  ID)
        app = np.array_equal(IGp, IDp)
        if raw and app:     n_match_both    += 1
        elif raw:           n_match_raw     += 1
        elif app:           n_match_app     += 1
        else:               n_match_neither += 1

        # parity prediction: raw should match iff b even
        if raw == (b % 2 == 0) and app == (b % 2 == 1):
            parity_consistent += 1

    print(f"[*] {N} trials, s != 0 enforced")
    print(f"    bad-index count histogram (b=0..): "
          f"{bad_hist[:10].tolist()} ...")
    print(f"    match raw only   : {n_match_raw}")
    print(f"    match append only: {n_match_app}")
    print(f"    match BOTH       : {n_match_both}")
    print(f"    match NEITHER    : {n_match_neither}   <-- counterexamples")
    print(f"    parity-prediction correct: {parity_consistent}/{N}")
    if n_match_neither == 0:
        print("[*] PASS: paper's parity patch holds on all trials.")
    else:
        print("[*] FAIL: found counterexample(s).")

    # Also test s==0 explicitly to confirm it's a genuine failure mode.
    print("\n[*] spot-check s=0 (paper says this fails, absorbed as 256/255):")
    A = rng.integers(0, 256, size=256, dtype=np.uint8); A[0] = 0
    a = A[1:]
    g = gf_inv(np.uint8(0) ^ gf_inv(a))             # = a  (since 0^-1=0)
    dinv = gf_inv(np.uint8(0) ^ a)                  # = a^{-1}
    IG = fingerprint_fb_single(g)
    ID = fingerprint_fb_single(dinv)
    IGp = fingerprint_fb_single(np.concatenate([g,    [np.uint8(0)]]))
    IDp = fingerprint_fb_single(np.concatenate([dinv, [np.uint8(0)]]))
    print(f"    raw match: {np.array_equal(IG,ID)}, "
          f"append match: {np.array_equal(IGp,IDp)}  (expect both False)")

if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    run(N=N)
