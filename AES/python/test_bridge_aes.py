#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
End-to-end bridge test (report.tex §3.1) with REAL 7-round AES.

For a random AES-128 key:
  - build a δ-set of 256 states at x_1 with x_1[0] = 0..255, rest fixed random;
  - encrypt forward through rounds 1..6 to ciphertexts C^(ω);
  - record the true intermediate a_ω := x_5[0]^(ω);
  - ONLINE: peel round 6 with the correct k_6[0,7,10,13], compute
        v_ω = MC^{-1}(x_6[col0])[0],   g_ω = (L^{-1}(v_0 ⊕ v_ω))^{-1},
    and fingerprint I({g_ω});
  - OFFLINE: d_ω = a_0 ⊕ a_ω, fingerprint I({d_ω^{-1}});
  - assert v_ω == S(a_ω) ⊕ u_5[0]   (bookkeeping sanity)
  - assert I_online == I_offline under exactly one parity (raw / append-0).

Everything vectorised over the 256 δ-set elements.
"""

import sys
import numpy as np
from test_invariants import EXP, LOG, gf_mul, gf_inv, gf_pow
from test_zeroinv import fingerprint_fb_single

# ---------------------------------------------------------------- AES parts

def _sbox_tables():
    inv = np.zeros(256, dtype=np.uint8)
    inv[1:] = gf_inv(np.arange(1, 256, dtype=np.uint8))     # 0^{-1}=0
    # affine: rotate-left XOR (standard AES)
    def aff(b):
        r = b
        for sh in (1, 2, 3, 4):
            r ^= ((b << sh) | (b >> (8 - sh))) & 0xFF
        return r ^ 0x63
    L = np.array([aff(x) for x in range(256)], dtype=np.uint8)   # L(x)=A·x⊕63
    # we want the *linear* part A alone, and its inverse, as LUTs
    A  = L ^ 0x63                                                # A[x]=A·x
    Ai = np.zeros(256, dtype=np.uint8)
    Ai[A] = np.arange(256, dtype=np.uint8)                       # inverse perm
    Sb  = L[inv]                                                 # S(t)=A·t^{-1}⊕63
    Sbi = np.zeros(256, dtype=np.uint8)
    Sbi[Sb] = np.arange(256, dtype=np.uint8)
    return Sb, Sbi, A, Ai, inv

SBOX, SBOX_INV, A_LIN, A_LIN_INV, FIELD_INV = _sbox_tables()
# sanity vs known values
assert SBOX[0x00] == 0x63 and SBOX[0x53] == 0xED and SBOX_INV[0x63] == 0x00

RCON = [0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1B,0x36]

def key_schedule(key16):
    """AES-128 key expansion -> list of 11 (or more) 16-byte round keys.
       We use rk[0]=k_{-1} (whitening), rk[1]=k_0, ..., rk[7]=k_6."""
    w = list(key16)
    for i in range(4, 4 * 11):
        t = w[4*(i-1):4*i]
        if i % 4 == 0:
            t = [SBOX[t[1]], SBOX[t[2]], SBOX[t[3]], SBOX[t[0]]]
            t[0] ^= RCON[i//4 - 1]
        w += [w[4*(i-4)+j] ^ t[j] for j in range(4)]
    rk = [np.array(w[16*r:16*(r+1)], dtype=np.uint8) for r in range(11)]
    return rk

# state layout: byte i is (row i%4, col i//4). Operate on (N,16) arrays.

SR_IDX  = np.array([ (r + 4*((c + r) % 4)) for c in range(4) for r in range(4)
                   ], dtype=np.intp).reshape(4,4).T.flatten()
# Derive: z[r,c] = y[r,(c+r)%4]  =>  z[idx(r,c)] = y[idx(r,(c+r)%4)]
# idx(r,c)=r+4c. Build SR_IDX s.t. z = y[SR_IDX].
def _build_sr():
    idx = np.empty(16, dtype=np.intp)
    for r in range(4):
        for c in range(4):
            idx[r + 4*c] = r + 4*((c + r) % 4)
    return idx
SR_IDX = _build_sr()
SRI_IDX = np.argsort(SR_IDX)

MC_MAT  = np.array([[2,3,1,1],[1,2,3,1],[1,1,2,3],[3,1,1,2]], dtype=np.uint8)
MCI_MAT = np.array([[14,11,13,9],[9,14,11,13],[13,9,14,11],[11,13,9,14]],
                   dtype=np.uint8)

def mix_columns(s, M=MC_MAT):
    """s: (N,16) -> (N,16)."""
    s4 = s.reshape(-1, 4, 4)                       # [N, col, row]? careful
    # byte i=(r,c)=r+4c, so s.reshape(N,4,4) with axis1=c, axis2=r? Let's do explicit.
    out = np.empty_like(s)
    for c in range(4):
        col = s[:, 4*c:4*c+4]                      # (N,4) rows 0..3
        for r in range(4):
            acc = np.zeros(s.shape[0], dtype=np.uint8)
            for k in range(4):
                acc ^= gf_mul(int(M[r, k]), col[:, k])
            out[:, 4*c + r] = acc
    return out

def aes_round(s, rk, last=False):
    s = SBOX[s]
    s = s[:, SR_IDX]
    if not last:
        s = mix_columns(s)
    return s ^ rk

def encrypt_from_x1(x1, rk):
    """x1: (N,16). rk: list from key_schedule. Returns (C, x5)."""
    s = x1.copy()
    x5 = None
    # rounds 1..6 (round 6 is last, no MC). rk index: round i uses rk[i+1].
    for rnd in range(1, 7):
        if rnd == 5:
            x5 = s.copy()                          # x_5 = input to round 5
        s = aes_round(s, rk[rnd + 1], last=(rnd == 6))
    return s, x5

# ------------------------------------------------------------ the test

def run(T=200, seed=0):
    rng = np.random.default_rng(seed)
    n_skip_s0 = n_raw = n_app = n_both = n_neither = 0
    bookkeeping_ok = True

    for trial in range(T):
        key = rng.integers(0, 256, size=16, dtype=np.uint8)
        rk  = key_schedule(key)
        # δ-set at x_1
        base = rng.integers(0, 256, size=16, dtype=np.uint8)
        x1 = np.tile(base, (256, 1))
        x1[:, 0] = np.arange(256, dtype=np.uint8)
        # encrypt
        C, x5 = encrypt_from_x1(x1, rk)
        a = x5[:, 0].astype(np.uint8)              # a_ω, ω=0..255
        s = int(a[0])
        if s == 0:
            n_skip_s0 += 1
            continue

        # --- ONLINE: peel round 6 on positions {0,7,10,13} ----------------
        k6 = rk[7]
        z6 = C ^ k6                                # (256,16)
        # y_6[col0] = z_6[0,13,10,7]
        y6c0 = z6[:, [0, 13, 10, 7]]
        x6c0 = SBOX_INV[y6c0]                      # (256,4) = x_6[0..3]
        # v = MC^{-1}(x_6 col0)[row 0] = 14·x6[0]⊕11·x6[1]⊕13·x6[2]⊕9·x6[3]
        v = (gf_mul(14, x6c0[:, 0]) ^ gf_mul(11, x6c0[:, 1]) ^
             gf_mul(13, x6c0[:, 2]) ^ gf_mul( 9, x6c0[:, 3]))
        # bookkeeping check: v_ω == S(a_ω) ⊕ u_5[0]
        u5_0 = (gf_mul(14, rk[6][0]) ^ gf_mul(11, rk[6][1]) ^
                gf_mul(13, rk[6][2]) ^ gf_mul( 9, rk[6][3]))
        if not np.array_equal(v, SBOX[a] ^ u5_0):
            bookkeeping_ok = False
            print(f"  !! trial {trial}: v != S(a)⊕κ  "
                  f"(mismatch at {int((v != (SBOX[a]^u5_0)).sum())} of 256)")
        # g_ω = (L^{-1}(v_0 ⊕ v_ω))^{-1}, ω=1..255, with 0^{-1}=0
        dv = v[0] ^ v[1:]                          # (255,)
        g  = FIELD_INV[A_LIN_INV[dv]]              # (255,)

        # --- OFFLINE: d_ω^{-1} ------------------------------------------
        d = a[0] ^ a[1:]                           # (255,)
        dinv = FIELD_INV[d]

        # fingerprints, both parities
        IG  = fingerprint_fb_single(g)
        ID  = fingerprint_fb_single(dinv)
        IGp = fingerprint_fb_single(np.concatenate([g,    [np.uint8(0)]]))
        IDp = fingerprint_fb_single(np.concatenate([dinv, [np.uint8(0)]]))

        raw = np.array_equal(IG,  ID)
        app = np.array_equal(IGp, IDp)
        if raw and app: n_both += 1
        elif raw:       n_raw  += 1
        elif app:       n_app  += 1
        else:           n_neither += 1

    print(f"[*] {T} random keys tested")
    print(f"    bookkeeping v=S(a)⊕κ : {'OK' if bookkeeping_ok else 'FAIL'}")
    print(f"    skipped (s=0)        : {n_skip_s0}")
    print(f"    match raw only       : {n_raw}")
    print(f"    match append-0 only  : {n_app}")
    print(f"    match BOTH           : {n_both}")
    print(f"    match NEITHER        : {n_neither}   <-- bridge failures")
    ok = (n_neither == 0) and bookkeeping_ok
    print(f"[*] {'PASS' if ok else 'FAIL'}: end-to-end bridge with real AES")
    return ok

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    sys.exit(0 if run(T=T) else 1)
