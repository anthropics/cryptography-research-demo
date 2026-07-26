#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Clumping / collision-entropy test for the I_{m,n} fingerprint on the REAL
offline distribution: {d_ω^{-1}} where d_ω = x_5[0]^(ω) ⊕ x_5[0]^(0) over a
δ-set pushed through 4 actual AES rounds (x_1 -> x_5).

Sampling: random x_1 base state + random independent round-key bytes k_1..k_4
(79 free bytes total). DFJ's 11 parameters are a compression of this; we
sample the full preimage, which induces the same distribution on {d_ω}.

Outputs the same prefix-H2 curve as test_clumping.py, plus the P_m=0 rate,
so we can see whether AES structure changes anything vs uniform multisets.
"""

import sys, time
import numpy as np
from collections import Counter

from test_invariants import gf_mul, gf_inv, ALL_EXP
from test_bridge_aes import SBOX, SR_IDX
from test_clumping_fast import (POW_ODD, SQ, MUL, SUBMASK_PAIRS,
                                EXPS, LOG32, EXP)

# -------------------------------------------------- fast AES round (batched)

Mc = {c: gf_mul(c, np.arange(256, dtype=np.uint8)) for c in (2, 3)}

def mix_columns_fast(s):
    """s: (N,16) uint8 -> (N,16). Uses 2x/3x LUTs only."""
    out = np.empty_like(s)
    for c in range(4):
        a = s[:, 4*c+0]; b = s[:, 4*c+1]
        d = s[:, 4*c+2]; e = s[:, 4*c+3]
        a2, b2, d2, e2 = Mc[2][a], Mc[2][b], Mc[2][d], Mc[2][e]
        out[:, 4*c+0] = a2 ^ b2 ^ b ^ d ^ e         # 2a⊕3b⊕c⊕d
        out[:, 4*c+1] = a ^ b2 ^ d2 ^ d ^ e         # a⊕2b⊕3c⊕d
        out[:, 4*c+2] = a ^ b ^ d2 ^ e2 ^ e         # a⊕b⊕2c⊕3d
        out[:, 4*c+3] = a2 ^ a ^ b ^ d ^ e2         # 3a⊕b⊕c⊕2d
    return out

def aes_rounds_1to4(x1, rk):
    """x1: (N,16); rk: (N,4,16) independent per-sample round keys k_1..k_4.
       Returns x5: (N,16). Full rounds (SB,SR,MC,ARK) ×4."""
    s = x1
    for r in range(4):
        s = SBOX[s]
        s = s[:, SR_IDX]
        s = mix_columns_fast(s)
        s = s ^ rk[:, r, :]
    return s

# ----------------------------------------- δ-set -> {d_ω^{-1}} -> fingerprint

def sample_dinv(B, rng):
    """Return dinv: (B,255) uint8 and a: (B,256) for diagnostics."""
    base = rng.integers(0, 256, size=(B, 16), dtype=np.uint8)
    rk   = rng.integers(0, 256, size=(B, 4, 16), dtype=np.uint8)
    # δ-set: replicate each sample 256×, vary byte 0
    x1 = np.repeat(base, 256, axis=0)                      # (B*256,16)
    x1[:, 0] = np.tile(np.arange(256, dtype=np.uint8), B)
    rkE = np.repeat(rk, 256, axis=0)                       # (B*256,4,16)
    x5 = aes_rounds_1to4(x1, rkE)
    a  = x5[:, 0].reshape(B, 256)                          # a_ω
    d  = a[:, :1] ^ a[:, 1:]                               # (B,255)
    return gf_inv(d), a

def fingerprint_batch(dinv):
    """dinv: (B,255) -> I: (B,12) uint8, P: (B,15) uint8. |V|=255 odd path."""
    B = dinv.shape[0]
    p_odd = np.zeros((B, 64), dtype=np.uint8)
    for w in range(255):
        p_odd ^= POW_ODD[dinv[:, w]].view(np.uint8).reshape(-1, 64)
    p = np.zeros((B, 128), dtype=np.uint8)
    p[:, 1::2] = p_odd
    for r in range(1, 64):
        p[:, 2*r] = SQ[p[:, r]]
    P = np.empty((B, 15), dtype=np.uint8)
    for j, m in enumerate(ALL_EXP):
        acc = p[:, m].copy()
        for (k, mk) in SUBMASK_PAIRS[m]:
            acc ^= MUL[p[:, k], p[:, mk]]
        P[:, j] = acc
    # fallback + log-domain ratios (same as test_clumping_fast)
    rank = np.where(P != 0, np.arange(15, dtype=np.int8)[None, :], np.int8(15))
    idx = np.argsort(rank, axis=1, kind="stable")[:, :13]
    selP = np.take_along_axis(P, idx, axis=1)
    selE = EXPS[idx]
    sentinel = selP[:, 12] == 0
    logP = LOG32[selP]
    m0, lPm0 = selE[:, 0], logP[:, 0]
    nj, lPnj = selE[:, 1:13], logP[:, 1:13]
    e = (lPm0[:, None] * nj - lPnj * m0[:, None]) % 255
    I = EXP[e]
    I[sentinel] = 0
    return I, P, sentinel

# ------------------------------------------------------------------- driver

def run(N=1_000_000, chunk=1024, seed=0):
    rng = np.random.default_rng(seed)
    print(f"[*] sampling N={N:,} δ-sets through 4 real AES rounds "
          f"(chunk={chunk}) ...")
    K = 6
    prefix_counters = [Counter() for _ in range(K)]
    full = Counter()
    nz_hist = np.zeros(16, dtype=np.int64)
    n_sentinel = 0
    t0 = time.time()
    done = 0
    while done < N:
        B = min(chunk, N - done)
        dinv, a = sample_dinv(B, rng)
        I, P, sen = fingerprint_batch(dinv)
        nz_hist += np.bincount((P != 0).sum(axis=1), minlength=16)
        n_sentinel += int(sen.sum())
        I = I[~sen]
        for row in I:
            full[bytes(row)] += 1
        for k in range(1, K + 1):
            for row in I[:, :k]:
                prefix_counters[k-1][bytes(row)] += 1
        done += B
        if done % (chunk * 32) == 0 or done == N:
            rate = done / (time.time() - t0)
            print(f"    {done:,}/{N:,}  {rate:.0f}/s  "
                  f"ETA {(N-done)/rate/60:.1f} min", end="\r")
    print()
    M = N - n_sentinel
    print(f"[*] done in {time.time()-t0:.1f}s; sentinels={n_sentinel}; M={M:,}")
    print(f"[*] nonzero-P_m count histogram (0..15): {nz_hist.tolist()}")
    print(f"    (cf. uniform: ~{N*(1-2**-8)**15:.0f} at 15, "
          f"~{N*2**-15.2:.1f} at <13)")

    # full
    distinct = len(full); maxc = max(full.values()) if full else 0
    print(f"\n[*] FULL 12-byte: distinct={distinct:,}  "
          f"M-distinct={M-distinct}  max mult={maxc}")
    if maxc > 1:
        for fp, c in full.most_common(5):
            print(f"    {fp.hex()}  x{c}")

    # prefix curve
    print(f"\n[*] PREFIX H2 curve (real-AES {{d_ω}} distribution):")
    print(f"    k   collisions C      H2_hat    ideal")
    for k in range(1, K + 1):
        C = sum(v*(v-1)//2 for v in prefix_counters[k-1].values())
        H2 = np.log2(M*(M-1)/(2*C)) if C else float("inf")
        print(f"    {k}   {C:>12,}   {H2:7.3f}   {k*np.log2(255):7.3f}")

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
    run(N=N, chunk=chunk)
