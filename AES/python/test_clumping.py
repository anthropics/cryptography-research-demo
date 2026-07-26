#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Clumping / collision-entropy test for the I_{m,n} fingerprint on
uniform-random 255-multisets over GF(256).

Goal: check that the *fingerprint construction itself* (the "novel piece")
does not concentrate mass, independent of the DFJ parameter structure.

  - Fast fingerprint via single-index power sums p_r (paper §3.2 Reduction 1),
    self-checked against the O(255^2) reference.
  - Stream N samples in chunks; maintain exact Counter on the full 12-byte
    fingerprint (distinct-count + heavy hitters) and per-prefix collision
    counts for k=1..6 bytes (collision-entropy curve).
"""

import sys
import numpy as np
from collections import Counter

from test_invariants import (EXP, LOG, gf_mul, gf_pow, gf_div,
                             ALL_EXP, N_SET, pairwise_power_sums)

# ------------------------------------------------------------------ fast P_m

MAX_M = max(ALL_EXP)                         # 127
R_ALL = np.arange(1, MAX_M + 1)              # we need p_1..p_127

# POW_R[v, r-1] = v^r  for v=0..255, r=1..127   -> (256, 127) uint8
POW_R = np.stack([gf_pow(np.arange(256, dtype=np.uint8), int(r))
                  for r in R_ALL], axis=1)

# For each m in ALL_EXP, precompute the list of submask pairs (k, m-k)
# with 0 < k < m-k and k ⊂ m (bitwise).
SUBMASK_PAIRS = {}
for m in ALL_EXP:
    pairs = []
    for k in range(1, m):
        if (k & m) == k and k < (m - k):
            pairs.append((k, m - k))
    SUBMASK_PAIRS[m] = pairs


def power_sums_fast(g):
    """
    g: (B, 255) uint8.  Returns P: (B, len(ALL_EXP)) uint8 via the p_r route.
    |V|=255 odd, so  P_m = XOR_{0<k<m-k, k⊂m} p_k * p_{m-k}  XOR  p_m.
    """
    B = g.shape[0]
    # parity histogram  h[b, v] = (# of w with g[b,w]==v) mod 2
    h = np.zeros((B, 256), dtype=bool)
    rows = np.arange(B)
    for w in range(g.shape[1]):
        h[rows, g[:, w]] ^= True
    # p[b, r-1] = XOR_{v: h[b,v]} v^r   -> (B, 127)
    #   = GF(2)-matrix product h (B,256) · POW_R (256,127) with XOR-add.
    # Realised as select-rows-then-XOR, per sample. Chunk B kept small by caller.
    # Trick: use np.where mask per row — but variable length. Simpler:
    #   p = XOR over v of h[:,v:v+1] * POW_R[v]  (broadcast), 256 iters.
    p = np.zeros((B, MAX_M), dtype=np.uint8)
    for v in range(256):
        mask = h[:, v]
        if mask.any():
            p[mask] ^= POW_R[v]
    # assemble P_m for each m
    P = np.zeros((B, len(ALL_EXP)), dtype=np.uint8)
    for j, m in enumerate(ALL_EXP):
        acc = p[:, m - 1].copy()                      # the lone p_m term
        for (k, mk) in SUBMASK_PAIRS[m]:
            acc ^= gf_mul(p[:, k - 1], p[:, mk - 1])
        P[:, j] = acc
    return P


def fingerprints_fast(g):
    """g: (B,255) -> I: (B,12) with the 15-exponent fallback rule."""
    P = power_sums_fast(g)
    B = g.shape[0]
    I = np.zeros((B, 12), dtype=np.uint8)
    exps = np.array(ALL_EXP, dtype=np.int64)
    core = P[:, :13]
    fast = np.all(core != 0, axis=1)
    if fast.any():
        Pm = core[fast, 0]
        for j, n in enumerate(N_SET):
            I[fast, j] = gf_div(gf_pow(Pm, n), gf_pow(core[fast, 1 + j], 7))
    for r in np.where(~fast)[0]:
        nz = np.where(P[r] != 0)[0]
        if len(nz) < 13:
            continue
        m = int(exps[nz[0]]); Pm = int(P[r, nz[0]])
        for j, npos in enumerate(nz[1:13]):
            n = int(exps[npos]); Pn = int(P[r, npos])
            I[r, j] = gf_div(gf_pow(Pm, n), gf_pow(Pn, m))
    return I, P


# --------------------------------------------------------------- self-check

def self_check(rng, B=256):
    g = rng.integers(0, 256, size=(B, 255), dtype=np.uint8)
    P_ref = pairwise_power_sums(g)[:, :len(ALL_EXP)]
    P_fast = power_sums_fast(g)
    ok = np.array_equal(P_ref, P_fast)
    print(f"[*] self-check fast P_m vs reference on {B} samples: "
          f"{'OK' if ok else 'FAIL'}")
    if not ok:
        diff = (P_ref != P_fast)
        print(f"    mismatched entries: {int(diff.sum())} / {diff.size}")
        sys.exit(1)


# --------------------------------------------------------------------- main

def run(N=1_000_000, chunk=4096, seed=0):
    rng = np.random.default_rng(seed)
    self_check(rng)

    print(f"[*] streaming N={N:,} random 255-multisets, chunk={chunk} ...")
    full_counter = Counter()
    sentinel = 0
    # prefix collision counts for k=1..6 bytes
    K = 6
    prefix_counters = [Counter() for _ in range(K)]

    done = 0
    while done < N:
        B = min(chunk, N - done)
        g = rng.integers(0, 256, size=(B, 255), dtype=np.uint8)
        I, _ = fingerprints_fast(g)
        # drop sentinels (<13 nonzero P_m)
        keep = np.any(I != 0, axis=1)
        sentinel += int((~keep).sum())
        I = I[keep]
        # full 12-byte
        for row in I:
            full_counter[bytes(row)] += 1
        # prefixes
        for k in range(1, K + 1):
            pk = I[:, :k]
            for row in pk:
                prefix_counters[k - 1][bytes(row)] += 1
        done += B
        if done % (chunk * 25) == 0 or done == N:
            print(f"    {done:,}/{N:,}", end="\r")
    print()

    M = N - sentinel
    print(f"[*] sentinel (<13 nonzero P_m) rows dropped: {sentinel} "
          f"(expect ~{N*2**-15.2:.1f})")
    print(f"[*] effective sample size M = {M:,}")

    # ---- full-fingerprint clumping ---------------------------------------
    distinct = len(full_counter)
    maxcount = max(full_counter.values())
    coll = sum(v * (v - 1) // 2 for v in full_counter.values())
    exp96 = M * (M - 1) / 2 / (255.0 ** 12)
    print(f"\n[*] FULL 12-byte fingerprint:")
    print(f"    distinct values : {distinct:,}  (M - distinct = {M-distinct})")
    print(f"    max multiplicity: {maxcount}")
    print(f"    pair-collisions : {coll}   "
          f"(birthday @95.9 bits expects {exp96:.2e})")
    if maxcount > 1:
        print(f"    top-5 heavy fingerprints:")
        for fp, c in full_counter.most_common(5):
            print(f"      {fp.hex()}  x{c}")

    # ---- prefix collision-entropy curve ----------------------------------
    print(f"\n[*] PREFIX collision-entropy curve "
          f"(H2_hat = log2(M(M-1)/(2C))):")
    print(f"    k  collisions C      H2_hat    ideal 8k*log2(255)/8")
    for k in range(1, K + 1):
        C = sum(v * (v - 1) // 2 for v in prefix_counters[k - 1].values())
        if C > 0:
            H2 = np.log2(M * (M - 1) / (2 * C))
        else:
            H2 = float("inf")
        ideal = k * np.log2(255)
        print(f"    {k}  {C:>12,}   {H2:7.2f}    {ideal:6.2f}")

    print("\n[*] done.")


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 4096
    run(N=N, chunk=chunk)
