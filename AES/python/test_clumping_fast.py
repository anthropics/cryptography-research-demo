#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
High-throughput clumping test for the I_{m,n} fingerprint.

Target: N ~ 3e8 random 255-multisets, enough for birthday signal at k=6
(~48 bits).  Strategy:
  - compute only the 64 odd-index power sums p_r via a single packed table,
    expand even indices by squaring, assemble all 15 P_m with a 256x256
    GF-mul LUT;
  - fallback fingerprint as before;
  - pack first 8 bytes of each fingerprint into one uint64, accumulate into
    a preallocated array;
  - one np.sort at the end, then count collisions for k=5,6,7 by right-shift.

All hot paths are numpy-vectorised (XOR/gather); no per-row Python.
"""

import sys, time
import numpy as np

from test_invariants import gf_pow, gf_div, ALL_EXP, N_SET
from test_clumping import SUBMASK_PAIRS, power_sums_fast

# ------------------------------------------------------------------ tables

ODD_R = np.arange(1, 128, 2)                               # 64 odd exponents
POW_ODD = np.stack([gf_pow(np.arange(256, dtype=np.uint8), int(r))
                    for r in ODD_R], axis=1)               # (256, 64) uint8
POW_ODD = np.ascontiguousarray(POW_ODD)

SQ  = gf_pow(np.arange(256, dtype=np.uint8), 2)            # (256,) v^2
MUL = np.empty((256, 256), dtype=np.uint8)                 # full GF-mul LUT
for a in range(256):
    MUL[a] = gf_pow(np.arange(256, dtype=np.uint8), 1)     # placeholder
# build MUL properly via log/exp (reuse gf_mul semantics without its overhead)
from test_invariants import EXP, LOG
for a in range(256):
    if a == 0:
        MUL[a] = 0
    else:
        la = int(LOG[a])
        row = np.zeros(256, dtype=np.uint8)
        nz = np.arange(1, 256)
        row[nz] = EXP[(la + LOG[nz]) % 255]
        MUL[a] = row

EXPS = np.array(ALL_EXP, dtype=np.int64)                   # (15,)
IDX_ODD = {int(r): i for i, r in enumerate(ODD_R)}

# -------------------------------------------------------- fingerprint kernel

LOG32 = LOG.astype(np.int64)                               # widen once

def fingerprint_chunk(g):
    """
    g: (B,255) uint8  ->  keys: (B',) uint64  (first 8 fp bytes, big-endian),
    where B' = B minus sentinels.  Also returns sentinel count.
    Fully vectorised; no per-row Python.
    """
    B = g.shape[0]

    # p_odd[b, i] = XOR_w g[b,w]^ODD_R[i]
    p_odd = np.zeros((B, 64), dtype=np.uint8)
    for w in range(255):
        p_odd ^= POW_ODD[g[:, w]]

    # expand to p[b, r] for r=0..127
    p = np.zeros((B, 128), dtype=np.uint8)
    p[:, ODD_R] = p_odd
    for r in range(1, 64):
        p[:, 2 * r] = SQ[p[:, r]]

    # assemble P_m for the 15 exponents
    P = np.empty((B, 15), dtype=np.uint8)
    for j, m in enumerate(ALL_EXP):
        acc = p[:, m].copy()
        for (k, mk) in SUBMASK_PAIRS[m]:
            acc ^= MUL[p[:, k], p[:, mk]]
        P[:, j] = acc

    # --- vectorised fallback: pick first 13 nonzero P's per row -------------
    # rank zeros to the back, keep original order among nonzeros
    rank = np.where(P != 0,
                    np.arange(15, dtype=np.int8)[None, :],
                    np.int8(15))
    idx = np.argsort(rank, axis=1, kind="stable")[:, :13]   # (B,13)
    selP = np.take_along_axis(P, idx, axis=1)               # (B,13)
    selE = EXPS[idx]                                        # (B,13) int64

    sentinel = selP[:, 12] == 0                             # <13 nonzero
    keep = ~sentinel

    # I_{m0,n_j} = EXP[(log Pm0 * n_j - log Pn_j * m0) mod 255]
    logP = LOG32[selP]                                      # (B,13) int64
    m0   = selE[:, 0]                                       # (B,)
    lPm0 = logP[:, 0]                                       # (B,)
    nj   = selE[:, 1:13]                                    # (B,12)
    lPnj = logP[:, 1:13]                                    # (B,12)
    e    = (lPm0[:, None] * nj - lPnj * m0[:, None]) % 255  # (B,12) int64
    I    = EXP[e]                                           # (B,12) uint8
    I[sentinel] = 0

    I8 = I[keep, :8].astype(np.uint64)
    key = (I8[:, 0] << 56 | I8[:, 1] << 48 | I8[:, 2] << 40 | I8[:, 3] << 32 |
           I8[:, 4] << 24 | I8[:, 5] << 16 | I8[:, 6] <<  8 | I8[:, 7])
    return key, int(sentinel.sum())


# ------------------------------------------------------------- self-check

def self_check(rng, B=512):
    g = rng.integers(0, 256, size=(B, 255), dtype=np.uint8)
    # reference P via test_clumping.power_sums_fast (already validated)
    P_ref = power_sums_fast(g)
    # recompute via this module's path (reassemble P without packing)
    p_odd = np.zeros((B, 64), dtype=np.uint8)
    for w in range(255):
        p_odd ^= POW_ODD[g[:, w]]
    p = np.zeros((B, 128), dtype=np.uint8)
    p[:, ODD_R] = p_odd
    for r in range(1, 64):
        p[:, 2 * r] = SQ[p[:, r]]
    P = np.empty((B, 15), dtype=np.uint8)
    for j, m in enumerate(ALL_EXP):
        acc = p[:, m].copy()
        for (k, mk) in SUBMASK_PAIRS[m]:
            acc ^= MUL[p[:, k], p[:, mk]]
        P[:, j] = acc
    ok = np.array_equal(P_ref, P)
    print(f"[*] self-check vs validated power_sums_fast on {B}: "
          f"{'OK' if ok else 'FAIL'}")
    if not ok:
        sys.exit(1)


# ----------------------------------------------------- sort-based collisions

def count_prefix_collisions(keys_sorted, k):
    """keys_sorted: sorted (M,) uint64 big-endian. k: prefix bytes (1..8)."""
    s = 8 * (8 - k)
    v = keys_sorted >> np.uint64(s)
    eq = v[1:] == v[:-1]
    if not eq.any():
        return 0
    # run-lengths of True in eq -> value-run length = ell+1 -> C(ell+1,2)
    d = np.diff(np.concatenate(([0], eq.view(np.int8), [0])))
    starts = np.where(d == 1)[0]
    ends   = np.where(d == -1)[0]
    ell = (ends - starts).astype(np.int64)
    return int((ell * (ell + 1) // 2).sum())


# -------------------------------------------------------------------- main

def run(N=300_000_000, chunk=65536, seed=0, kmin=5, kmax=7):
    rng = np.random.default_rng(seed)
    self_check(rng)

    print(f"[*] allocating key buffer: {N*8/2**30:.2f} GiB ...")
    keys = np.empty(N, dtype=np.uint64)
    pos = 0
    sentinels = 0

    print(f"[*] streaming N={N:,} random 255-multisets, chunk={chunk} ...")
    t0 = time.time()
    done = 0
    while done < N:
        B = min(chunk, N - done)
        g = rng.integers(0, 256, size=(B, 255), dtype=np.uint8)
        k, s = fingerprint_chunk(g)
        keys[pos:pos + len(k)] = k
        pos += len(k)
        sentinels += s
        done += B
        if done % (chunk * 64) == 0 or done == N:
            rate = done / (time.time() - t0)
            eta = (N - done) / rate
            print(f"    {done:>13,}/{N:,}  "
                  f"{rate/1e6:.2f} M/s  ETA {eta/60:.1f} min", end="\r")
    print()
    M = pos
    keys = keys[:M]
    dt = time.time() - t0
    print(f"[*] generated {M:,} fingerprints in {dt:.1f}s "
          f"({M/dt/1e6:.2f} M samples/s); sentinels dropped: {sentinels}")

    print(f"[*] sorting {M:,} uint64 keys ...")
    t1 = time.time()
    keys.sort(kind="quicksort")
    print(f"    sort: {time.time()-t1:.1f}s")

    print(f"\n[*] PREFIX collision-entropy (k={kmin}..{kmax}):")
    print(f"    k   collisions C        H2_hat     ideal k*log2(255)")
    for k in range(kmin, kmax + 1):
        C = count_prefix_collisions(keys, k)
        if C > 0:
            H2 = np.log2(M) + np.log2(M - 1) - np.log2(2.0 * C)
        else:
            H2 = float("inf")
        ideal = k * np.log2(255)
        exp_C = M * (M - 1) / 2 / (255.0 ** k)
        print(f"    {k}   {C:>14,}   {H2:8.3f}    {ideal:7.3f}   "
              f"(expect ~{exp_C:.1f})")

    # full-8-byte: any exact duplicates?
    dup8 = count_prefix_collisions(keys, 8)
    print(f"\n[*] exact 8-byte duplicates: {dup8} "
          f"(expect ~{M*(M-1)/2/255.0**8:.2e})")
    print("[*] done.")


if __name__ == "__main__":
    N     = int(sys.argv[1]) if len(sys.argv) > 1 else 300_000_000
    chunk = int(sys.argv[2]) if len(sys.argv) > 2 else 65536
    run(N=N, chunk=chunk)
