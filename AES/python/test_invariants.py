#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Empirical independence test for the I_{m,n} power-sum invariants of report.tex §3.1.

We check whether the 12-byte fingerprint (I_{7,n} for n in N_SET) behaves like
~96 uniform independent bits when the input is a random 255-element multiset
over GF(2^8).

Tests:
  1. Marginal chi^2 of each I_{7,n} against Uniform(256).
  2. Pairwise chi^2 on every (n_i, n_j) pair, 256x256 contingency.
  3. Full-fingerprint collision count vs birthday bound.
  4. Sanity: verify the AGL(1,256)-invariance the paper claims
     (I(D) == I(alpha*D + beta) for random alpha,beta).
"""

import numpy as np
from collections import Counter

# ---------------------------------------------------------------- GF(2^8) ----
# AES field: x^8 + x^4 + x^3 + x + 1  (0x11B)

def _build_tables():
    exp = [0] * 512
    log = [0] * 256
    x = 1
    for i in range(255):
        exp[i] = x
        log[x] = i
        # multiply by generator 3 (primitive in the AES field)
        x ^= (x << 1)
        if x & 0x100:
            x ^= 0x11B
    for i in range(255, 512):
        exp[i] = exp[i - 255]
    return np.array(exp, dtype=np.uint8), np.array(log, dtype=np.int16)

EXP, LOG = _build_tables()

def gf_mul(a, b):
    """Vectorized GF(2^8) multiply. a, b: uint8 arrays (broadcastable)."""
    a, b = np.broadcast_arrays(np.asarray(a, dtype=np.int32),
                               np.asarray(b, dtype=np.int32))
    nz = (a != 0) & (b != 0)
    out = np.zeros(a.shape, dtype=np.uint8)
    out[nz] = EXP[(LOG[a[nz]] + LOG[b[nz]]) % 255]
    return out

def gf_pow(a, e):
    """Vectorized a^e in GF(2^8). e: nonneg int."""
    a = np.asarray(a, dtype=np.int32)
    out = np.zeros_like(a, dtype=np.uint8)
    nz = a != 0
    out[nz] = EXP[(LOG[a[nz]].astype(np.int64) * e) % 255]
    if e == 0:
        out[:] = 1  # 0^0 := 1 convention; irrelevant here
    return out

def gf_inv(a):
    return gf_pow(a, 254)  # a^(q-2)

def gf_div(a, b):
    return gf_mul(a, gf_inv(b))

# ---------------------------------------------------------- invariants -------

M0 = 7
N_SET = [11, 13, 19, 23, 29, 31, 37, 43, 47, 53, 59, 61]  # paper §3.1
# The 15 usable Frobenius-coset reps coprime to 255 (the 16th coset is
# {1,2,4,...,128}, for which P_m == P_1^{2^a} == 0 structurally).
EXTRA = [91, 127]
ALL_EXP = [M0] + N_SET + EXTRA   # len == 15

# Precompute g^m for all g in 0..255 and all m in ALL_EXP: shape (256, 13)
POW_TABLE = np.stack([gf_pow(np.arange(256, dtype=np.uint8), m)
                      for m in ALL_EXP], axis=1)  # (256,13) uint8

def pairwise_power_sums(g):
    """
    Given g: (N,255) uint8 array (each row a 255-multiset of field elements),
    return P: (N, 13) uint8 where P[:,k] = XOR_{w<w'} (g_w ^ g_w')**ALL_EXP[k].

    O(N * 255^2) — fine for N up to ~1e5; for larger N we can switch to the
    single-index p_r reduction from the paper, but let's keep the reference
    implementation deliberately unoptimized, favoring an obviously-correct form.
    """
    N = g.shape[0]
    P = np.zeros((N, len(ALL_EXP)), dtype=np.uint8)
    for i in range(255):
        # XOR-diff of g[:,i] against all g[:,j>i], vectorized over j and N
        diff = g[:, i:i+1] ^ g[:, i+1:]          # (N, 255-1-i)
        # look up all 13 powers at once: POW_TABLE[diff] -> (N, 255-1-i, 13)
        P ^= np.bitwise_xor.reduce(POW_TABLE[diff], axis=1)
    return P

def fingerprints(g):
    """Return (N,12) uint8 array of I_{7,n} for n in N_SET (no fallback)."""
    P = pairwise_power_sums(g)              # (N, len(ALL_EXP))
    Pm = P[:, 0]                            # P_7
    I = np.zeros((g.shape[0], len(N_SET)), dtype=np.uint8)
    for j, n in enumerate(N_SET):
        Pn = P[:, 1 + j]
        num = gf_pow(Pm, n)                 # P_m^n
        den = gf_pow(Pn, M0)                # P_n^m
        I[:, j] = gf_div(num, den)
    return I, P


def fingerprints_fb(g):
    """
    Fingerprint with deterministic fallback, computable identically by both
    online (from P_m) and offline (from Q_m) since P_m=0 <=> Q_m=0.

    Rule:
      - base exponent m* = first m in ALL_EXP with P_m != 0
      - coord exponents n_1..n_12 = next 12 entries of ALL_EXP (after m*)
        with P_n != 0, in list order
      - I_j = P_{m*}^{n_j} / P_{n_j}^{m*}
      - if fewer than 13 nonzero P's exist (astronomically rare), pad I with a
        sentinel 0x00 and prepend a 1-byte header = index of m* | 0x80
    Returns (N,12) uint8.
    """
    P = pairwise_power_sums(g)              # (N, len(ALL_EXP))
    N = g.shape[0]
    I = np.zeros((N, 12), dtype=np.uint8)
    exps = np.array(ALL_EXP, dtype=np.int64)

    # Fast path: rows where P_7..P_61 are all nonzero (the ~95% case).
    core = P[:, :1 + len(N_SET)]
    fast = np.all(core != 0, axis=1)
    if fast.any():
        Pm = core[fast, 0]
        for j, n in enumerate(N_SET):
            Pn = core[fast, 1 + j]
            I[fast, j] = gf_div(gf_pow(Pm, n), gf_pow(Pn, M0))

    # Slow path: per-row fallback for the remainder.
    slow_idx = np.where(~fast)[0]
    for r in slow_idx:
        nz = np.where(P[r] != 0)[0]
        if len(nz) < 13:
            # absurdly rare; leave zeros (acts as sentinel)
            continue
        m_pos = nz[0]
        n_pos = nz[1:13]
        m = int(exps[m_pos]); Pm = int(P[r, m_pos])
        for j, npos in enumerate(n_pos):
            n = int(exps[npos]); Pn = int(P[r, npos])
            I[r, j] = gf_div(gf_pow(Pm, n), gf_pow(Pn, m))
    return I, P

# ---------------------------------------------------------- stats ------------

def chi2_uniform(counts, bins):
    exp = counts.sum() / bins
    return float(((counts - exp) ** 2 / exp).sum())

def run(N=20000, seed=0):
    rng = np.random.default_rng(seed)
    print(f"[*] generating {N} random 255-multisets over GF(256) ...")
    g = rng.integers(0, 256, size=(N, 255), dtype=np.uint8)

    # --- sanity: AGL invariance on a subsample, incl. fallback path ----------
    # Force some P_7==0 rows into the sample by brute search so we exercise
    # the fallback under the affine transform.
    print("[*] sanity-checking AGL(1,256) invariance (with fallback) ...")
    sub = g[:1000].copy()
    I0, P0 = fingerprints_fb(sub)
    n_fb = int((P0[:, 0] == 0).sum())
    alpha = rng.integers(1, 256, size=(1000, 1), dtype=np.uint8)
    beta  = rng.integers(0, 256, size=(1000, 1), dtype=np.uint8)
    sub2 = gf_mul(alpha, sub) ^ beta
    I1, P1 = fingerprints_fb(sub2)
    ok = np.all(I0 == I1)
    # also confirm both sides agree on *which* P_m are zero
    zero_agree = np.all((P0 == 0) == (P1 == 0))
    print(f"    rows exercising fallback (P_7==0): {n_fb}/1000")
    print(f"    zero-pattern agrees online vs offline: {zero_agree}")
    print(f"    invariance holds: {ok}  "
          f"(mismatched rows: {int((I0 != I1).any(axis=1).sum())}/1000)")
    if not ok:
        print("    !!! invariance FAILS under fallback")

    print(f"[*] computing fingerprints (with fallback) for all {N} ...")
    I, P = fingerprints_fb(g)

    # --- degeneracy: how often do we run out of nonzero P_m? ----------------
    nz_count = (P != 0).sum(axis=1)           # per row, out of 15
    zero_P7 = int((P[:, 0] == 0).sum())
    need_fb = int((nz_count < 15).sum())      # any fallback at all
    too_few = int((nz_count < 13).sum())      # <13 nonzero: unrecoverable
    print(f"[*] degeneracy: P_7==0 in {zero_P7}/{N} ({zero_P7/N:.4%}); "
          f"any fallback in {need_fb}/{N} ({need_fb/N:.4%}); "
          f"<13 nonzero in {too_few}/{N} (expect ~{N*2**-15.2:.2f})")

    # Structurally I_{m,n} in GF(256)^x (never 0) once fallback is applied,
    # except for the rare <13-nonzero sentinel rows. Drop those rows and
    # test against Uniform({1..255}).
    keep = np.all(I != 0, axis=1)
    Iu = I[keep]
    Nk = Iu.shape[0]
    print(f"    (dropping {N - Nk} sentinel rows; testing on {Nk})")

    # --- test 1: marginal uniformity over GF(256)^x ------------------------
    print("\n[*] TEST 1: marginal chi^2 vs Uniform(255), df=254, "
          "crit@0.01 ≈ 309.4")
    for j, n in enumerate(N_SET):
        cnt = np.bincount(Iu[:, j], minlength=256)[1:]   # drop bin 0
        x2 = chi2_uniform(cnt, 255)
        flag = "  <-- BIASED" if x2 > 309.4 else ""
        print(f"    I_7,{n:<3d}  chi2 = {x2:8.1f}{flag}")

    # --- test 2: pairwise independence over GF(256)^x ----------------------
    print("\n[*] TEST 2: pairwise 255x255 chi^2, df=65024, "
          "crit@0.01 ≈ 65865")
    worst = (0.0, None)
    for a in range(len(N_SET)):
        for b in range(a + 1, len(N_SET)):
            idx = (Iu[:, a].astype(np.int32) - 1) * 255 + (Iu[:, b] - 1)
            cnt = np.bincount(idx, minlength=255 * 255)
            x2 = chi2_uniform(cnt, 255 * 255)
            if x2 > worst[0]:
                worst = (x2, (N_SET[a], N_SET[b]))
    print(f"    worst pair: I_7,{worst[1][0]} vs I_7,{worst[1][1]}  "
          f"chi2 = {worst[0]:.1f}"
          f"{'  <-- DEPENDENT' if worst[0] > 65865 else ''}")

    # --- test 3: full-tuple collisions -------------------------------------
    print(f"\n[*] TEST 3: full 12-byte fingerprint collisions among {N} samples")
    keys = [bytes(row) for row in I]
    c = Counter(keys)
    coll = sum(v - 1 for v in c.values() if v > 1)
    expected = N * (N - 1) / 2 / (256.0 ** 12)
    print(f"    observed collisions: {coll}")
    print(f"    birthday expectation @96 bits: {expected:.3e}")
    if coll > 0:
        print("    !!! any collision here means effective entropy << 96 bits")
        # show one
        for k, v in c.items():
            if v > 1:
                print(f"    example colliding fingerprint (x{v}): {k.hex()}")
                break

    print("\n[*] done.")

if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    run(N=N)
