#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Wrong-key randomization test for the I_{m,n} fingerprint.

For a random AES-128 key and its correct δ-set at x_1:
  - compute the correct offline fingerprint I_off (both parities);
  - for many WRONG guesses of k_6[0,7,10,13], compute the online fingerprint
    I_on and check:
      (a) ghost keys   : any wrong guess that matches I_off?
      (b) near-miss    : for 1-byte-off k_6, is I_on far from I_off?
      (c) clustering   : do wrong-key I_on collide with each other more than
                         birthday-expected?

Reuses the AES + bridge machinery from test_bridge_aes.py.
"""

import sys
import numpy as np
from collections import Counter

from test_invariants import gf_mul, gf_inv
from test_bridge_aes import (SBOX, SBOX_INV, A_LIN_INV, FIELD_INV,
                             key_schedule, encrypt_from_x1)
from test_zeroinv import fingerprint_fb_single

K6_POS = np.array([0, 13, 10, 7], dtype=np.intp)   # z_6 positions s.t. y_6[0..3]=z_6[K6_POS]


def online_fingerprints(C, k6_guesses):
    """
    C: (256,16) ciphertexts of the δ-set.
    k6_guesses: (G,4) uint8, each a guess for k_6[K6_POS].
    Returns (G, 2, 12) uint8: for each guess, [raw, append-0] fingerprints.
    """
    G = k6_guesses.shape[0]
    # z_6[K6_POS] = C[K6_POS] ⊕ guess  -> broadcast over 256 δ-set elements
    Cpos = C[:, K6_POS]                            # (256,4)
    out = np.empty((G, 2, 12), dtype=np.uint8)
    for gi in range(G):
        z6 = Cpos ^ k6_guesses[gi]                 # (256,4)
        x6 = SBOX_INV[z6]                          # (256,4) = x_6[0,1,2,3]
        v  = (gf_mul(14, x6[:, 0]) ^ gf_mul(11, x6[:, 1]) ^
              gf_mul(13, x6[:, 2]) ^ gf_mul( 9, x6[:, 3]))
        dv = v[0] ^ v[1:]
        g  = FIELD_INV[A_LIN_INV[dv]]              # (255,)
        out[gi, 0] = fingerprint_fb_single(g)
        out[gi, 1] = fingerprint_fb_single(np.concatenate([g, [np.uint8(0)]]))
    return out


def run(N_keys=5, N_random_wrong=100_000, seed=0):
    rng = np.random.default_rng(seed)
    print(f"[*] wrong-key randomization test: {N_keys} keys, "
          f"{N_random_wrong:,} random + 1020 near-miss k_6 guesses each\n")

    for trial in range(N_keys):
        key = rng.integers(0, 256, size=16, dtype=np.uint8)
        rk  = key_schedule(key)
        base = rng.integers(0, 256, size=16, dtype=np.uint8)
        x1 = np.tile(base, (256, 1)); x1[:, 0] = np.arange(256, dtype=np.uint8)
        C, x5 = encrypt_from_x1(x1, rk)
        a = x5[:, 0]
        if a[0] == 0:
            print(f"  key {trial}: s=0, skip")
            continue
        # offline
        d = a[0] ^ a[1:]
        dinv = FIELD_INV[d]
        I_off = np.stack([
            fingerprint_fb_single(dinv),
            fingerprint_fb_single(np.concatenate([dinv, [np.uint8(0)]])),
        ])                                          # (2,12)

        k6_true = rk[7][K6_POS]

        # --- (a) correct guess sanity ------------------------------------
        I_right = online_fingerprints(C, k6_true[None, :])[0]  # (2,12)
        right_match = any(np.array_equal(I_right[p], I_off[p]) for p in (0,1))
        assert right_match, "correct k_6 failed to match — harness bug"

        # --- near-miss guesses: flip one of the 4 bytes to each wrong value
        near = []
        for pos in range(4):
            for val in range(256):
                if val == k6_true[pos]:
                    continue
                g = k6_true.copy(); g[pos] = val
                near.append(g)
        near = np.array(near, dtype=np.uint8)       # (1020,4)

        # --- random wrong guesses (reject the true one if drawn) ----------
        rnd = rng.integers(0, 256, size=(N_random_wrong, 4), dtype=np.uint8)
        same = np.all(rnd == k6_true, axis=1)
        rnd[same] ^= 1                              # nudge away (negligible)

        all_wrong = np.concatenate([near, rnd], axis=0)
        I_on = online_fingerprints(C, all_wrong)    # (W,2,12)
        W = I_on.shape[0]

        # (a) ghost keys: any wrong guess matching I_off under either parity?
        ghost = 0
        for p in (0, 1):
            ghost += int(np.all(I_on[:, p, :] == I_off[p], axis=1).sum())
        # (b) Hamming distance (bytes) of near-miss raw fp to correct raw fp
        #     (compare under whichever parity the correct key matched)
        p_ok = 0 if np.array_equal(I_right[0], I_off[0]) else 1
        hd_near = (I_on[:1020, p_ok, :] != I_off[p_ok]).sum(axis=1)  # (1020,)
        # (c) clustering among wrong-key fingerprints (use raw parity)
        fps = [bytes(r) for r in I_on[:, 0, :]]
        cnt = Counter(fps)
        coll = sum(v*(v-1)//2 for v in cnt.values())
        exp_coll = W*(W-1)/2 / 255.0**12

        print(f"  key {trial}: W={W:,} wrong guesses")
        print(f"    ghost matches to I_off         : {ghost}"
              f"{'   <-- !!' if ghost else ''}")
        print(f"    near-miss Hamming dist to I_off: "
              f"min={hd_near.min()}, mean={hd_near.mean():.2f}/12, "
              f"max={hd_near.max()}   (uniform-random expects ~11.95)")
        print(f"    wrong-key self-collisions      : {coll} "
              f"(birthday @95.9b expects {exp_coll:.2e})")
        # prefix-H2 on wrong-key fps, k=1,2
        for k in (1, 2):
            c = Counter(bytes(r[:k]) for r in I_on[:, 0, :])
            Ck = sum(v*(v-1)//2 for v in c.values())
            H2 = np.log2(W*(W-1)/(2*Ck)) if Ck else float('inf')
            print(f"    wrong-key H2(k={k})             : {H2:.3f} "
                  f"(ideal {k*np.log2(255):.3f})")
        print()

    print("[*] done.")


if __name__ == "__main__":
    N_keys = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    N_rand = int(sys.argv[2]) if len(sys.argv) > 2 else 100_000
    run(N_keys=N_keys, N_random_wrong=N_rand)
