#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Self-contained verification of the I_{m,n} power-sum invariant from report.tex §3.1.

Setup (paper notation, 7-round AES, rounds 0..6, whitening key k_{-1}):
  - Build a δ-set: 256 plaintexts whose round-1 states agree on bytes 1..15
    and whose x_1[0] ranges over all 256 values.
  - Encrypt each through 7 rounds.
  - "Offline" side knows the internal differences d_ω := x_5[0]_0 ⊕ x_5[0]_ω.
  - "Online" side peels the last round with the *correct* key bytes to get
    column 0 of x_6, forms v_ω := MC^{-1}(x_6)[0], then
        g_ω := ( L^{-1}(v_0 ⊕ v_ω) )^{-1}        (GF(256) inverse, 0^{-1}:=0)
  - Bridge identity (paper):  g_ω = s^2 · d_ω^{-1} ⊕ s,   s := x_5[0]_0.
  - Power sums:
        P_m = ⊕_{1≤ω<ω'≤255} (g_ω ⊕ g_ω')^m        (online)
        Q_m = ⊕_{1≤ω<ω'≤255} (d_ω^{-1} ⊕ d_ω'^{-1})^m   (offline)
        I_{m,n} = P_m^n / P_n^m   should equal   Q_m^n / Q_n^m.

This script checks, for a random key / random δ-set:
  (a) the elementwise bridge identity (away from the 0^{-1} corner case),
  (b) I_{m,n}(online) == I_{m,n}(offline) for m=7, n∈{11,13,19,23,29,31,37,43,47,53,59,61},
      using the paper's "add-0" parity patch so exactly one of the two
      offline fingerprints matches.
"""

import os

# ----------------------------------------------------------------------
# GF(2^8) with AES polynomial x^8+x^4+x^3+x+1  (0x11b)
# ----------------------------------------------------------------------
def _build_tables():
    exp = [0] * 512
    log = [0] * 256
    x = 1
    for i in range(255):
        exp[i] = x
        log[x] = i
        # multiply by generator 3
        x ^= (x << 1) ^ (0x11b if x & 0x80 else 0)
        x &= 0xFF
    for i in range(255, 512):
        exp[i] = exp[i - 255]
    return exp, log

EXP, LOG = _build_tables()

def gf_mul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]

def gf_inv(a):
    return 0 if a == 0 else EXP[255 - LOG[a]]

def gf_pow(a, e):
    if a == 0:
        return 0 if e != 0 else 1
    return EXP[(LOG[a] * e) % 255]

# ----------------------------------------------------------------------
# AES S-box  S(t) = L(t^{-1}) ⊕ 0x63,  and the GF(2)-linear map L / L^{-1}
# ----------------------------------------------------------------------
def _affine_L(b):
    r = 0
    for i in range(8):
        bit = ((b >> i) ^ (b >> ((i + 4) & 7)) ^ (b >> ((i + 5) & 7))
               ^ (b >> ((i + 6) & 7)) ^ (b >> ((i + 7) & 7))) & 1
        r |= bit << i
    return r

L_TBL  = [_affine_L(b) for b in range(256)]
LINV   = [0] * 256
for b in range(256):
    LINV[L_TBL[b]] = b

SBOX  = [L_TBL[gf_inv(b)] ^ 0x63 for b in range(256)]
SINV  = [0] * 256
for b in range(256):
    SINV[SBOX[b]] = b

assert SBOX[0x53] == 0xed and SBOX[0x00] == 0x63

# ----------------------------------------------------------------------
# AES round functions.  State is a flat list of 16 bytes, column-major:
#   index = row + 4*col,  so bytes 0..3 are column 0.
# ----------------------------------------------------------------------
def sub_bytes(s):     return [SBOX[b] for b in s]
def inv_sub_bytes(s): return [SINV[b] for b in s]

def shift_rows(s):
    out = [0] * 16
    for r in range(4):
        for c in range(4):
            out[r + 4 * c] = s[r + 4 * ((c + r) & 3)]
    return out

def inv_shift_rows(s):
    out = [0] * 16
    for r in range(4):
        for c in range(4):
            out[r + 4 * c] = s[r + 4 * ((c - r) & 3)]
    return out

def _mix_col(c, coeffs):
    a0, a1, a2, a3 = c
    m = coeffs
    return [
        gf_mul(m[0], a0) ^ gf_mul(m[1], a1) ^ gf_mul(m[2], a2) ^ gf_mul(m[3], a3),
        gf_mul(m[3], a0) ^ gf_mul(m[0], a1) ^ gf_mul(m[1], a2) ^ gf_mul(m[2], a3),
        gf_mul(m[2], a0) ^ gf_mul(m[3], a1) ^ gf_mul(m[0], a2) ^ gf_mul(m[1], a3),
        gf_mul(m[1], a0) ^ gf_mul(m[2], a1) ^ gf_mul(m[3], a2) ^ gf_mul(m[0], a3),
    ]

def mix_columns(s):
    out = [0] * 16
    for c in range(4):
        out[4*c:4*c+4] = _mix_col(s[4*c:4*c+4], (2, 3, 1, 1))
    return out

def inv_mix_columns(s):
    out = [0] * 16
    for c in range(4):
        out[4*c:4*c+4] = _mix_col(s[4*c:4*c+4], (14, 11, 13, 9))
    return out

def xor_state(a, b): return [x ^ y for x, y in zip(a, b)]

# ----------------------------------------------------------------------
# AES-128 key schedule → round keys rk[0..R]  (rk[0] is the paper's k_{-1})
# ----------------------------------------------------------------------
RCON = [0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36]

def key_schedule(master, n_rounds):
    rk = [list(master)]
    for i in range(n_rounds):
        prev = rk[-1]
        t = prev[12:16]
        t = [SBOX[t[1]] ^ RCON[i], SBOX[t[2]], SBOX[t[3]], SBOX[t[0]]]
        nk = []
        for c in range(4):
            col = [prev[4*c+r] ^ t[r] for r in range(4)]
            nk += col
            t = col
        rk.append(nk)
    return rk  # rk[0]=k_{-1}, rk[1]=k_0, ..., rk[n_rounds]=k_{n_rounds-1}

# ----------------------------------------------------------------------
# 7-round AES encrypt, returning every intermediate x_i (round inputs).
# Paper: x_0 = P ⊕ k_{-1};  rounds 0..6;  last round (6) has no MC.
# ----------------------------------------------------------------------
NROUNDS = 7

def encrypt_trace(pt, rk):
    xs = []
    x = xor_state(pt, rk[0])          # x_0
    xs.append(x)
    for r in range(NROUNDS - 1):      # rounds 0..5
        x = xor_state(mix_columns(shift_rows(sub_bytes(x))), rk[r + 1])
        xs.append(x)                  # x_{r+1}
    # last round, no MC
    ct = xor_state(shift_rows(sub_bytes(x)), rk[NROUNDS])
    return ct, xs                     # xs[i] = x_i for i=0..6

# ----------------------------------------------------------------------
# Build the δ-set: fix a random x_1 with x_1[0] varying over 0..255,
# invert round 0 to get the 256 plaintexts.
# ----------------------------------------------------------------------
def build_delta_set(rk):
    base_x1 = list(os.urandom(16))
    pts = []
    for i in range(256):
        x1 = list(base_x1)
        x1[0] = i
        # invert round 0:  x_1 = MC(SR(SB(x_0))) ⊕ k_0  ⇒
        x0 = inv_sub_bytes(inv_shift_rows(inv_mix_columns(xor_state(x1, rk[1]))))
        pt = xor_state(x0, rk[0])
        pts.append(pt)
    return pts

# ----------------------------------------------------------------------
# Online side: from ciphertexts + correct last-round key bytes,
# peel round 6 to get column 0 of x_6, then v = MC^{-1}(x_6)[0].
# ----------------------------------------------------------------------
def online_v(ct, k6):
    # need x_6[0..3].  y_6 = SR^{-1}(C ⊕ k6);  column 0 of y_6 uses
    # (C⊕k6) at flat indices 0, 13, 10, 7.
    z = [ct[j] ^ k6[j] for j in (0, 13, 10, 7)]
    x6_col0 = [SINV[b] for b in z]         # = x_6[0..3]
    a0, a1, a2, a3 = x6_col0
    return gf_mul(14, a0) ^ gf_mul(11, a1) ^ gf_mul(13, a2) ^ gf_mul(9, a3)

# ----------------------------------------------------------------------
# Power-sum invariant I_{m,n}
# ----------------------------------------------------------------------
def pairwise_power_sum(vals, m):
    """⊕_{ω<ω'} (vals[ω] ⊕ vals[ω'])^m   over indices 1..255."""
    acc = 0
    n = len(vals)
    for i in range(1, n):
        vi = vals[i]
        for j in range(i + 1, n):
            acc ^= gf_pow(vi ^ vals[j], m)
    return acc

def fingerprint(vals, m, ns):
    """Return tuple ( I_{m,n} for n in ns ).  vals indexed 0..; index 0 ignored."""
    exps = sorted(set(ns) | {m})
    P = {e: pairwise_power_sum(vals, e) for e in exps}
    out = []
    for n in ns:
        num = gf_pow(P[m], n)
        den = gf_pow(P[n], m)
        out.append(gf_mul(num, gf_inv(den)))
    return tuple(out)

# --- §"Degenerate P_m" fingerprint: first-13-nonzero along E ----------
E_LIST = (7, 11, 13, 19, 23, 29, 31, 37, 43, 47, 53, 59, 61, 91, 127)

# Fast P_m via report.tex eq. (lucas):
#   P_m = ⊕_{0<k<m-k, k⊆m} p_k p_{m-k}  ⊕  p_m      (m odd, |V| odd)
# where p_r = ⊕_ω g_ω^r.  We pack all p_r (r=0..127) into one 1024-bit int.
_POW_ROW = [int.from_bytes(bytes(gf_pow(g, r) for r in range(128)), 'little')
            for g in range(256)]

def _single_index_psums(vals):
    acc = 0
    for v in vals[1:]:
        acc ^= _POW_ROW[v]
    b = acc.to_bytes(128, 'little')
    return b  # p[r] = b[r]

def _Pm_from_p(p, m):
    # m odd; |V| parity = p[0]
    pm = p[m] if (p[0] & 1) else 0
    k = (m - 1) & m
    while k:
        j = m - k
        if k < j:
            pm ^= gf_mul(p[k], p[j])
        k = (k - 1) & m
    return pm

def all_Pm(vals, exps=E_LIST):
    p = _single_index_psums(vals)
    return {m: _Pm_from_p(p, m) for m in exps}

def fingerprint_nz(vals, _naive=False):
    """Paper's degenerate-P_m-safe fingerprint.
    Compute P_e for all e∈E; let m0 = first e with P_e≠0 and n_1..n_12 the
    next twelve nonzero ones; return (I_{m0,n_j})_j ∈ (GF(256)^×)^12.
    Returns None if fewer than 13 of the 15 P_e are nonzero (prob ~2^-15)."""
    P = ({e: pairwise_power_sum(vals, e) for e in E_LIST}
         if _naive else all_Pm(vals))
    nz = [e for e in E_LIST if P[e] != 0]
    if len(nz) < 13:
        return None
    m0, ns = nz[0], nz[1:13]
    return tuple(gf_mul(gf_pow(P[m0], n), gf_inv(gf_pow(P[n], m0))) for n in ns)

# ----------------------------------------------------------------------
# One trial
# ----------------------------------------------------------------------
def trial(seed=None, verbose=True):
    key = list(os.urandom(16)) if seed is None else list(seed)
    rk = key_schedule(key, NROUNDS)
    pts = build_delta_set(rk)

    cts  = []
    x5_0 = []
    for pt in pts:
        ct, xs = encrypt_trace(pt, rk)
        cts.append(ct)
        x5_0.append(xs[5][0])

    # --- offline side: differences d_ω and their inverses -------------
    s = x5_0[0]
    if s == 0:
        # Bridge identity degenerates (g_ω = a_ω, not s²d_ω⁻¹⊕s = 0).
        # Paper absorbs this as a 256/255 data factor; the attack misses
        # this δ-set.  Skip rather than assert.
        if verbose:
            print(f"key = {bytes(key).hex()}")
            print("s = 0  →  invariant inapplicable (paper: 256/255 data factor). SKIP")
        return None
    d     = [s ^ a for a in x5_0]            # d[0] = 0
    dinv  = [gf_inv(di) for di in d]

    # --- online side: v_ω, then g_ω ----------------------------------
    k6 = rk[NROUNDS]
    v  = [online_v(ct, k6) for ct in cts]
    g  = [gf_inv(LINV[v[0] ^ vw]) for vw in v]   # g[0] = 0

    # --- (a) bridge identity check -----------------------------------
    # bad indices: a_ω ∈ {0, s}  (where S-box 0^{-1}:=0 fiat breaks the algebra)
    bad = [w for w in range(1, 256) if x5_0[w] == 0 or x5_0[w] == s]
    s2 = gf_mul(s, s)
    mism = [w for w in range(1, 256)
            if w not in bad and g[w] != (gf_mul(s2, dinv[w]) ^ s)]
    if verbose:
        print(f"key = {bytes(key).hex()}")
        print(f"s = x_5[0]_0 = 0x{s:02x}   #bad (a_ω∈{{0,s}}, ω≥1) = {len(bad)}")
        print(f"bridge identity g_ω = s^2·d_ω^-1 ⊕ s  (good ω): "
              f"{'OK' if not mism else f'FAIL at {mism[:5]}'}")
    assert not mism, "bridge identity failed"

    m  = 7
    ns = [11, 13, 19, 23, 29, 31, 37, 43, 47, 53, 59, 61]

    # --- (b) CORE CLAIM: I_{m,n} matches on the good indices ----------
    # (both sides restricted to ω where the bridge identity is exact;
    #  this isolates the algebra of §3.1 from the 0^{-1} corner case)
    good = [w for w in range(1, 256) if w not in bad]
    g_good    = [0] + [g[w]    for w in good]   # leading 0 is the ignored ref slot
    dinv_good = [0] + [dinv[w] for w in good]
    I_on_core  = fingerprint(g_good,    m, ns)
    I_off_core = fingerprint(dinv_good, m, ns)
    core_ok = (I_on_core == I_off_core)
    if verbose:
        print(f"I_{{7,n}} (good ω only) online  = {['%02x'%b for b in I_on_core]}")
        print(f"I_{{7,n}} (good ω only) offline = {['%02x'%b for b in I_off_core]}")
        print(f"CORE invariant (good indices): {'OK' if core_ok else 'FAIL'}")
    assert core_ok, "I_{m,n} core invariant FAILED"

    # --- (c) full δ-set, attack-faithful ------------------------------
    # Online does NOT know which ω are bad.  Per the paper's 0^{-1} patch,
    # BOTH sides compute two fingerprints — "raw" and "add-0" (append one
    # synthetic zero: online g=0, offline d^{-1}=0) — and online accepts if
    # either pair matches.  Claim: exactly one matches, determined by the
    # parity of #bad = |{ω≥1 : a_ω ∈ {0,s}}|.
    I_on  = (fingerprint(g,          m, ns), fingerprint(g    + [0], m, ns))
    I_off = (fingerprint(dinv,       m, ns), fingerprint(dinv + [0], m, ns))
    hit = [p for p in (0, 1) if I_on[p] == I_off[p]]
    if verbose:
        print(f"full δ-set:  raw match={0 in hit}  add-0 match={1 in hit}  "
              f"#bad={len(bad)}")
    assert hit, "I_{m,n} invariant FAILED on full δ-set (neither parity)"
    assert hit == [len(bad) & 1], \
        f"parity prediction wrong: #bad={len(bad)} but hit={hit}"

    # --- (d) degenerate-P_m-safe fingerprint (first-13-nonzero) -------
    # Online/offline must select the same exponent subset (since
    # P_m = s^{2m} Q_m and s≠0 ⇒ zero patterns coincide).
    for parity in (0, 1):
        gp = g    + ([0] if parity else [])
        dp = dinv + ([0] if parity else [])
        zp_g = tuple(1 if pairwise_power_sum(gp, e) == 0 else 0 for e in E_LIST)
        zp_d = tuple(1 if pairwise_power_sum(dp, e) == 0 else 0 for e in E_LIST)
        if parity == (len(bad) & 1):
            assert zp_g == zp_d, f"zero-pattern mismatch at parity {parity}"
            fon, foff = fingerprint_nz(gp), fingerprint_nz(dp)
            assert fon is not None and fon == foff, "nz-fingerprint mismatch"
            assert 0 not in fon, "nz-fingerprint contains a zero byte"
    if verbose:
        print("first-13-nonzero fingerprint: zero-patterns agree, values match, all bytes ∈ GF(256)^×")
    return core_ok, hit[0], len(bad)

def online_fingerprints(cts, k6_guess):
    """Given 256 ciphertexts and a guess for k_6 at positions {0,7,10,13},
    compute v_ω, g_ω, and return (fp_raw, fp_add0) via fingerprint_nz.
    This is exactly what the online attacker evaluates per key guess;
    note it never touches u_5[0]."""
    v = [online_v(ct, k6_guess) for ct in cts]
    g = [gf_inv(LINV[v[0] ^ vw]) for vw in v]
    return fingerprint_nz(g), fingerprint_nz(g + [0])

_MCINV_COL0 = (14, 9, 13, 11)   # MC^{-1}[:,0]
DIAG = (0, 5, 10, 15)

def build_delta_set_from_guess(km1_guess, base_other):
    """Attacker-side δ-set construction.  Given a *guess* for k_{-1}[0,5,10,15],
    pick 256 plaintexts (from a diagonal structure) such that, under that guess,
    x_1 varies only in byte 0.  Concretely z_0[col0] = t·MC^{-1}[:,0] for
    t=0..255 ⇒ w_0[col0] varies only in row 0.  base_other supplies the 12
    fixed off-diagonal plaintext bytes."""
    pts = []
    for t in range(256):
        pt = list(base_other)
        for r, pos in enumerate(DIAG):
            pt[pos] = SINV[gf_mul(t, _MCINV_COL0[r])] ^ km1_guess[pos]
        pts.append(pt)
    return pts

def wrong_km1_test(verbose=True):
    """Sweep each byte of k_{-1}[0,5,10,15] (others correct, k_6 correct).
    Finding: the I_{m,n} invariant is *independent of k_{-1}* — the bridge
    identity g_ω = s²d_ω⁻¹⊕s relates online g (needs only correct k_6) to the
    true-internal d for *whatever* 256 plaintexts were encrypted.  So every
    k_{-1} guess matches its own true-internal offline point.  Filtering wrong
    k_{-1} comes from table *membership* (wrong guess ⇒ x_1 active in >1 byte
    ⇒ d-sequence not in the 11-byte-parametrized DFJ family), not from the
    algebra — that needs a (scaled-down) offline table, which is the next step."""
    import random
    key = list(os.urandom(16))
    rk  = key_schedule(key, NROUNDS)
    k6  = rk[NROUNDS]
    km1 = rk[0]
    base_other = list(os.urandom(16))
    if verbose:
        print(f"\n--- k_{{-1}} sweep (key={bytes(key).hex()}) ---")
        print(f"true k_{{-1}}[{DIAG}] = {[km1[p] for p in DIAG]}")

    def run(km1_guess):
        pts = build_delta_set_from_guess(km1_guess, base_other)
        cts, x5_0, x1 = [], [], []
        for pt in pts:
            ct, xs = encrypt_trace(pt, rk)
            cts.append(ct); x5_0.append(xs[5][0]); x1.append(xs[1])
        # how many x_1 bytes actually vary?
        varying = [j for j in range(16)
                   if any(x1[t][j] != x1[0][j] for t in range(1, 256))]
        s = x5_0[0]
        if s == 0:
            return None, varying
        dinv = [gf_inv(s ^ a) for a in x5_0]
        off  = (fingerprint_nz(dinv), fingerprint_nz(dinv + [0]))
        on   = online_fingerprints(cts, k6)
        return (on[0] == off[0] or on[1] == off[1]), varying

    m, vary = run(km1)
    assert m and vary == [0], f"correct k_{{-1}}: match={m}, x_1 varies in {vary}"
    if verbose:
        print(f"  correct guess: match=True, x_1 varies only in byte {vary}")

    for p in DIAG:
        nmatch = 0
        sample_vary = None
        for b in range(256):
            g = list(km1); g[p] = b
            m, vary = run(g)
            if m: nmatch += 1
            if b == (km1[p] ^ 1): sample_vary = vary
        if verbose:
            print(f"  sweep k_{{-1}}[{p:2d}]: online==true-internal offline for "
                  f"{nmatch}/256 guesses; wrong guess makes x_1 vary in {sample_vary}")
    if verbose:
        print("  ⇒ I_{m,n} algebra does NOT filter k_{-1}; table membership does.")

def wrong_key_test(n_random=2000, verbose=True):
    """For one random key/δ-set:
      - sweep each of the 4 guessed k_6 bytes over all 256 values (others
        correct) and count online↔offline fingerprint matches;
      - try n_random fully random 4-byte guesses.
    Expect exactly 1 match per sweep (the true byte) and 0 random matches."""
    import random
    key = list(os.urandom(16))
    rk  = key_schedule(key, NROUNDS)
    pts = build_delta_set(rk)
    cts, x5_0 = [], []
    for pt in pts:
        ct, xs = encrypt_trace(pt, rk)
        cts.append(ct); x5_0.append(xs[5][0])
    s = x5_0[0]
    if s == 0:
        return wrong_key_test(n_random, verbose)  # reroll (2^-8)
    dinv = [gf_inv(s ^ a) for a in x5_0]
    off  = (fingerprint_nz(dinv), fingerprint_nz(dinv + [0]))
    k6   = rk[NROUNDS]
    POS  = (0, 7, 10, 13)

    def hits(guess):
        on = online_fingerprints(cts, guess)
        return on[0] == off[0] or on[1] == off[1]

    # sanity: true key matches
    assert hits(k6), "correct k_6 guess failed to match!"
    if verbose:
        print(f"\n--- wrong-key test (key={bytes(key).hex()}) ---")
        print(f"correct k_6[{POS}] = {[k6[p] for p in POS]} → match ✓")

    for p in POS:
        matches = []
        for b in range(256):
            g = list(k6); g[p] = b
            if hits(g):
                matches.append(b)
        ok = matches == [k6[p]]
        if verbose:
            print(f"  sweep k_6[{p:2d}]: matches at {matches}  "
                  f"{'OK (only true byte)' if ok else 'UNEXPECTED'}")
        assert ok

    rnd_hits = 0
    for _ in range(n_random):
        g = list(k6)
        for p in POS:
            g[p] = random.randrange(256)
        if g[POS[0]]==k6[POS[0]] and g[POS[1]]==k6[POS[1]] and \
           g[POS[2]]==k6[POS[2]] and g[POS[3]]==k6[POS[3]]:
            continue
        rnd_hits += hits(g)
    if verbose:
        print(f"  {n_random} random wrong 4-byte guesses: {rnd_hits} matches "
              f"(expect 0; each ~2^-96 vs this one offline point)")
    assert rnd_hits == 0

def distribution_test(n_samples=20000):
    """Empirical uniformity of fingerprint_nz over random 255-element inputs
    (models wrong-key online g-lists).  Reports per-byte χ² over GF(256)^×
    and pairwise collision rate vs the ideal 255^{-12}."""
    import random, math
    hist = [[0] * 256 for _ in range(12)]
    fps = []
    skipped = 0
    for _ in range(n_samples):
        vals = [0] + [random.randrange(256) for _ in range(255)]
        fp = fingerprint_nz(vals)
        if fp is None:
            skipped += 1
            continue
        fps.append(fp)
        for j, b in enumerate(fp):
            hist[j][b] += 1
    N = len(fps)
    print(f"\n--- distribution test: {N} fingerprints ({skipped} skipped, <13 nonzero P_e) ---")
    exp = N / 255.0
    for j in range(12):
        h = hist[j]
        assert h[0] == 0, f"byte {j}: saw value 0 (should be impossible)"
        chi2 = sum((h[v] - exp) ** 2 / exp for v in range(1, 256))
        lo, hi = min(h[v] for v in range(1, 256)), max(h[v] for v in range(1, 256))
        print(f"  byte[{j:2d}]  χ²(254 dof)={chi2:7.1f}   min/exp/max = {lo}/{exp:.1f}/{hi}")
    # pairwise collision rate
    from collections import Counter
    cnt = Counter(fps)
    coll = sum(c * (c - 1) // 2 for c in cnt.values())
    pairs = N * (N - 1) // 2
    ideal = pairs * 255.0 ** -12
    print(f"  full-tuple collisions: {coll}  (ideal ≈ {ideal:.2e} over {pairs:.2e} pairs)")
    # effective per-byte collision prob (treat bytes as independent)
    pcoll = [sum(c * (c - 1) // 2 for c in Counter(fp[j] for fp in fps).values()) / pairs
             for j in range(12)]
    print(f"  per-byte collision prob: "
          f"min={min(pcoll):.5f} max={max(pcoll):.5f} ideal=1/255={1/255:.5f}")
    bits = -sum(math.log2(p) for p in pcoll)
    print(f"  naive independent-byte entropy estimate ≈ {bits:.1f} bits (structural max ≈ 95.9)")

if __name__ == "__main__":
    import sys
    ntrials = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    for t in range(ntrials):
        print(f"\n=== trial {t} ===")
        trial()
    print(f"\nAll {ntrials} trials: bridge identity OK, core I_{{m,n}} OK, "
          f"full-δ-set add-0 patch OK, first-13-nonzero fingerprint OK.")
    # cross-check fast P_m against naive pairwise sum once
    import random
    tv = [0] + [random.randrange(256) for _ in range(255)]
    assert fingerprint_nz(tv) == fingerprint_nz(tv, _naive=True), "fast P_m disagrees with naive"
    assert fingerprint_nz(tv + [0]) == fingerprint_nz(tv + [0], _naive=True)
    print("fast P_m (Lucas reduction) matches naive pairwise sum ✓")
    if len(sys.argv) > 2 and sys.argv[2] == "wrong":
        for _ in range(int(sys.argv[3]) if len(sys.argv) > 3 else 3):
            wrong_key_test()
    if len(sys.argv) > 2 and sys.argv[2] == "km1":
        wrong_km1_test()
    if len(sys.argv) > 2 and sys.argv[2] == "dist":
        distribution_test(int(sys.argv[3]) if len(sys.argv) > 3 else 20000)
