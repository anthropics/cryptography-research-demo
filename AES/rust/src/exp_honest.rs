// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! E4: the honest bridge on real 7-round AES-128 (planted key).
use crate::aes::*;
use crate::bridge::*;
use crate::gf::Gf;
use crate::rng::Rng;
use rayon::prelude::*;

#[derive(Default, Clone)]
pub struct Stats {
    pub n: u64,
    pub s_zero: u64,
    pub b_hist: Vec<u64>,           // histogram of b = A + B
    pub raw_pm_match: u64,          // P_m == s^{2m} Q_m for all 15 m (raw parity)
    pub app_pm_match: u64,          // appended-parity variant matches
    pub both_pm_match: u64,
    pub neither_pm_match: u64,
    pub raw_match_when_b_even: u64, pub b_even: u64,
    pub app_match_when_b_odd: u64, pub b_odd: u64,
    pub imn_defined: u64,           // >= 13 nonzero P_m online (raw)
    pub imn_key_match_raw: u64,     // I_on == I_off (raw)
    pub imn_key_match_app: u64,     // appended
    pub imn_key_match_either: u64,
    pub imn_zero_pattern_agree: u64, // online/offline see the same P_m zero pattern (raw, matching parity)
    pub imn_zero_pattern_checked: u64,
    pub alg1_raw_pass: u64,         // exact multiset {g} == {s^2 d^-1 + s}
    pub alg1_app_pass: u64,         // with the appended zero
    pub alg1_either_pass: u64,
    pub corrected_pass: u64,        // multiset equality after deleting values {0, s}
    pub s_lit_correct: u64,         // s = (P_7/Q_7)^{1/14} (paper's step 2, raw)
    pub s_corr_correct: u64,        // parity-corrected version
    pub s_undef: u64,               // Q_7 = 0 (or P_7 = 0) so step 2 undefined
    pub keyless_multiset_ok: u64,   // {a ^ a_w} == {d(e)} from the 24 state bytes
    pub peel_ok: u64,               // ciphertext peel reproduces x_5[0]
    pub chi_raw_match: u64, pub chi_app_match: u64, pub chi_either_match: u64, pub chi_undef: u64,
    pub app_needed: u64,            // b odd and s != 0: the appended variant is the one that must match
    pub app_needed_undef: u64,      // ... but the appended online fingerprint is undefined (<13 nonzero P_m)
    pub app_undef_s1zero: u64,      // ... and S_1 (= p_1) of the online multiset is 0
    pub s1_zero_online: u64,        // S_1 = 0 for the online multiset (independent of b)
    pub raw_undef: u64,             // raw online fingerprint undefined
    pub s_corr_defined: u64,
}

impl Stats {
    fn new() -> Stats { let mut s = Stats::default(); s.b_hist = vec![0; 40]; s }
    fn merge(mut a: Stats, b: Stats) -> Stats {
        a.n += b.n; a.s_zero += b.s_zero;
        for i in 0..40 { a.b_hist[i] += b.b_hist[i]; }
        a.raw_pm_match += b.raw_pm_match; a.app_pm_match += b.app_pm_match;
        a.both_pm_match += b.both_pm_match; a.neither_pm_match += b.neither_pm_match;
        a.raw_match_when_b_even += b.raw_match_when_b_even; a.b_even += b.b_even;
        a.app_match_when_b_odd += b.app_match_when_b_odd; a.b_odd += b.b_odd;
        a.imn_defined += b.imn_defined;
        a.imn_key_match_raw += b.imn_key_match_raw; a.imn_key_match_app += b.imn_key_match_app;
        a.imn_key_match_either += b.imn_key_match_either;
        a.imn_zero_pattern_agree += b.imn_zero_pattern_agree; a.imn_zero_pattern_checked += b.imn_zero_pattern_checked;
        a.alg1_raw_pass += b.alg1_raw_pass; a.alg1_app_pass += b.alg1_app_pass;
        a.alg1_either_pass += b.alg1_either_pass; a.corrected_pass += b.corrected_pass;
        a.s_lit_correct += b.s_lit_correct; a.s_corr_correct += b.s_corr_correct; a.s_undef += b.s_undef;
        a.keyless_multiset_ok += b.keyless_multiset_ok; a.peel_ok += b.peel_ok;
        a.chi_raw_match += b.chi_raw_match; a.chi_app_match += b.chi_app_match;
        a.chi_either_match += b.chi_either_match; a.chi_undef += b.chi_undef;
        a.app_needed += b.app_needed; a.app_needed_undef += b.app_needed_undef; a.app_undef_s1zero += b.app_undef_s1zero;
        a.s1_zero_online += b.s1_zero_online; a.raw_undef += b.raw_undef; a.s_corr_defined += b.s_corr_defined;
        a
    }
}

/// Keyless offline computation of d(e) = Delta x_5[0] as a function of the
/// 24 anchor bytes (x_2^(0)[0..3], x_3^(0)[0..15], x_4^(0)[0,5,10,15]) and
/// the round-1 output difference e = Delta z_1[0].
pub fn keyless_d(gf: &Gf, x2c0: &[u8; 4], x3: &State, x4d: &[u8; 4], e: u8) -> u8 {
    // Delta x_2 column 0 = MC*(e,0,0,0) = (2e, e, e, 3e)
    let dx2 = [gf.mul2[e as usize], e, e, gf.mul3[e as usize]];
    let mut dy2 = [0u8; 16];
    for j in 0..4 { dy2[j] = gf.sbox[x2c0[j] as usize] ^ gf.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
    let dz2 = shift_rows(&dy2);
    let dx3 = mix_columns(gf, &dz2);
    let mut dy3 = [0u8; 16];
    for j in 0..16 { dy3[j] = gf.sbox[x3[j] as usize] ^ gf.sbox[(x3[j] ^ dx3[j]) as usize]; }
    let dz3 = shift_rows(&dy3);
    let dx4 = mix_columns(gf, &dz3);
    // x_5[0] = 2*y_4[0] + 3*y_4[5] + y_4[10] + y_4[15] (+ key)
    let diag = [0usize, 5, 10, 15];
    let coef = [2u8, 3, 1, 1];
    let mut d = 0u8;
    for r in 0..4 {
        let j = diag[r];
        let dy = gf.sbox[x4d[r] as usize] ^ gf.sbox[(x4d[r] ^ dx4[j]) as usize];
        d ^= gf.mul(coef[r], dy);
    }
    d
}

fn one_instance(gf: &Gf, pt: &PowTab, rng: &mut Rng, st: &mut Stats) {
    let key = rng.state16();
    let rk = expand_key(gf, &key);
    let base_x1 = rng.state16();
    let mut a = [0u8; 256];       // a_v = x_5[0]^(v)
    let mut w = [0u8; 256];       // v_v = MC^{-1}(x_6)[0] from ciphertext peel
    let pos = [0usize, 13, 10, 7];
    let k6 = rk[7];
    let u5 = inv_mix_columns(gf, &rk[6]);
    let mut peel_ok = true;
    let mut anchor_x2 = [0u8; 4]; let mut anchor_x3 = [0u8; 16]; let mut anchor_x4 = [0u8; 4];
    for v in 0..256usize {
        let mut x1 = base_x1; x1[0] = base_x1[0] ^ (v as u8);
        let p = plaintext_from_x1(gf, &rk, &x1);
        let t = encrypt_trace(gf, &rk, &p, 7);
        a[v] = t.xs[5][0];
        if v == 0 {
            for j in 0..4 { anchor_x2[j] = t.xs[2][j]; }
            anchor_x3 = t.xs[3];
            let dg = [0usize, 5, 10, 15];
            for j in 0..4 { anchor_x4[j] = t.xs[4][dg[j]]; }
        }
        // peel round 6 with true k6 at {0,13,10,7}
        let mut x6c0 = [0u8; 4];
        for j in 0..4 { x6c0[j] = gf.inv_sbox[(t.c[pos[j]] ^ k6[pos[j]]) as usize]; }
        w[v] = gf.mul14[x6c0[0] as usize] ^ gf.mul11[x6c0[1] as usize] ^ gf.mul13[x6c0[2] as usize] ^ gf.mul9[x6c0[3] as usize];
        if gf.inv_sbox[(w[v] ^ u5[0]) as usize] != a[v] { peel_ok = false; }
    }
    st.n += 1;
    if peel_ok { st.peel_ok += 1; }
    let s = a[0];
    if s == 0 { st.s_zero += 1; }
    // online / offline value multisets
    let mut g = [0u8; 255]; let mut dinv = [0u8; 255];
    let mut cnt_a = 0u32; let mut cnt_b = 0u32;
    for v in 1..256usize {
        g[v - 1] = gf.inv[gf.linv[(w[0] ^ w[v]) as usize] as usize];
        dinv[v - 1] = gf.inv[(a[0] ^ a[v]) as usize];
        if a[v] == 0 { cnt_a += 1; }
        if a[v] == a[0] { cnt_b += 1; }
    }
    let b = cnt_a + cnt_b;
    st.b_hist[(b as usize).min(39)] += 1;
    // predicted online values from the offline formula
    let s2 = gf.mul(s, s);
    let mut pred = [0u8; 255];
    for i in 0..255 { pred[i] = gf.mul(s2, dinv[i]) ^ s; }
    // Alg 1 exact multiset test (raw)
    let mut gs = g; gs.sort_unstable();
    let mut ps = pred; ps.sort_unstable();
    let raw_pass = gs == ps;
    if raw_pass { st.alg1_raw_pass += 1; }
    // appended: online adds value 0, offline formula on appended d^-1 = 0 gives value s
    let mut ga: Vec<u8> = g.to_vec(); ga.push(0); ga.sort_unstable();
    let mut pa: Vec<u8> = pred.to_vec(); pa.push(s); pa.sort_unstable();
    let app_pass = ga == pa;
    if app_pass { st.alg1_app_pass += 1; }
    if raw_pass || app_pass { st.alg1_either_pass += 1; }
    // corrected test: delete all values in {0, s} from both, then compare
    let gc: Vec<u8> = g.iter().cloned().filter(|&x| x != 0 && x != s).collect();
    let pc: Vec<u8> = pred.iter().cloned().filter(|&x| x != 0 && x != s).collect();
    let mut gc2 = gc; gc2.sort_unstable(); let mut pc2 = pc; pc2.sort_unstable();
    if gc2 == pc2 { st.corrected_pass += 1; }
    // power sums / P_m
    let p_on = power_sums(gf, pt, &g, 0);
    let p_off = power_sums(gf, pt, &dinv, 0);
    let mut p_on_app = p_on; p_on_app[0] ^= 1;
    let mut p_off_app = p_off; p_off_app[0] ^= 1;
    let pms_on = all_pms(gf, &p_on); let pms_off = all_pms(gf, &p_off);
    let pms_on_app = all_pms(gf, &p_on_app); let pms_off_app = all_pms(gf, &p_off_app);
    let mut raw_ok = true; let mut app_ok = true;
    for i in 0..15 {
        let m = E15[i] as u32;
        let s2m = gf.pow(s, 2 * m);
        if pms_on[i] != gf.mul(s2m, pms_off[i]) { raw_ok = false; }
        if pms_on_app[i] != gf.mul(s2m, pms_off_app[i]) { app_ok = false; }
    }
    if raw_ok { st.raw_pm_match += 1; }
    if app_ok { st.app_pm_match += 1; }
    if raw_ok && app_ok { st.both_pm_match += 1; }
    if !raw_ok && !app_ok { st.neither_pm_match += 1; }
    if b % 2 == 0 { st.b_even += 1; if raw_ok { st.raw_match_when_b_even += 1; } }
    else { st.b_odd += 1; if app_ok { st.app_match_when_b_odd += 1; } }
    // 12-byte fingerprints
    let f_on = imn_fingerprint(gf, &pms_on); let f_off = imn_fingerprint(gf, &pms_off);
    let f_on_app = imn_fingerprint(gf, &pms_on_app); let f_off_app = imn_fingerprint(gf, &pms_off_app);
    if f_on.is_some() { st.imn_defined += 1; } else { st.raw_undef += 1; }
    if p_on[1] == 0 { st.s1_zero_online += 1; }
    if b % 2 == 1 && s != 0 {
        st.app_needed += 1;
        if f_on_app.is_none() {
            st.app_needed_undef += 1;
            if p_on[1] == 0 { st.app_undef_s1zero += 1; }
        }
    }
    let raw_key = matches!((&f_on, &f_off), (Some(x), Some(y)) if x.bytes == y.bytes);
    let app_key = matches!((&f_on_app, &f_off_app), (Some(x), Some(y)) if x.bytes == y.bytes);
    if raw_key { st.imn_key_match_raw += 1; }
    if app_key { st.imn_key_match_app += 1; }
    if raw_key || app_key { st.imn_key_match_either += 1; }
    // zero-pattern agreement (matching parity)
    if s != 0 {
        st.imn_zero_pattern_checked += 1;
        let (a1, a2) = if b % 2 == 0 { (&pms_on, &pms_off) } else { (&pms_on_app, &pms_off_app) };
        if (0..15).all(|i| (a1[i] == 0) == (a2[i] == 0)) { st.imn_zero_pattern_agree += 1; }
    }
    // s recovery
    let (p7, q7) = (pms_on[0], pms_off[0]); // E15[0] = 7
    if q7 == 0 || p7 == 0 { st.s_undef += 1; }
    else {
        let s_lit = gf.root(gf.div(p7, q7), 14);
        if s_lit == s { st.s_lit_correct += 1; }
        // corrected: choose parity that matched
        let (pp, qq) = if b % 2 == 0 { (pms_on[0], pms_off[0]) } else { (pms_on_app[0], pms_off_app[0]) };
        if qq != 0 && pp != 0 { st.s_corr_defined += 1; let s_c = gf.root(gf.div(pp, qq), 14); if s_c == s { st.s_corr_correct += 1; } }
    }
    // keyless d computation from the 24 anchor bytes
    let mut ds: Vec<u8> = (1..=255u8).map(|e| keyless_d(gf, &anchor_x2, &anchor_x3, &anchor_x4, e)).collect();
    let mut real: Vec<u8> = (1..256usize).map(|v| a[0] ^ a[v]).collect();
    ds.sort_unstable(); real.sort_unstable();
    if ds == real { st.keyless_multiset_ok += 1; }
    // chi-star agreement raw / appended
    let c_on = chi_star(gf, pt, &g, 0); let c_off = chi_star(gf, pt, &dinv, 0);
    let c_on_a = chi_star(gf, pt, &g, 1); let c_off_a = chi_star(gf, pt, &dinv, 1);
    match (c_on, c_off, c_on_a, c_off_a) {
        (Some(x), Some(y), Some(xa), Some(ya)) => {
            let r = x.vec == y.vec; let ap = xa.vec == ya.vec;
            if r { st.chi_raw_match += 1; }
            if ap { st.chi_app_match += 1; }
            if r || ap { st.chi_either_match += 1; }
        }
        _ => { st.chi_undef += 1; }
    }
}

pub fn e4(gf: &Gf, instances: u64) -> String {
    let pt = PowTab::new(gf);
    let nthreads = rayon::current_num_threads() as u64;
    let per = (instances + nthreads - 1) / nthreads;
    let stats: Stats = (0..nthreads).into_par_iter().map(|t| {
        let mut rng = Rng::new(0xE4E4_0000 + t);
        let mut st = Stats::new();
        for _ in 0..per { one_instance(gf, &pt, &mut rng, &mut st); }
        st
    }).reduce(Stats::new, Stats::merge);
    let n = stats.n as f64;
    let pct = |x: u64| 100.0 * x as f64 / n;
    let mut out = String::new();
    out += &format!("E4 honest bridge on real 7-round AES-128: {} instances (each: fresh key, delta-set of 256 plaintexts built at x_1 via known k_-1)\n", stats.n);
    out += &format!("  ciphertext peel + equivalent key reproduces x_5[0]: {} / {} instances\n", stats.peel_ok, stats.n);
    out += &format!("  keyless offline d(e) from the 24 anchor bytes reproduces the multiset {{a ^ a_w}}: {} / {} instances\n", stats.keyless_multiset_ok, stats.n);
    out += &format!("  s = a_0 = 0 (bridge fails): {} = {:.4}% (paper: 2^-8 = 0.3906%)\n", stats.s_zero, pct(stats.s_zero));
    // b histogram vs Poisson(255*2/256)
    out += "  bad-slot count b = A + B histogram (paper: expect ~2 on average):\n";
    let lam: f64 = 255.0 * 2.0 / 256.0;
    let mut mean = 0f64;
    for k in 0..12 {
        let obs = stats.b_hist[k] as f64 / n;
        let mut fac = 1f64; for i in 1..=k { fac *= i as f64; }
        let poi = (-lam).exp() * lam.powi(k as i32) / fac;
        mean += k as f64 * stats.b_hist[k] as f64;
        out += &format!("     b={k}: observed {:.5}  Poisson({lam:.4}) predicts {:.5}\n", obs, poi);
    }
    out += &format!("     mean b = {:.5}\n", mean / n);
    out += &format!("  parity claim: raw fingerprint matches (P_m = s^2m Q_m for all 15 m) in {:.4}% of instances; appended-zero variant matches in {:.4}%; both {:.4}%; neither {:.4}%\n",
        pct(stats.raw_pm_match), pct(stats.app_pm_match), pct(stats.both_pm_match), pct(stats.neither_pm_match));
    out += &format!("     among b even ({} instances): raw matched {}; among b odd ({}): appended matched {}  (paper: exactly one matches unless s=0)\n",
        stats.b_even, stats.raw_match_when_b_even, stats.b_odd, stats.app_match_when_b_odd);
    out += &format!("  12-byte I_(m,n) fingerprint defined online (>=13 nonzero P_m): {:.5}%  (paper: fewer with prob <= 2^-15.2 = 0.003%)\n", pct(stats.imn_defined));
    out += &format!("  table-key equality I_on == I_off: raw {:.4}%, appended {:.4}%, either {:.4}%  (=> the true entry is hit iff s != 0)\n",
        pct(stats.imn_key_match_raw), pct(stats.imn_key_match_app), pct(stats.imn_key_match_either));
    out += &format!("  online/offline P_m zero patterns agree (matching parity, s!=0): {} / {}\n", stats.imn_zero_pattern_agree, stats.imn_zero_pattern_checked);
    out += &format!("  ALGORITHM 1 (line 1305) exact multiset test {{g_w}} == {{s^2 d^-1 + s}} on the TRUE hit: pass rate raw {:.3}%, appended {:.3}%, either {:.3}%\n     (Poisson model prediction: 30.9% raw / 52.5% either; i.e. Alg. 1 as written rejects the true hit ~48% of the time)\n",
        pct(stats.alg1_raw_pass), pct(stats.alg1_app_pass), pct(stats.alg1_either_pass));
    out += &format!("  corrected test (multiset equality after deleting the values {{0, s}} from both sides): pass rate {:.4}%\n", pct(stats.corrected_pass));
    let s_defined = stats.n - stats.s_undef;
    out += &format!("  s recovery s = (P_7/Q_7)^(1/14): paper's literal step 2 (raw P_7,Q_7) correct in {:.3}% of defined cases ({} undefined, P_7 or Q_7 = 0); parity-corrected version correct in {:.4}% of its defined cases ({} instances) [remaining failures = the s=0 instances]\n",
        100.0 * stats.s_lit_correct as f64 / s_defined as f64, stats.s_undef, 100.0 * stats.s_corr_correct as f64 / stats.s_corr_defined as f64, stats.s_corr_defined);
    out += &format!("  even-|D| correlated-zeros defect: online S_1 = 0 in {:.4}% of instances (expect 1/256 = 0.3906%); the appended fingerprint was needed (b odd, s!=0) in {} instances, of which it was UNDEFINED (<13 nonzero P_m) in {} = {:.4}%; of those, S_1 (=p_1 online) was 0 in {}\n     (the five weight-3 exponents 7,11,13,19,37 vanish together at even |D| when S_1=0, leaving 10 < 13 nonzero P_m -- the paper's line 1169-1170 assumes independent zeros, prob <= 2^-15.2)\n",
        pct(stats.s1_zero_online), stats.app_needed, stats.app_needed_undef, 100.0 * stats.app_needed_undef as f64 / stats.app_needed.max(1) as f64, stats.app_undef_s1zero);
    out += &format!("  raw (odd-|D|) online fingerprint undefined: {} of {} (paper's independence estimate 2^-15.2 -> {:.2} expected)\n", stats.raw_undef, stats.n, n * 2f64.powf(-15.2));
    out += &format!("  chi-star canonical fingerprint agrees online/offline: raw {:.4}%, appended {:.4}%, either {:.4}% (undefined: {})\n",
        pct(stats.chi_raw_match), pct(stats.chi_app_match), pct(stats.chi_either_match), stats.chi_undef);
    out
}

/// debug: find instances where the appended P_m identity holds but the
/// 12-byte fingerprints differ.
pub fn e4_debug(gf: &Gf) -> String {
    let pt = PowTab::new(gf);
    let mut rng = Rng::new(0xDEB);
    let mut out = String::new();
    let mut found = 0;
    for _ in 0..200_000 {
        let key = rng.state16(); let rk = expand_key(gf, &key);
        let base_x1 = rng.state16();
        let mut a = [0u8; 256]; let mut w = [0u8; 256];
        let pos = [0usize, 13, 10, 7]; let k6 = rk[7];
        for v in 0..256usize {
            let mut x1 = base_x1; x1[0] = base_x1[0] ^ (v as u8);
            let p = plaintext_from_x1(gf, &rk, &x1);
            let t = encrypt_trace(gf, &rk, &p, 7);
            a[v] = t.xs[5][0];
            let mut x6c0 = [0u8; 4];
            for j in 0..4 { x6c0[j] = gf.inv_sbox[(t.c[pos[j]] ^ k6[pos[j]]) as usize]; }
            w[v] = gf.mul14[x6c0[0] as usize] ^ gf.mul11[x6c0[1] as usize] ^ gf.mul13[x6c0[2] as usize] ^ gf.mul9[x6c0[3] as usize];
        }
        let s = a[0]; if s == 0 { continue; }
        let mut g = [0u8; 255]; let mut dinv = [0u8; 255];
        let (mut ca, mut cb) = (0, 0);
        for v in 1..256usize { g[v-1] = gf.inv[gf.linv[(w[0]^w[v]) as usize] as usize]; dinv[v-1] = gf.inv[(a[0]^a[v]) as usize];
            if a[v] == 0 { ca += 1; } if a[v] == a[0] { cb += 1; } }
        let b = ca + cb; if b % 2 == 0 { continue; }
        let p_on = power_sums(gf, &pt, &g, 0); let p_off = power_sums(gf, &pt, &dinv, 0);
        let mut pa = p_on; pa[0] ^= 1; let mut qa = p_off; qa[0] ^= 1;
        let pms_on = all_pms(gf, &pa); let pms_off = all_pms(gf, &qa);
        let mut ok = true;
        for i in 0..15 { if pms_on[i] != gf.mul(gf.pow(s, 2 * E15[i] as u32), pms_off[i]) { ok = false; } }
        if !ok { continue; }
        let f1 = imn_fingerprint(gf, &pms_on); let f2 = imn_fingerprint(gf, &pms_off);
        let matched = matches!((&f1, &f2), (Some(x), Some(y)) if x.bytes == y.bytes);
        if !matched {
            found += 1;
            out += &format!("case {found}: s={s:#04x} b={b} A={ca} B={cb}\n  P'on = {:?}\n  P'off= {:?}\n  f_on = {:?}\n  f_off= {:?}\n",
                pms_on, pms_off, f1.map(|f| f.bytes), f2.map(|f| f.bytes));
            if found >= 3 { break; }
        }
    }
    if found == 0 { out += "no discrepancy found in 200k instances\n"; }
    out
}
