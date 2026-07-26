// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! E1 bridge identity, E2 sigma_xi group facts, E3 exponent/coset facts,
//! E8 DDT statistics, E9 AES sanity + differential-shape checks.
use crate::aes::*;
use crate::bridge::*;
use crate::gf::Gf;
use crate::rng::Rng;

pub fn e1(gf: &Gf) -> String {
    // Bridge identity: online g = (L^{-1}(v0 ^ vw))^{-1} with v = S(a) ^ kappa,
    // vw = S(e) ^ kappa; formula: g = s^2 d^{-1} + s, s = a, d = a ^ e.
    let mut ok = 0u64; let mut total = 0u64;
    for a in 1..=255u8 { for e in 1..=255u8 {
        if a == e { continue; }
        total += 1;
        for kappa in [0u8, 0x5a, 0xff] {
            let v0 = gf.sbox[a as usize] ^ kappa;
            let vw = gf.sbox[e as usize] ^ kappa;
            let g = gf.inv[gf.linv[(v0 ^ vw) as usize] as usize];
            let d = a ^ e;
            let pred = gf.mul(gf.mul(a, a), gf.inv[d as usize]) ^ a;
            if kappa == 0 { if g == pred { ok += 1; } } else { assert_eq!(g == pred, true); }
        }
    }}
    // corner cases: e = 0 (a_omega = 0) and e = a (a_omega = s = a)
    let mut typea = String::new(); let mut typeb = String::new();
    let a = 0x37u8;
    {
        let e = 0u8; // a_w = 0: online g uses 0^{-1} := 0
        let v0 = gf.sbox[a as usize]; let vw = gf.sbox[e as usize];
        let g = gf.inv[gf.linv[(v0 ^ vw) as usize] as usize];
        let d = a ^ e; let pred = gf.mul(gf.mul(a, a), gf.inv[d as usize]) ^ a;
        typea = format!("a_w=0: online g = {:#04x} (s = a = {:#04x}), formula gives {:#04x}", g, a, pred);
    }
    {
        let e = a; // a_w = s: d = 0
        let v0 = gf.sbox[a as usize]; let vw = gf.sbox[e as usize];
        let g = gf.inv[gf.linv[(v0 ^ vw) as usize] as usize];
        let d = a ^ e; let pred = gf.mul(gf.mul(a, a), gf.inv[d as usize]) ^ a;
        typeb = format!("a_w=s: online g = {:#04x}, formula gives {:#04x} (= s)", g, pred);
    }
    // check universality of the bad-slot value swap over all a
    let mut swap_ok = true;
    for a in 1..=255u8 {
        let ga = gf.inv[gf.linv[(gf.sbox[a as usize] ^ gf.sbox[0]) as usize] as usize]; // e=0
        let pa = gf.mul(gf.mul(a, a), gf.inv[a as usize]) ^ a; // d=a
        if ga != a || pa != 0 { swap_ok = false; }
        // e=a: online v0^vw = 0 -> g = 0; formula d=0 -> 0^{-1}=0 -> pred = a
        if gf.inv[gf.linv[0] as usize] != 0 { swap_ok = false; }
    }
    format!(
        "E1 bridge identity: {ok}/{total} pairs (a,e) in GF(256)^*x GF(256)^*, a!=e satisfy g = s^2 d^-1 + s (also re-checked with kappa in {{0x5a,0xff}}: kappa cancels).\n\
         paper claims 64770/64770 = 254*255 = {}.\n\
         bad slots: {typea}\n           {typeb}\n\
         universal swap check (online value at a_w=0 is s, formula 0; at a_w=s online 0, formula s) for all 255 s: {swap_ok}\n",
        254 * 255)
}

pub fn e2(gf: &Gf) -> String {
    // sigma_xi(v) = S(xi ^ v) ^ S(xi) as a permutation of the 255 nonzero v.
    let cyc_type = |perm: &[u8; 256]| -> Vec<usize> {
        let mut seen = [false; 256]; let mut out = vec![];
        for v in 1..=255usize { if seen[v] { continue; }
            let mut c = 0; let mut x = v;
            while !seen[x] { seen[x] = true; x = perm[x] as usize; c += 1; }
            out.push(c); }
        out.sort_by(|a, b| b.cmp(a)); out
    };
    let sigma = |xi: u8| -> [u8; 256] {
        let mut p = [0u8; 256];
        for v in 0..256usize { p[v] = gf.sbox[(xi as usize) ^ v] ^ gf.sbox[xi as usize]; }
        p
    };
    let s0 = sigma(0);
    let ct = cyc_type(&s0);
    let mut all_odd = true;
    let mut n_even = 0;
    for xi in 0..=255u8 {
        let ct2 = cyc_type(&sigma(xi));
        let transpositions: usize = ct2.iter().map(|c| c - 1).sum();
        if transpositions % 2 == 0 { all_odd = false; n_even += 1; }
    }
    // pow of sigma_0 giving a 7-cycle: 13446 = lcm(166,81) = 166*81 (coprime)
    let l = 166 * 81;
    format!(
        "E2 sigma_xi facts: cycle type of sigma_0 on the 255 nonzero indices = {:?} (paper: (166,81,7,1)); lcm(166,81) = {} (paper: sigma_0^13446 is a 7-cycle: 13446 = {}).\n\
         all 256 sigma_xi are odd permutations: {} ({} even ones found).\n",
        ct, l, l, all_odd, n_even)
}

pub fn e3(gf: &Gf) -> String {
    // cyclotomic cosets of exponents coprime to 255 under x -> 2x mod 255
    let gcd = |a: u32, b: u32| { let (mut a, mut b) = (a, b); while b != 0 { let t = a % b; a = b; b = t; } a };
    let mut seen = [false; 256];
    let mut reps = vec![];
    let mut ncop = 0;
    for m in 1..255u32 {
        if gcd(m, 255) != 1 { continue; }
        ncop += 1;
        if seen[m as usize] { continue; }
        let mut sz = 0; let mut x = m;
        while !seen[x as usize] { seen[x as usize] = true; x = (2 * x) % 255; sz += 1; }
        reps.push((m, sz));
    }
    // P_1 == 0 identity and pm formula cross-check on random data
    let mut rng = Rng::new(42);
    let pt = PowTab::new(gf);
    let mut mismatch = 0; let mut p1nonzero_odd = 0; let mut checks = 0;
    let mut sq_fail = 0;
    for trial in 0..2000 {
        let n = if trial % 2 == 0 { 255 } else { 256 };
        let vals: Vec<u8> = (0..n).map(|_| rng.byte()).collect();
        let p = power_sums(gf, &pt, &vals, 0);
        for &m in &[1u8, 3, 5, 7, 11, 13, 19, 23, 29, 31, 37, 43, 47, 53, 59, 61, 91, 127] {
            let a = pm_from_p(gf, m, &p);
            let b = pm_direct(gf, m, &vals);
            checks += 1;
            if a != b { mismatch += 1; }
            if m == 1 && n == 255 && b != 0 { p1nonzero_odd += 1; }
        }
        // Frobenius: P_{2m} = P_m^2 (direct)
        for &m in &[7u8, 11, 13] {
            let a = pm_direct(gf, 2 * m, &vals);
            let b = gf.mul(pm_direct(gf, m, &vals), pm_direct(gf, m, &vals));
            if a != b { sq_fail += 1; }
        }
    }
    let mut s = String::new();
    s += &format!("E3 exponent facts:\n  #{{m in [1,255): gcd(m,255)=1}} = {} (phi(255)=128); Frobenius cosets: {} (paper says exactly 15). representatives (rep,size): {:?}\n",
        ncop, reps.len(), reps);
    s += &format!("  P_m combine formula (Lucas/Reduction 1, with the [wt(m)=1] term) vs direct O(n^2) sum on 2000 random multisets (|V| in {{255,256}}), 18 exponents each: mismatches = {mismatch} / {checks}\n");
    s += &format!("  P_1 (m=1) nonzero cases at |V|=255: {p1nonzero_odd} (paper: identically 0 by parity cancellation)\n");
    s += &format!("  Frobenius P_2m = P_m^2 failures: {sq_fail}\n");
    // 12-exponent 'optimal' set gcd check
    let opt = [7u8, 11, 13, 15, 23, 27, 39, 43, 31, 47, 55, 59];
    let bad: Vec<(u8, u32)> = opt.iter().map(|&m| (m, gcd(m as u32, 255))).filter(|&(_, g)| g != 1).collect();
    s += &format!("  paper's 'optimal' 12-exponent set (line 2108): elements sharing a factor with 255: {:?}\n", bad);
    // I_{m,2n} = I_{m,n}^2 identity
    let mut idfail = 0;
    for _ in 0..500 {
        let vals: Vec<u8> = (0..255).map(|_| rng.byte()).collect();
        let p = power_sums(gf, &pt, &vals, 0);
        let (m, n) = (7u8, 11u8);
        let pm7 = pm_from_p(gf, m, &p); let pn = pm_from_p(gf, n, &p);
        if pn == 0 { continue; }
        let imn = gf.div(gf.pow(pm7, n as u32), gf.pow(pn, m as u32));
        // I_{m,2n} using P_{2n} = P_n^2
        let p2n = gf.mul(pn, pn);
        let im2n = gf.div(gf.pow(pm7, 2 * n as u32), gf.pow(p2n, m as u32));
        if im2n != gf.mul(imn, imn) { idfail += 1; }
    }
    s += &format!("  I_(m,2n) = I_(m,n)^2 identity failures over 500 samples: {idfail}\n");
    s
}

pub fn e8(gf: &Gf) -> String {
    // DDT of the AES S-box: rows sum to 256, entries in {0,2,4}, spectrum per row
    let mut hist = [[0u32; 5]; 256]; // hist[delta][entry value 0/2/4]
    let mut spectrum_ok = true;
    for delta in 1..256usize {
        let mut row = [0u32; 256];
        for x in 0..256usize {
            let g = gf.sbox[x] ^ gf.sbox[x ^ delta];
            row[g as usize] += 1;
        }
        for gamma in 0..256usize {
            match row[gamma] { 0 => hist[delta][0] += 1, 2 => hist[delta][2] += 1, 4 => hist[delta][4] += 1, _ => spectrum_ok = false }
        }
    }
    let mut counts_ok = true;
    for delta in 1..256usize {
        if hist[delta][0] != 129 || hist[delta][2] != 126 || hist[delta][4] != 1 { counts_ok = false; }
    }
    // c_iSB closed form and W-cache Monte Carlo
    let closed = 4.0 * (127.0f64 / 256.0).powi(3);
    // Monte Carlo: random anti-diagonal ciphertext difference dc[0..4] and delta z_5[0]:
    // required delta x_6 column = MC applied to (delta,0,0,0) = delta*(2,1,1,3)
    let mut rng = Rng::new(7);
    let (mut num, mut den) = (0f64, 0f64); // sum |U_r| built vs candidates served
    let coefs = [2u8, 1, 1, 3];
    let mut ddt = vec![[0u8; 256]; 256];
    for d in 0..256usize { for x in 0..256usize { ddt[d][(gf.sbox[x] ^ gf.sbox[x ^ d]) as usize] += 1; } }
    let (mut alive, mut cands_per_pair, mut trials) = (0f64, 0f64, 0f64);
    for _ in 0..500_000 {
        // random nonzero input differences at the 4 anti-diagonal ciphertext bytes
        let mut dc = [0u8; 4];
        for r in 0..4 { dc[r] = 1 + (rng.byte() % 255); }
        trials += 1.0;
        // for each of the 255 output-difference groups delta (= Delta z_5[0]):
        // |U_r(delta)| = ddt[dc_r][coef_r*delta]; a group is alive iff all four are nonzero
        for delta in 1..=255u8 {
            let mut prod = 1f64; let mut su = 0f64;
            for r in 0..4 {
                let out = gf.mul(coefs[r], delta);
                let n = ddt[dc[r] as usize][out as usize] as f64;
                prod *= n; su += n;
            }
            if prod > 0.0 { alive += 1.0; num += su; den += prod; }
            cands_per_pair += prod;
        }
    }
    // The paper's model: E[sum_r |U_r|]/E[#cands] = 4*(256/127)/(256/127)^4 over surviving groups
    let model = 4.0 * (256.0f64 / 127.0) / (256.0f64 / 127.0).powi(4);
    format!(
        "E8 DDT statistics: entries in {{0,2,4}} for every nonzero row: {spectrum_ok}; every row has exactly 129 zeros, 126 twos, one four: {counts_ok}\n\
         c_iSB^on closed form 4*(127/256)^3 = {closed:.6} = 4*(256/127)/(256/127)^4 = {model:.6} (paper: 0.48837)\n\
         Monte-Carlo over 5e5 random pairs (nonzero anti-diagonal ciphertext differences): E[sum_r |U_r|]/E[prod_r |U_r|] over surviving delta-groups = {:.5}; surviving groups per pair = {:.3} (of 255); u_6 candidates per pair = {:.2} (paper: ~2^8)\n\
         P[a base has >=1 four-solution position among 20] = 1-(126/127)^20 = {:.4} (paper: ~14.6%); expected 4-count per position given a solution = 1/127 = {:.5}\n",
        num / den, alive / trials, cands_per_pair / trials, 1.0 - (126.0f64 / 127.0).powi(20), 1.0 / 127.0)
}

pub fn e9(gf: &Gf) -> String {
    let mut s = String::new();
    // FIPS-197 Appendix B example, AES-128, 10 rounds
    let key: State = [0x2b,0x7e,0x15,0x16,0x28,0xae,0xd2,0xa6,0xab,0xf7,0x15,0x88,0x09,0xcf,0x4f,0x3c];
    let pt: State = [0x32,0x43,0xf6,0xa8,0x88,0x5a,0x30,0x8d,0x31,0x31,0x98,0xa2,0xe0,0x37,0x07,0x34];
    let expect: State = [0x39,0x25,0x84,0x1d,0x02,0xdc,0x09,0xfb,0xdc,0x11,0x85,0x97,0x19,0x6a,0x0b,0x32];
    let rk = expand_key(gf, &key);
    let ct = encrypt(gf, &rk, &pt, 10);
    s += &format!("E9 AES sanity: FIPS-197 App.B AES-128 (10 rounds): got {} expected {} -> {}\n",
        hex(&ct), hex(&expect), if ct == expect { "MATCH" } else { "MISMATCH" });
    // FIPS-197 App. C.1 vector
    let key2: State = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
    let pt2: State = [0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,0x88,0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xff];
    let ex2: State = [0x69,0xc4,0xe0,0xd8,0x6a,0x7b,0x04,0x30,0xd8,0xcd,0xb7,0x80,0x70,0xb4,0xc5,0x5a];
    let rk2 = expand_key(gf, &key2);
    let ct2 = encrypt(gf, &rk2, &pt2, 10);
    s += &format!("   FIPS-197 App.C.1: got {} expected {} -> {}\n", hex(&ct2), hex(&ex2), if ct2 == ex2 { "MATCH" } else { "MISMATCH" });
    // MC/MC^{-1}, SR/SR^{-1} round trip and MC^{-1} row-0 coefficients (14,11,13,9)
    let mut rng = Rng::new(11);
    let mut rt_ok = true;
    for _ in 0..1000 {
        let st = rng.state16();
        if inv_mix_columns(gf, &mix_columns(gf, &st)) != st { rt_ok = false; }
        if inv_shift_rows(&shift_rows(&st)) != st { rt_ok = false; }
    }
    s += &format!("   MC/SR round trips ok: {rt_ok}\n");
    // differential shape checks on real 7-round AES
    // (a) plaintext structure varying diagonal {0,5,10,15} => x_1 differences confined to column 0
    let key = rng.state16(); let rk = expand_key(gf, &key);
    let mut colok = true;
    for _ in 0..1000 {
        let base = rng.state16();
        let mut p1 = base; let mut p2 = base;
        for &j in &[0usize, 5, 10, 15] { p1[j] = rng.byte(); p2[j] = rng.byte(); }
        let t1 = encrypt_trace(gf, &rk, &p1, 7); let t2 = encrypt_trace(gf, &rk, &p2, 7);
        let d1 = xor(&t1.xs[1], &t2.xs[1]);
        for j in 4..16 { if d1[j] != 0 { colok = false; } }
    }
    s += &format!("   plaintexts differing only on the diagonal {{0,5,10,15}} => Delta x_1 confined to column 0: {colok}\n");
    // (b) states x_5 differing in byte 0 only => ciphertext difference exactly at {0,7,10,13}
    let mut supp_ok = true; let mut always4 = true;
    for _ in 0..2000 {
        // pick x_5 states directly by choosing x_1 ... simpler: choose plaintext, get x_5, modify byte 0, run remaining rounds
        let p = rng.state16();
        let t = encrypt_trace(gf, &rk, &p, 7);
        let mut x5b = t.xs[5]; x5b[0] ^= 1 + (rng.byte() % 255);
        // continue rounds 5,6 from x5b
        let cont = |x5: &State| -> State {
            let y5 = sub_bytes(gf, x5); let z5 = shift_rows(&y5); let w5 = mix_columns(gf, &z5); let x6 = xor(&w5, &rk[6]);
            let y6 = sub_bytes(gf, &x6); let z6 = shift_rows(&y6); xor(&z6, &rk[7])
        };
        let c1 = cont(&t.xs[5]); let c2 = cont(&x5b);
        assert_eq!(c1, t.c);
        let dc = xor(&c1, &c2);
        for j in 0..16 {
            let active = dc[j] != 0;
            let should = j == 0 || j == 7 || j == 10 || j == 13;
            if active && !should { supp_ok = false; }
            if !active && should { always4 = false; }
        }
    }
    s += &format!("   Delta x_5 = e_0*delta => Delta C supported on {{0,7,10,13}}: {supp_ok}; all four always active (MC branch number): {always4}\n");
    // (c) equivalent-key peel: v_w := MC^{-1}(x_6)[0] and x_5[0] = S^{-1}(v_w ^ u_5[0]) with u_5 = MC^{-1}(k_5)
    let mut peel_ok = true;
    let u5 = inv_mix_columns(gf, &rk[6]);
    for _ in 0..1000 {
        let p = rng.state16();
        let t = encrypt_trace(gf, &rk, &p, 7);
        // peel with true k_6 at anti-diagonal {0,13,10,7}
        let c = t.c; let k6 = &rk[7];
        let mut x6col0 = [0u8; 4];
        let pos = [0usize, 13, 10, 7];
        for j in 0..4 { x6col0[j] = gf.inv_sbox[(c[pos[j]] ^ k6[pos[j]]) as usize]; }
        for j in 0..4 { if x6col0[j] != t.xs[6][j] { peel_ok = false; } }
        let v = gf.mul14[x6col0[0] as usize] ^ gf.mul11[x6col0[1] as usize] ^ gf.mul13[x6col0[2] as usize] ^ gf.mul9[x6col0[3] as usize];
        let x50 = gf.inv_sbox[(v ^ u5[0]) as usize];
        if x50 != t.xs[5][0] { peel_ok = false; }
    }
    s += &format!("   ciphertext peel with k_6[0,13,10,7] gives x_6[0..3], and S^-1(MC^-1(x_6)[0] ^ u_5[0]) = x_5[0] with u_5 = MC^-1(k_5): {peel_ok}\n");
    s
}
