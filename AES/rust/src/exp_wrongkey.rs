// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! E5: wrong-key randomization over all 2^32 guesses of k_6[0,7,10,13],
//! plus a reduced-scale table lookup (2^24 random offline entries) counting
//! false-positive hits.
use crate::aes::*;
use crate::bridge::*;
use crate::gf::Gf;
use crate::rng::Rng;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::Mutex;

pub fn e5(gf: &Gf, log2_guesses: u32) -> String {
    let pt = PowTab::new(gf);
    let mut rng = Rng::new(0xE5);
    let key = rng.state16();
    let rk = expand_key(gf, &key);
    // honest delta-set at x_1
    let base_x1 = rng.state16();
    let mut c_all = vec![[0u8; 16]; 256];
    let mut a = [0u8; 256];
    for v in 0..256usize {
        let mut x1 = base_x1; x1[0] = base_x1[0] ^ (v as u8);
        let p = plaintext_from_x1(gf, &rk, &x1);
        let t = encrypt_trace(gf, &rk, &p, 7);
        c_all[v] = t.c; a[v] = t.xs[5][0];
    }
    let s = a[0];
    // true offline entry: values d^-1
    let dinv: Vec<u8> = (1..256usize).map(|v| gf.inv[(a[0] ^ a[v]) as usize]).collect();
    let p_off = power_sums(gf, &pt, &dinv, 0);
    let mut p_off_app = p_off; p_off_app[0] ^= 1;
    let f_off = imn_fingerprint(gf, &all_pms(gf, &p_off)).expect("offline fingerprint undefined");
    let f_off_app = imn_fingerprint(gf, &all_pms(gf, &p_off_app)).expect("offline app fingerprint undefined");
    // per-byte peel tables T_j[gg][v] = coef_j * InvS(C^{(v)}[pos_j] ^ gg)
    let pos = [0usize, 13, 10, 7];
    let coef = [14u8, 11, 13, 9];
    let mut ttab = vec![vec![[0u8; 256]; 256]; 4]; // ttab[j][gg][v]
    for j in 0..4 { for gg in 0..256usize { for v in 0..256usize {
        ttab[j][gg][v] = gf.mul(coef[j], gf.inv_sbox[(c_all[v][pos[j]] as usize) ^ gg]);
    }}}
    let mut invlinv = [0u8; 256];
    for x in 0..256usize { invlinv[x] = gf.inv[gf.linv[x] as usize]; }
    // random offline table (reduced scale): 2^24 fingerprints of random 255-multisets, plus the true entry
    let table_log = 24u32;
    let table: HashSet<[u8; 12]> = {
        let entries: Vec<[u8; 12]> = (0..(1u64 << table_log)).into_par_iter().map(|i| {
            let mut r = Rng::new(0xEE00_0000 + i);
            loop {
                let vals: Vec<u8> = (0..255).map(|_| r.byte()).collect();
                let p = power_sums(gf, &pt, &vals, 0);
                if let Some(f) = imn_fingerprint(gf, &all_pms(gf, &p)) { return f.bytes; }
            }
        }).collect();
        let mut hs: HashSet<[u8; 12]> = entries.into_iter().collect();
        hs.insert(f_off.bytes); hs.insert(f_off_app.bytes);
        hs
    };
    let true_k6 = [rk[7][0], rk[7][13], rk[7][10], rk[7][7]];
    let n_guess = 1u64 << log2_guesses;
    let outer = 1u64 << (log2_guesses.saturating_sub(24)); // g0 range
    // stats
    let ham_hist = Mutex::new(vec![0u64; 13]);
    let hits_raw = Mutex::new(0u64); let hits_app = Mutex::new(0u64);
    let table_hits = Mutex::new(0u64); let true_found = Mutex::new(false);
    let undef = Mutex::new(0u64);
    // 3-byte prefix histogram for collision entropy (raw fingerprint, first 3 bytes)
    let hist3: Vec<Mutex<Vec<u32>>> = Vec::new(); let _ = hist3;
    let global3: Mutex<Vec<u32>> = Mutex::new(vec![0u32; 1 << 24]);
    (0..outer).into_par_iter().for_each(|g0| {
        let g0 = g0 as usize;
        let mut lh = vec![0u64; 13];
        let (mut lhr, mut lha, mut lth, mut lundef) = (0u64, 0u64, 0u64, 0u64);
        let mut ltrue = false;
        let mut lhist3 = vec![0u32; 1 << 24];
        let inner: u64 = if log2_guesses >= 24 { 1 << 24 } else { 1 << log2_guesses };
        for gi in 0..inner {
            let g1 = ((gi >> 16) & 0xff) as usize; let g2 = ((gi >> 8) & 0xff) as usize; let g3 = (gi & 0xff) as usize;
            // online values
            let t0 = &ttab[0][g0]; let t1 = &ttab[1][g1]; let t2 = &ttab[2][g2]; let t3 = &ttab[3][g3];
            let w0 = t0[0] ^ t1[0] ^ t2[0] ^ t3[0];
            // accumulate power sums directly
            let mut acc = [0u64; 8];
            let mut vals = [0u8; 255];
            for v in 1..256usize {
                let wv = t0[v] ^ t1[v] ^ t2[v] ^ t3[v];
                let gv = invlinv[(w0 ^ wv) as usize];
                vals[v - 1] = gv;
                let row = &pt.row[gv as usize];
                for wi in 0..8 { acc[wi] ^= u64::from_le_bytes(row[8 * wi..8 * wi + 8].try_into().unwrap()); }
            }
            let mut p = [0u8; 128];
            let mut odd = [0u8; 64];
            for wi in 0..8 { odd[8 * wi..8 * wi + 8].copy_from_slice(&acc[wi].to_le_bytes()); }
            p[0] = 1; // |V| = 255
            for j in 0..64 { p[2 * j + 1] = odd[j]; }
            for r in (2..128).step_by(2) { p[r] = gf.mul(p[r / 2], p[r / 2]); }
            let pms = all_pms(gf, &p);
            let mut p_app = p; p_app[0] = 0;
            let pms_app = all_pms(gf, &p_app);
            let f = imn_fingerprint(gf, &pms);
            let fa = imn_fingerprint(gf, &pms_app);
            let is_true = g0 as u8 == true_k6[0] && g1 as u8 == true_k6[1] && g2 as u8 == true_k6[2] && g3 as u8 == true_k6[3];
            if let Some(fr) = &f {
                if fr.bytes == f_off.bytes { lhr += 1; }
                let ham = (0..12).filter(|&i| fr.bytes[i] != f_off.bytes[i]).count();
                lh[ham] += 1;
                if table.contains(&fr.bytes) { lth += 1; if is_true { ltrue = true; } }
                let idx = ((fr.bytes[0] as usize) << 16) | ((fr.bytes[1] as usize) << 8) | (fr.bytes[2] as usize);
                lhist3[idx] += 1;
            } else { lundef += 1; }
            if let Some(far) = &fa {
                if far.bytes == f_off_app.bytes { lha += 1; if is_true { ltrue = true; } }
                if table.contains(&far.bytes) { lth += 1; if is_true { ltrue = true; } }
            }
            let _ = &vals;
        }
        {
            let mut h = ham_hist.lock().unwrap(); for i in 0..13 { h[i] += lh[i]; }
            *hits_raw.lock().unwrap() += lhr; *hits_app.lock().unwrap() += lha;
            *table_hits.lock().unwrap() += lth; *undef.lock().unwrap() += lundef;
            if ltrue { *true_found.lock().unwrap() = true; }
            let mut g = global3.lock().unwrap();
            for i in 0..(1 << 24) { g[i] = g[i].wrapping_add(lhist3[i]); }
        }
    });
    let ham = ham_hist.into_inner().unwrap();
    let h3 = global3.into_inner().unwrap();
    // collision entropy of 1,2,3-byte prefixes from the 3-byte histogram
    let mut c1 = [0u64; 256]; let mut c2 = vec![0u64; 65536];
    let mut coll3 = 0u128; let mut m: u64 = 0;
    for i in 0..(1usize << 24) {
        let cnt = h3[i] as u64; m += cnt;
        c1[i >> 16] += cnt; c2[i >> 8] += cnt;
        coll3 += (cnt as u128) * ((cnt.saturating_sub(1)) as u128) / 2;
    }
    let coll1: u128 = c1.iter().map(|&c| (c as u128) * ((c.saturating_sub(1)) as u128) / 2).sum();
    let coll2: u128 = c2.iter().map(|&c| (c as u128) * ((c.saturating_sub(1)) as u128) / 2).sum();
    let h2 = |c: u128| { let mm = m as f64; ((mm * (mm - 1.0)) / 2.0 / (c as f64)).log2() };
    let n_wrong = (n_guess - 1) as f64;
    // Binom(12, 254/255) predictions
    let mut binom = String::new();
    let pdiff = 254.0f64 / 255.0;
    for h in 6..=12usize {
        let comb = binomial(12, h);
        let e = comb * pdiff.powi(h as i32) * (1.0 - pdiff).powi(12 - h as i32) * n_wrong;
        let obs = ham[h];
        let z = if e > 0.0 { (obs as f64 - e) / e.sqrt() } else { 0.0 };
        binom += &format!("     h={h}: observed {obs:>13}  Binom(12,254/255)*N predicts {:>15.1}  z={:+.2}\n", e, z);
    }
    format!(
        "E5 wrong-key randomization over 2^{log2_guesses} guesses of k_6[0,7,10,13] (fixed key, honest delta-set; s = {:#04x}):\n\
         \x20 exact 12-byte fingerprint hits vs the true entry: raw {} , appended {} (the correct guess is included and must give 1 hit under one parity)\n\
         \x20 fingerprint undefined (<13 nonzero P_m): {}\n\
         \x20 byte-Hamming distance to the true entry's raw fingerprint (h<=5 rows are all zero except the true guess at h=0):\n{}\
         \x20     rows h=0..5: {:?}\n\
         \x20 prefix collision entropy over the {} defined online fingerprints: k=1: {:.4} (ideal 7.9944), k=2: {:.4} (ideal 15.9887), k=3: {:.4} (ideal 23.9831); paper reports 7.9944, 15.9887, 23.9831.\n\
         \x20 reduced-scale table test: {} online fingerprints (raw+app) looked up in a table of 2^{table_log} random-multiset fingerprints + the true entry (both parities): total hits = {} (expected FP = 2^{} * 2^{table_log} / 2^96 = 2^{}); true (guess,entry) hit found: {}\n",
        s, hits_raw.into_inner().unwrap(), hits_app.into_inner().unwrap(), undef.into_inner().unwrap(), binom,
        &ham[0..6], m, h2(coll1), h2(coll2), h2(coll3), 2 * n_guess,
        table_hits.into_inner().unwrap(), log2_guesses + 1, (log2_guesses as i64 + 1 + table_log as i64 - 96), true_found.into_inner().unwrap())
}

fn binomial(n: usize, k: usize) -> f64 {
    let mut r = 1.0f64;
    for i in 0..k { r = r * (n - i) as f64 / (i + 1) as f64; }
    r
}

/// E5c: wrong-key randomization for the chi-star fingerprint (13-byte prefix).
pub fn e5_chi(gf: &Gf, log2_guesses: u32) -> String {
    let pt = PowTab::new(gf);
    let mut rng = Rng::new(0xE5C);
    let key = rng.state16();
    let rk = expand_key(gf, &key);
    let base_x1 = rng.state16();
    let mut c_all = vec![[0u8; 16]; 256];
    let mut a = [0u8; 256];
    for v in 0..256usize {
        let mut x1 = base_x1; x1[0] = base_x1[0] ^ (v as u8);
        let p = plaintext_from_x1(gf, &rk, &x1);
        let t = encrypt_trace(gf, &rk, &p, 7);
        c_all[v] = t.c; a[v] = t.xs[5][0];
    }
    let s = a[0];
    let dinv: Vec<u8> = (1..256usize).map(|v| gf.inv[(a[0] ^ a[v]) as usize]).collect();
    let f_off = chi_star(gf, &pt, &dinv, 0).expect("chi off undefined");
    let f_off_app = chi_star(gf, &pt, &dinv, 1).expect("chi off app undefined");
    let pos = [0usize, 13, 10, 7];
    let coef = [14u8, 11, 13, 9];
    let mut ttab = vec![vec![[0u8; 256]; 256]; 4];
    for j in 0..4 { for gg in 0..256usize { for v in 0..256usize {
        ttab[j][gg][v] = gf.mul(coef[j], gf.inv_sbox[(c_all[v][pos[j]] as usize) ^ gg]);
    }}}
    let mut invlinv = [0u8; 256];
    for x in 0..256usize { invlinv[x] = gf.inv[gf.linv[x] as usize]; }
    let true_k6 = [rk[7][0], rk[7][13], rk[7][10], rk[7][7]];
    let outer = 1u64 << (log2_guesses.saturating_sub(24));
    let ham_hist = Mutex::new(vec![0u64; 14]);
    let hits = Mutex::new(0u64);
    let global3: Mutex<Vec<u32>> = Mutex::new(vec![0u32; 1 << 24]);
    (0..outer).into_par_iter().for_each(|g0| {
        let g0 = g0 as usize;
        let mut lh = vec![0u64; 14];
        let mut lhits = 0u64;
        let mut lhist3 = vec![0u32; 1 << 24];
        let inner: u64 = if log2_guesses >= 24 { 1 << 24 } else { 1 << log2_guesses };
        let mut vals = [0u8; 255];
        for gi in 0..inner {
            let g1 = ((gi >> 16) & 0xff) as usize; let g2 = ((gi >> 8) & 0xff) as usize; let g3 = (gi & 0xff) as usize;
            let t0 = &ttab[0][g0]; let t1 = &ttab[1][g1]; let t2 = &ttab[2][g2]; let t3 = &ttab[3][g3];
            let w0 = t0[0] ^ t1[0] ^ t2[0] ^ t3[0];
            for v in 1..256usize { vals[v - 1] = invlinv[(w0 ^ t0[v] ^ t1[v] ^ t2[v] ^ t3[v]) as usize]; }
            let is_true = g0 as u8 == true_k6[0] && g1 as u8 == true_k6[1] && g2 as u8 == true_k6[2] && g3 as u8 == true_k6[3];
            if let Some(f) = chi_star(gf, &pt, &vals, 0) {
                if f.vec[..13] == f_off.vec[..13] { lhits += 1; }
                let ham = (0..13).filter(|&i| f.vec[i] != f_off.vec[i]).count();
                lh[ham] += 1;
                let idx = ((f.vec[0] as usize) << 16) | ((f.vec[1] as usize) << 8) | (f.vec[2] as usize);
                lhist3[idx] += 1;
            }
            if let Some(fa) = chi_star(gf, &pt, &vals, 1) {
                if fa.vec[..13] == f_off_app.vec[..13] { lhits += 1; }
            }
            let _ = is_true;
        }
        {
            let mut h = ham_hist.lock().unwrap(); for i in 0..14 { h[i] += lh[i]; }
            *hits.lock().unwrap() += lhits;
            let mut g = global3.lock().unwrap();
            for i in 0..(1 << 24) { g[i] = g[i].wrapping_add(lhist3[i]); }
        }
    });
    let ham = ham_hist.into_inner().unwrap();
    let h3 = global3.into_inner().unwrap();
    let mut c1 = [0u64; 256]; let mut c2 = vec![0u64; 65536];
    let mut coll3 = 0u128; let mut m: u64 = 0;
    for i in 0..(1usize << 24) {
        let cnt = h3[i] as u64; m += cnt;
        c1[i >> 16] += cnt; c2[i >> 8] += cnt;
        coll3 += (cnt as u128) * ((cnt.saturating_sub(1)) as u128) / 2;
    }
    let coll1: u128 = c1.iter().map(|&c| (c as u128) * ((c.saturating_sub(1)) as u128) / 2).sum();
    let coll2: u128 = c2.iter().map(|&c| (c as u128) * ((c.saturating_sub(1)) as u128) / 2).sum();
    let h2 = |c: u128| { let mm = m as f64; ((mm * (mm - 1.0)) / 2.0 / (c as f64)).log2() };
    format!(
        "E5c chi-star wrong-key run over 2^{log2_guesses} guesses (s = {:#04x}):\n\
         \x20 exact 13-byte prefix hits vs the true entry (raw+appended parities): {} (the correct guess must contribute exactly 1)\n\
         \x20 byte-Hamming distance over the 13-byte prefix, rows h=0..7: {:?}\n\
         \x20 prefix collision entropy of the online chi-star vectors: k=1: {:.4}, k=2: {:.4}, k=3: {:.4} (independence prediction 7.7905, 15.5811, 23.3716)\n",
        s, hits.into_inner().unwrap(), &ham[0..8], h2(coll1), h2(coll2), h2(coll3))
}
