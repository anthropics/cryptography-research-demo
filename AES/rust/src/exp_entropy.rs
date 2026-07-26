// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! E6: prefix collision entropy of the I_{m,n} and chi-star fingerprints
//! over uniformly random multisets (odd |D|=255 and even/appended |D|=256).
//! E7: chi-star AGL-invariance and degenerate-branch frequencies.
use crate::bridge::*;
use crate::gf::Gf;
use crate::rng::Rng;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct EntropyResult {
    pub m: u64,
    pub h2: Vec<(usize, u64, f64, f64)>, // (k, C_k, H2, sigma)
    pub degenerate: u64,
    pub extra: Vec<(String, u64)>,
}

/// which fingerprint to measure
#[derive(Clone, Copy, PartialEq)]
pub enum Kind { Imn, ChiStar }

pub fn measure(gf: &Gf, pt: &PowTab, kind: Kind, appended: bool, n_samples: u64, seed: u64, off: usize, order: &[usize]) -> EntropyResult {
    let extra_zeros = if appended { 1 } else { 0 };
    let mut buf: Vec<u64> = vec![0u64; n_samples as usize];
    let degenerate = AtomicU64::new(0);
    let chi_level = [const { AtomicU64::new(0) }; 8];
    let chi_branch2 = AtomicU64::new(0);
    let chi_two_roots = AtomicU64::new(0);
    let chunk = 1_000_000usize;
    buf.par_chunks_mut(chunk).enumerate().for_each(|(ci, ch)| {
        let mut rng = Rng::new(seed ^ ((ci as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)));
        let mut vals = [0u8; 255];
        for slot in ch.iter_mut() {
            loop {
                for v in vals.iter_mut() { *v = rng.byte(); }
                match kind {
                    Kind::Imn => {
                        let p = power_sums(gf, pt, &vals, extra_zeros);
                        let pms = all_pms(gf, &p);
                        match imn_fingerprint(gf, &pms) {
                            Some(f) => { let mut b = [0u8; 8]; if !order.is_empty() { for (i, &oi) in order.iter().take(8).enumerate() { b[i] = f.bytes[oi]; } } else { let hi = (off + 8).min(12); b[..hi - off].copy_from_slice(&f.bytes[off..hi]); } *slot = u64::from_be_bytes(b); break; }
                            None => { degenerate.fetch_add(1, Ordering::Relaxed); }
                        }
                    }
                    Kind::ChiStar => {
                        match chi_star(gf, pt, &vals, extra_zeros) {
                            Some(c) => {
                                if c.chain_level < 8 { chi_level[c.chain_level].fetch_add(1, Ordering::Relaxed); }
                                if c.branch == 2 { chi_branch2.fetch_add(1, Ordering::Relaxed); }
                                if c.two_roots { chi_two_roots.fetch_add(1, Ordering::Relaxed); }
                                let mut b = [0u8; 8]; b.copy_from_slice(&c.vec[off..off + 8]); *slot = u64::from_be_bytes(b); break;
                            }
                            None => { degenerate.fetch_add(1, Ordering::Relaxed); }
                        }
                    }
                }
            }
        }
    });
    buf.par_sort_unstable();
    // prefix collision counts for k = 1..8 in one pass
    let mut coll = [0u128; 8];
    let mut runlen = [1u64; 8];
    for i in 1..buf.len() {
        let x = buf[i - 1] ^ buf[i];
        let common_bytes = if x == 0 { 8 } else { (x.leading_zeros() / 8) as usize };
        for k in 0..8 {
            if k < common_bytes { runlen[k] += 1; }
            else { let l = runlen[k] as u128; coll[k] += l * (l - 1) / 2; runlen[k] = 1; }
        }
    }
    for k in 0..8 { let l = runlen[k] as u128; coll[k] += l * (l - 1) / 2; }
    let m = buf.len() as f64;
    let mut h2 = vec![];
    for k in 0..8 {
        let c = coll[k];
        if c == 0 { h2.push((k + 1, 0, f64::INFINITY, 0.0)); continue; }
        let est = ((m * (m - 1.0)) / 2.0 / (c as f64)).log2();
        let sigma = 1.0 / ((c as f64).sqrt() * std::f64::consts::LN_2);
        h2.push((k + 1, c as u64, est, sigma));
    }
    let mut extra = vec![];
    if kind == Kind::ChiStar {
        for i in 0..6 { extra.push((format!("alpha* from P_{} (chain level {})", CHI_CHAIN[i], i), chi_level[i].load(Ordering::Relaxed))); }
        extra.push(("beta*: even n, S_1=0 -> min over 256 beta".to_string(), chi_branch2.load(Ordering::Relaxed)));
        extra.push(("beta*: even n, two Artin-Schreier roots (min-hash of 2)".to_string(), chi_two_roots.load(Ordering::Relaxed)));
    }
    EntropyResult { m: buf.len() as u64, h2, degenerate: degenerate.load(Ordering::Relaxed), extra }
}

pub fn e6_report(name: &str, res: &EntropyResult, ideal_per_byte: f64, ideal_note: &str) -> String {
    let mut s = format!("{name}: M = {} samples ({} degenerate resamples). ideal per byte = {ideal_per_byte:.4} bits ({ideal_note}).\n", res.m, res.degenerate);
    s += "     k   C_k(colliding pairs)      H2(k) [bits]     ideal        gap\n";
    for &(k, c, h, sig) in &res.h2 {
        if c == 0 { s += &format!("     {k}   {:>18}   (no collisions)\n", 0); continue; }
        s += &format!("     {k}   {:>18}   {:>10.4} +- {:.4}   {:>8.4}   {:>+8.4}\n", c, h, sig, ideal_per_byte * k as f64, h - ideal_per_byte * k as f64);
    }
    for (label, cnt) in &res.extra { s += &format!("     [{label}]: {cnt}\n"); }
    s
}

pub fn e7(gf: &Gf, samples: u64) -> String {
    let pt = PowTab::new(gf);
    let nthreads = rayon::current_num_threads() as u64;
    let per = samples / nthreads;
    // (chi mismatches raw, chi undefined, I mismatches raw, appended-chi mismatches [expected: not invariant])
    let (fail_odd, undef, imn_fail, app_fail, fail_even) = (0..nthreads).into_par_iter().map(|t| {
        let mut rng = Rng::new(0xE7 ^ t.wrapping_mul(1315423911));
        let (mut fo, mut un, mut fi, mut fa, mut fe) = (0u64, 0u64, 0u64, 0u64, 0u64);
        let mut vals = [0u8; 255];
        let mut vals2 = [0u8; 255];
        for _ in 0..per {
            for v in vals.iter_mut() { *v = rng.byte(); }
            let alpha = 1 + (rng.byte() % 255); let beta = rng.byte();
            for i in 0..255 { vals2[i] = gf.mul(alpha, vals[i]) ^ beta; }
            let a = chi_star(gf, &pt, &vals, 0); let b = chi_star(gf, &pt, &vals2, 0);
            match (a, b) { (Some(x), Some(y)) => { if x.vec != y.vec { fo += 1 } }, _ => un += 1 }
            let p1 = power_sums(gf, &pt, &vals, 0); let p2 = power_sums(gf, &pt, &vals2, 0);
            let f1 = imn_fingerprint(gf, &all_pms(gf, &p1)); let f2 = imn_fingerprint(gf, &all_pms(gf, &p2));
            if let (Some(x), Some(y)) = (f1, f2) { if x.bytes != y.bytes { fi += 1 } }
            // appended (synthetic zero appended to both) is NOT literally AGL-invariant:
            let aa = chi_star(gf, &pt, &vals, 1); let ba = chi_star(gf, &pt, &vals2, 1);
            if let (Some(x), Some(y)) = (aa, ba) { if x.vec != y.vec { fa += 1 } }
            // genuine even-|V| invariance: 256-element multisets transformed as a whole
            let mut e1 = [0u8; 256]; let mut e2 = [0u8; 256];
            for i in 0..256 { e1[i] = rng.byte(); e2[i] = gf.mul(alpha, e1[i]) ^ beta; }
            let ea = chi_star(gf, &pt, &e1, 0); let eb = chi_star(gf, &pt, &e2, 0);
            match (ea, eb) { (Some(x), Some(y)) => { if x.vec != y.vec { fe += 1 } }, _ => un += 1 }
        }
        (fo, un, fi, fa, fe)
    }).reduce(|| (0, 0, 0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3, a.4 + b.4));
    format!(
        "E7 AGL(1,256) invariance V -> alpha V + beta over {} random (V, alpha, beta), |V|=255:\n\
         \x20 chi-star canonical vector mismatches: {} (undefined: {}); I_(m,n) 12-byte fingerprint mismatches: {}\n\
         \x20 chi-star canonical vector mismatches for genuine 256-element multisets (even |V|, both roots / min-hash, S_1=0 branch): {}\n\
         \x20 [control] with a synthetic zero appended to both sides (which is NOT an AGL-covariant operation), chi-star mismatches: {} of {} -- as expected the appended variant is not literally invariant; it works only via the bad-slot parity trick tested in E4.\n",
        per * nthreads, fail_odd, undef, imn_fail, fail_even, app_fail, per * nthreads)
}

/// E12: direct query-vs-table false-positive count on 5-byte fingerprint
/// prefixes, to validate the collision-entropy -> FP-rate model.
pub fn e12(gf: &Gf, log2_n: u32) -> String {
    let pt = PowTab::new(gf);
    let n = 1u64 << log2_n;
    let mut out = String::new();
    for &appended in &[false, true] {
        let ez = if appended { 1 } else { 0 };
        let gen = |seed: u64| -> Vec<u64> {
            let mut buf = vec![0u64; n as usize];
            buf.par_chunks_mut(1_000_000).enumerate().for_each(|(ci, ch)| {
                let mut rng = Rng::new(seed ^ ((ci as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)));
                let mut vals = [0u8; 255];
                for slot in ch.iter_mut() {
                    loop {
                        for v in vals.iter_mut() { *v = rng.byte(); }
                        let p = power_sums(gf, &pt, &vals, ez);
                        if let Some(f) = imn_fingerprint(gf, &all_pms(gf, &p)) {
                            let mut b = [0u8; 8]; b[..5].copy_from_slice(&f.bytes[..5]);
                            *slot = u64::from_be_bytes(b); break;
                        }
                    }
                }
            });
            buf
        };
        let mut table = gen(0x1200 + ez as u64);
        let queries = gen(0x3400 + ez as u64);
        table.par_sort_unstable();
        // count query hits: for each query prefix, number of equal table prefixes
        let hits: u64 = queries.par_iter().map(|&q| {
            let lo = table.partition_point(|&x| x < q);
            let hi = table.partition_point(|&x| x <= q);
            (hi - lo) as u64
        }).sum();
        // predicted from measured H2(5): odd 39.972, even 36.825 (E6a/b at N=3e9)
        let h2 = if appended { 36.8252 } else { 39.9724 };
        let pred = (n as f64) * (n as f64) * 2f64.powf(-h2);
        let ideal = (n as f64) * (n as f64) * 2f64.powf(-5.0 * 255f64.log2());
        out += &format!("E12 {} 5-byte prefix: {} queries x {} table entries -> {} hits; prediction from measured H2(5)={} bits: {:.0}; naive 2^-40 (5 x log2 255) prediction: {:.0}\n",
            if appended { "even/appended" } else { "odd" }, n, n, hits, h2, pred, ideal);
    }
    out
}
