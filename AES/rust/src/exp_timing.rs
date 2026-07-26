// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! E10: wall-clock cost of one 255-element table entry (scalar Rust,
//! single core) for the I_{m,n}/packed-POW path and the chi-star path.
use crate::bridge::*;
use crate::gf::Gf;
use crate::rng::Rng;
use std::time::Instant;

pub fn e10(gf: &Gf) -> String {
    let pt = PowTab::new(gf);
    let mut rng = Rng::new(0xE10);
    let iters = 200_000usize;
    // pre-generate delta-sets of 255 values (as if peeled: v_w bytes) and reference
    let mut sets = vec![[0u8; 255]; 1024];
    for s in sets.iter_mut() { for x in s.iter_mut() { *x = rng.byte(); } }
    let mut invlinv = [0u8; 256];
    for x in 0..256usize { invlinv[x] = gf.inv[gf.linv[x] as usize]; }
    // (a) I_{m,n}: per element 1 lookup (invlinv) + 64-byte pow-row XOR; combine per entry
    let t0 = Instant::now();
    let mut sink = 0u8;
    for it in 0..iters {
        let vals = &sets[it & 1023];
        let w0 = vals[0];
        let mut acc = [0u64; 8];
        for i in 0..255 {
            let g = invlinv[(w0 ^ vals[i]) as usize];
            let row = &pt.row[g as usize];
            for w in 0..8 { acc[w] ^= u64::from_le_bytes(row[8 * w..8 * w + 8].try_into().unwrap()); }
        }
        let mut p = [0u8; 128];
        p[0] = 1;
        for j in 0..64 { let bytes = acc[j / 8].to_le_bytes(); p[2 * j + 1] = bytes[j % 8]; }
        for r in (2..128).step_by(2) { p[r] = gf.mul(p[r / 2], p[r / 2]); }
        let pms = all_pms(gf, &p);
        if let Some(f) = imn_fingerprint(gf, &pms) { sink ^= f.bytes[0]; }
    }
    let dt_imn = t0.elapsed();
    // (b) chi-star: per element chi bit flip + 4-byte odd-power accumulate; canonicalize per entry
    let t1 = Instant::now();
    for it in 0..iters {
        let vals = &sets[it & 1023];
        let w0 = vals[0];
        let mut vv = [0u8; 255];
        for i in 0..255 { vv[i] = invlinv[(w0 ^ vals[i]) as usize]; }
        if let Some(c) = chi_star(gf, &pt, &vv, 0) { sink ^= c.vec[0]; }
    }
    let dt_chi = t1.elapsed();
    // (c) just the chi accumulate (bit flips + S_1,S_3,S_5,S_7 via a 4-byte pow row) without canonicalization
    let mut pow4 = [[0u8; 4]; 256];
    for v in 0..256usize { pow4[v] = [v as u8, gf.pow(v as u8, 3), gf.pow(v as u8, 5), gf.pow(v as u8, 7)]; }
    let t2 = Instant::now();
    for it in 0..iters {
        let vals = &sets[it & 1023];
        let w0 = vals[0];
        let mut chi = [0u64; 4];
        let mut sacc = 0u32;
        for i in 0..255 {
            let g = invlinv[(w0 ^ vals[i]) as usize];
            chi[(g >> 6) as usize] ^= 1u64 << (g & 63);
            sacc ^= u32::from_le_bytes(pow4[g as usize]);
        }
        sink ^= chi[0] as u8 ^ (sacc as u8);
    }
    let dt_chi_acc = t2.elapsed();
    let ns = |d: std::time::Duration| d.as_nanos() as f64 / iters as f64;
    // rough clock estimate
    format!(
        "E10 timing (single core, scalar Rust, entries of 255 elements; sink={sink}):\n\
         \x20 I_(m,n) packed-POW path (1 lookup + 64-byte row XOR per element, 15 P_m + 12 ratios per entry): {:.0} ns/entry = {:.2} ns/element\n\
         \x20 chi accumulate only (bit flip + 4-byte power row per element): {:.0} ns/entry = {:.2} ns/element\n\
         \x20 full chi-star canonicalization (as implemented here, incl. a scalar 256-step permutation): {:.0} ns/entry\n\
         \x20 (paper: 83 cycles/entry claimed for a GFNI+AVX-512 kernel; at 3 GHz that is 27.7 ns/entry, i.e. 0.11 ns per element for 255 lookups + 5 wide XORs + solve + permutation)\n",
        ns(dt_imn), ns(dt_imn) / 255.0, ns(dt_chi_acc), ns(dt_chi_acc) / 255.0, ns(dt_chi))
}
