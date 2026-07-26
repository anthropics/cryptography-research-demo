// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! SR(7,4,4,4): 7-round small-scale AES with a real 4x4 state and 4-bit
//! words (Cid-Murphy-Robshaw style: inversion in GF(2^4) + GF(2)-affine
//! S-box, ShiftRows, MixColumns circulant (2,3,1,1) over GF(16)).
//!
//! E13a: cipher self-tests + constructive right pairs + differential shape.
//! E13b: honest-bridge statistics (parity trick, Alg.1 rates) at 4-bit words.
//! E13c: end-to-end -- one genuine right pair, the online delta-set with the
//!       true 8 key nibbles, and a STREAMING pass over the FULL basic
//!       offline table (all outer parameter tuples x all DDT branch choices):
//!       count entries, count fingerprint collisions with the true query and
//!       with wrong-key queries, run the multiset-consistency dispatch on
//!       every collision, and check the true entry is present and survives.
use crate::rng::Rng;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

// ---------------------------------------------------------------- GF(16)
pub struct Gf16 {
    pub log: [u8; 16], pub exp: [u8; 32], pub inv: [u8; 16],
    pub sbox: [u8; 16], pub inv_sbox: [u8; 16], pub linv: [u8; 16],
    pub mc: [[u8; 4]; 4], pub imc: [[u8; 4]; 4],
    pub mtab: [[u8; 16]; 16], pub pow357: [[u8; 3]; 16],
}
fn xt(a: u8) -> u8 { let hi = a & 0x8; let mut r = (a << 1) & 0xf; if hi != 0 { r ^= 0x3 } r }
fn mul_slow(mut a: u8, mut b: u8) -> u8 { let mut p = 0; while b != 0 { if b & 1 != 0 { p ^= a } a = xt(a); b >>= 1; } p }
impl Gf16 {
    pub fn new() -> Gf16 {
        let mut log = [0u8; 16]; let mut exp = [0u8; 32];
        let mut x = 1u8;
        for i in 0..15 { exp[i] = x; log[x as usize] = i as u8; x = mul_slow(x, 2); }
        assert_eq!(x, 1);
        for i in 15..32 { exp[i] = exp[i - 15]; }
        let mut inv = [0u8; 16];
        for a in 1..16 { inv[a] = exp[15 - log[a] as usize]; }
        // affine part: circulant of (1,1,1,0) over GF(2), constant 6
        let rows = [0b0111u8, 0b1110, 0b1101, 0b1011];
        let mut ltab = [0u8; 16];
        for y in 0..16u8 {
            let mut r = 0u8;
            for i in 0..4 { let bit = ((rows[i] & y).count_ones() & 1) as u8; r |= bit << i; }
            ltab[y as usize] = r;
        }
        let mut linv = [0u8; 16];
        for y in 0..16 { linv[ltab[y] as usize] = y as u8; }
        for y in 1..16u8 { assert_ne!(ltab[y as usize], 0, "L not injective"); }
        let mut sbox = [0u8; 16]; let mut inv_sbox = [0u8; 16];
        for a in 0..16 { sbox[a] = ltab[inv[a] as usize] ^ 0x6; }
        for a in 0..16 { inv_sbox[sbox[a] as usize] = a as u8; }
        let mc = [[2u8, 3, 1, 1], [1, 2, 3, 1], [1, 1, 2, 3], [3, 1, 1, 2]];
        // invert mc over GF(16) by Gauss-Jordan
        let g = |a: u8, b: u8| -> u8 { if a == 0 || b == 0 { 0 } else { exp[log[a as usize] as usize + log[b as usize] as usize] } };
        let mut m = [[0u8; 8]; 4];
        for i in 0..4 { for j in 0..4 { m[i][j] = mc[i][j]; } m[i][4 + i] = 1; }
        for c in 0..4 {
            let mut p = c; while m[p][c] == 0 { p += 1; }
            m.swap(c, p);
            let iv = inv[m[c][c] as usize];
            for j in 0..8 { m[c][j] = g(m[c][j], iv); }
            for r in 0..4 { if r != c && m[r][c] != 0 { let f = m[r][c]; for j in 0..8 { m[r][j] ^= g(f, m[c][j]); } } }
        }
        let mut imc = [[0u8; 4]; 4];
        for i in 0..4 { for j in 0..4 { imc[i][j] = m[i][4 + j]; } }
        let mut mtab = [[0u8; 16]; 16];
        for a in 0..16u8 { for b in 0..16u8 { mtab[a as usize][b as usize] = mul_slow(a, b); } }
        let mut pow357 = [[0u8; 3]; 16];
        for v in 0..16u8 { pow357[v as usize] = [mul_slow(mul_slow(v, v), v), mul_slow(mul_slow(mul_slow(v, v), mul_slow(v, v)), v), mul_slow(mul_slow(mul_slow(v, v), mul_slow(v, v)), mul_slow(mul_slow(v, v), v))]; }
        Gf16 { log, exp, inv, sbox, inv_sbox, linv, mc, imc, mtab, pow357 }
    }
    #[inline] pub fn mul(&self, a: u8, b: u8) -> u8 { self.mtab[a as usize][b as usize] }
    #[inline] pub fn div(&self, a: u8, b: u8) -> u8 { if a == 0 { 0 } else { self.exp[((self.log[a as usize] as i32 - self.log[b as usize] as i32).rem_euclid(15)) as usize] } }
    #[inline] pub fn pow(&self, a: u8, e: u32) -> u8 { if e == 0 { 1 } else if a == 0 { 0 } else { self.exp[((self.log[a as usize] as u32 * (e % 15)) % 15) as usize] } }
}

// ---------------------------------------------------------------- cipher
pub type St = [u8; 16];
pub const ROUNDS: usize = 7;
fn sr(s: &St) -> St { let mut r = [0u8; 16]; for c in 0..4 { for row in 0..4 { r[4 * c + row] = s[4 * ((c + row) % 4) + row] } } r }
fn isr(s: &St) -> St { let mut r = [0u8; 16]; for c in 0..4 { for row in 0..4 { r[4 * ((c + row) % 4) + row] = s[4 * c + row] } } r }
fn mc(f: &Gf16, s: &St) -> St { let mut r = [0u8; 16]; for c in 0..4 { for i in 0..4 { let mut v = 0; for j in 0..4 { v ^= f.mul(f.mc[i][j], s[4 * c + j]); } r[4 * c + i] = v; } } r }
fn imc(f: &Gf16, s: &St) -> St { let mut r = [0u8; 16]; for c in 0..4 { for i in 0..4 { let mut v = 0; for j in 0..4 { v ^= f.mul(f.imc[i][j], s[4 * c + j]); } r[4 * c + i] = v; } } r }
fn sb(f: &Gf16, s: &St) -> St { let mut r = *s; for x in r.iter_mut() { *x = f.sbox[*x as usize] } r }
fn isb(f: &Gf16, s: &St) -> St { let mut r = *s; for x in r.iter_mut() { *x = f.inv_sbox[*x as usize] } r }
fn xor(a: &St, b: &St) -> St { let mut r = [0u8; 16]; for i in 0..16 { r[i] = a[i] ^ b[i] } r }

pub fn expand_key(f: &Gf16, key: &St) -> Vec<St> {
    let mut w = vec![[0u8; 4]; 4 * (ROUNDS + 2)];
    for i in 0..4 { for j in 0..4 { w[i][j] = key[4 * i + j]; } }
    let mut rc = 1u8;
    for i in 4..w.len() {
        let mut t = w[i - 1];
        if i % 4 == 0 {
            t = [f.sbox[t[1] as usize] ^ rc, f.sbox[t[2] as usize], f.sbox[t[3] as usize], f.sbox[t[0] as usize]];
            rc = f.mul(rc, 2);
        }
        for j in 0..4 { w[i][j] = w[i - 4][j] ^ t[j]; }
    }
    (0..ROUNDS + 2).map(|r| { let mut k = [0u8; 16]; for c in 0..4 { for j in 0..4 { k[4 * c + j] = w[4 * r + c][j]; } } k }).collect()
}
pub struct Trace { pub xs: Vec<St>, pub c: St }
pub fn enc(f: &Gf16, rk: &[St], p: &St) -> Trace {
    let mut xs = vec![]; let mut x = xor(p, &rk[0]); let mut c = [0u8; 16];
    for i in 0..ROUNDS {
        xs.push(x);
        let z = sr(&sb(f, &x));
        if i + 1 < ROUNDS { x = xor(&mc(f, &z), &rk[i + 1]); } else { c = xor(&z, &rk[ROUNDS]); }
    }
    Trace { xs, c }
}
fn plaintext_from_x1(f: &Gf16, rk: &[St], x1: &St) -> St {
    let x0 = isb(f, &isr(&imc(f, &xor(x1, &rk[1]))));
    xor(&x0, &rk[0])
}

// ---------------------------------------------------------------- power sums / chi over GF(16)
fn power_sums16_fast(f: &Gf16, vals: &[u8], extra_zero: bool) -> [u8; 16] {
    // only p[0..7]: p1,p3,p5,p7 accumulated; even ones by squaring
    let mut p = [0u8; 16];
    p[0] = ((vals.len() + extra_zero as usize) & 1) as u8;
    let (mut a1, mut a3, mut a5, mut a7) = (0u8, 0u8, 0u8, 0u8);
    for &v in vals { let t = f.pow357[v as usize]; a1 ^= v; a3 ^= t[0]; a5 ^= t[1]; a7 ^= t[2]; }
    p[1] = a1; p[3] = a3; p[5] = a5; p[7] = a7;
    p[2] = f.mul(a1, a1); p[4] = f.mul(p[2], p[2]); p[6] = f.mul(a3, a3);
    p
}
fn power_sums16(f: &Gf16, vals: &[u8], extra_zero: bool) -> [u8; 16] {
    // p[r] for r = 0..14 (p[0] = |V| mod 2), plus we keep p[8..14]
    let mut p = [0u8; 16];
    p[0] = ((vals.len() + extra_zero as usize) & 1) as u8;
    for r in 1..15u32 { let mut a = 0u8; for &v in vals { a ^= f.pow(v, r); } p[r as usize] = a; }
    p
}
fn pm16(f: &Gf16, m: u8, p: &[u8; 16]) -> u8 {
    // pairwise sum P_m via submask expansion; exponents of value powers reduced mod 15 (v^15 = v^0=1 handled by pow)
    let m32 = m as u32; let mut acc = 0u8;
    let mut k = m32;
    loop {
        k = (k.wrapping_sub(1)) & m32;
        let comp = m32 ^ k;
        if k != 0 && k < comp { acc ^= f.mul(p[k as usize], p[comp as usize]); }
        if k == 0 { break; }
    }
    if p[0] & 1 == 1 { acc ^= p[m as usize]; }
    if m32.count_ones() == 1 { acc ^= p[m as usize]; }
    acc
}
/// canonical chi-star over GF(16): if P_7 != 0 use the (alpha*, beta*) frame
/// (odd n) with brute-min over beta if S1 is unusable; else brute-min over
/// the whole AGL(1,16) (240 elements). Both are class functions.
fn chi_vec(vals: &[u8], extra_zero: bool) -> u16 {
    let mut c = 0u16; for &v in vals { c ^= 1u16 << v; }
    if extra_zero { c ^= 1; }
    c
}
fn permute16(f: &Gf16, chi: u16, a: u8, b: u8) -> u16 {
    // out[a d + b] = chi[d]
    let mut out = 0u16;
    let ainv = f.inv[a as usize] as usize;
    let row = &f.mtab[ainv];
    for t in 0..16u8 { let d = row[(t ^ b) as usize]; out |= ((chi >> d) & 1) << t; }
    out
}
/// frame-based chi-star only (odd |V|); None when P_7 = 0 (frame undefined)
pub fn chi_star16_frame(f: &Gf16, vals: &[u8]) -> Option<(u16, u8)> {
    let p = power_sums16_fast(f, vals, false);
    let p7 = pm16(f, 7, &p);
    if p7 == 0 { return None; }
    let alpha = f.pow(f.inv[p7 as usize], 13);
    let beta = f.mul(alpha, p[1]);
    Some((permute16(f, chi_vec(vals, false), alpha, beta), p7))
}
pub fn chi_star16(f: &Gf16, vals: &[u8], extra_zero: bool) -> u16 {
    let p = power_sums16_fast(f, vals, extra_zero);
    let chi = chi_vec(vals, extra_zero);
    let n_odd = p[0] & 1 == 1;
    let p7 = pm16(f, 7, &p);
    if p7 != 0 && n_odd {
        // alpha* = P_7^{-1/7}; 7^{-1} mod 15 = 13
        let alpha = f.pow(f.inv[p7 as usize], 13);
        let beta = f.mul(alpha, p[1]);
        return permute16(f, chi, alpha, beta);
    }
    // brute-min over the whole group (class function)
    let mut best = u16::MAX;
    for a in 1..16u8 { for b in 0..16u8 { let v = permute16(f, chi, a, b); if v < best { best = v } } }
    best
}

// ---------------------------------------------------------------- keyless offline d
fn keyless_d16(f: &Gf16, x2c0: &[u8; 4], x3: &St, x4d: &[u8; 4], e: u8) -> u8 {
    let dx2 = [f.mul(f.mc[0][0], e), f.mul(f.mc[1][0], e), f.mul(f.mc[2][0], e), f.mul(f.mc[3][0], e)];
    let mut dy2 = [0u8; 16];
    for j in 0..4 { dy2[j] = f.sbox[x2c0[j] as usize] ^ f.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
    let dx3 = mc(f, &sr(&dy2));
    let mut dy3 = [0u8; 16];
    for j in 0..16 { dy3[j] = f.sbox[x3[j] as usize] ^ f.sbox[(x3[j] ^ dx3[j]) as usize]; }
    let dx4 = mc(f, &sr(&dy3));
    let diag = [0usize, 5, 10, 15];
    let mut d = 0u8;
    for r in 0..4 {
        let dy = f.sbox[x4d[r] as usize] ^ f.sbox[(x4d[r] ^ dx4[diag[r]]) as usize];
        d ^= f.mul(f.mc[0][r], dy); // x5[0] = MC row 0 of (y4 diagonal as z4 col 0)
    }
    d
}

// ---------------------------------------------------------------- E13a: sanity + right pairs
pub struct Pair { pub key: St, pub p0: St, pub p1: St }

/// Construct a genuine right pair for the 4-round core inside 7-round SR
/// under key `key`, using key knowledge (meet in the middle at round 3).
pub fn make_right_pair(f: &Gf16, rng: &mut Rng, key: &St) -> (Pair, u64) {
    let base = rng.next_u64();
    (0..u64::MAX).into_par_iter().find_map_any(|i| {
        let mut r = Rng::new(base ^ i.wrapping_mul(0x9e37_79b9));
        make_right_pair_seq(f, &mut r, key, 200_000)
    }).unwrap()
}
pub fn make_right_pair_seq(f: &Gf16, rng: &mut Rng, key: &St, max_attempts: u64) -> Option<(Pair, u64)> {
    let rk = expand_key(f, key);
    let mut ddt: Vec<Vec<Vec<u8>>> = vec![vec![Vec::new(); 16]; 16];
    for din in 1..16usize { for x in 0..16u8 { let dout = f.sbox[x as usize] ^ f.sbox[(x as usize) ^ din]; ddt[din][dout as usize].push(x); } }
    let mut attempts = 0u64;
    loop {
        attempts += 1;
        if attempts > max_attempts { return None; }
        let e = 1 + (rng.byte() % 15) as u8;
        let mut dvec = [0u8; 4]; for r in 0..4 { dvec[r] = 1 + (rng.byte() % 15) as u8; }
        // forward: dx2 col0, choose dy2 col0 among reachable output diffs
        let dx2 = [f.mul(f.mc[0][0], e), f.mul(f.mc[1][0], e), f.mul(f.mc[2][0], e), f.mul(f.mc[3][0], e)];
        let mut dy2 = [0u8; 16]; let mut ok = true;
        for j in 0..4 {
            let opts: Vec<u8> = (1..16u8).filter(|&g| !ddt[dx2[j] as usize][g as usize].is_empty()).collect();
            if opts.is_empty() { ok = false; break; }
            dy2[j] = opts[(rng.byte() as usize) % opts.len()];
        }
        if !ok { continue; }
        let dx3 = mc(f, &sr(&dy2));
        // backward: dx4 = dvec on the diagonal
        let mut dx4 = [0u8; 16]; let diag = [0usize, 5, 10, 15];
        for r in 0..4 { dx4[diag[r]] = dvec[r]; }
        let dy3 = isr(&imc(f, &dx4));
        // round-3 DDT compat
        if (0..16).any(|j| dx3[j] == 0 || dy3[j] == 0 || ddt[dx3[j] as usize][dy3[j] as usize].is_empty()) { continue; }
        // per column c of x3: search x2 nibbles feeding it. z2 col c has y2 positions p_r = 4*((c+r)%4)+r for row r.
        // x3 col c = MC(z2 col c) ^ rk[3] col c ; x2 col0 nibble in this column: row r0 = (4-c)%4 (position r0, col 0)
        let mut x2 = [0u8; 16]; let mut x3 = [0u8; 16]; let mut fail = false;
        for c in 0..4 {
            let r0 = (4 - c) % 4; // row whose z2 source is x2 column 0
            let mut found = None;
            // enumerate: x2 col0 nibble at row r0 from round-2 DDT solutions; other 3 rows: y2 free (16 each)
            let sols_col0 = &ddt[dx2[r0] as usize][dy2[r0] as usize];
            'outer: for &x2v in sols_col0.iter() {
                let y_r0 = f.sbox[x2v as usize];
                for a in 0..16u8 { for b in 0..16u8 { for cc in 0..16u8 {
                    let mut zcol = [0u8; 4]; let mut it = [a, b, cc].into_iter();
                    for row in 0..4 { zcol[row] = if row == r0 { y_r0 } else { it.next().unwrap() }; }
                    // x3 col c = MC * zcol ^ rk[3] col c ; must lie in DDT solution sets of round 3
                    let mut good = true; let mut col = [0u8; 4];
                    for i in 0..4 {
                        let mut v = 0u8; for j in 0..4 { v ^= f.mul(f.mc[i][j], zcol[j]); }
                        v ^= rk[3][4 * c + i];
                        let pos = 4 * c + i;
                        if !ddt[dx3[pos] as usize][dy3[pos] as usize].contains(&v) { good = false; break; }
                        col[i] = v;
                    }
                    if good { found = Some((x2v, zcol, col)); break 'outer; }
                }}}
            }
            match found {
                None => { fail = true; break; }
                Some((x2v, zcol, col)) => {
                    x2[r0] = x2v; // column 0, row r0
                    for row in 0..4 { let src = 4 * ((c + row) % 4) + row; if row != r0 { x2[src] = f.inv_sbox[zcol[row] as usize]; } }
                    for i in 0..4 { x3[4 * c + i] = col[i]; }
                }
            }
        }
        if fail { continue; }
        // consistency check of the derived states
        let x3_chk = xor(&mc(f, &sr(&sb(f, &x2))), &rk[3]);
        if x3_chk != x3 { continue; }
        let mut x2p = x2; for j in 0..4 { x2p[j] ^= dx2[j]; }
        let mut x3p = x3; for j in 0..16 { x3p[j] ^= dx3[j]; }
        let x4 = xor(&mc(f, &sr(&sb(f, &x3))), &rk[4]);
        let x4p = xor(&mc(f, &sr(&sb(f, &x3p))), &rk[4]);
        let dx4c = xor(&x4, &x4p);
        if dx4c != dx4 { continue; }
        // round 4: need dw4 single active nibble at position 0
        let dw4 = xor(&mc(f, &sr(&sb(f, &x4))), &mc(f, &sr(&sb(f, &x4p))));
        if dw4[0] == 0 || (1..16).any(|j| dw4[j] != 0) { continue; }
        // x1 from x2 (invert round 1 with rk[2]... x2 = MC(SR(SB(x1))) ^ rk[2])
        let x1 = isb(f, &isr(&imc(f, &xor(&x2, &rk[2]))));
        let x1p = isb(f, &isr(&imc(f, &xor(&x2p, &rk[2]))));
        let dx1 = xor(&x1, &x1p);
        if dx1[0] == 0 || (1..16).any(|j| dx1[j] != 0) { continue; }
        let p0 = plaintext_from_x1(f, &rk, &x1);
        let p1 = plaintext_from_x1(f, &rk, &x1p);
        return Some((Pair { key: *key, p0, p1 }, attempts));
    }
}

pub fn e13a(f: &Gf16) -> String {
    let mut out = String::new();
    // S-box facts
    let mut ddt = [[0u32; 16]; 16];
    for d in 1..16usize { for x in 0..16usize { ddt[d][(f.sbox[x] ^ f.sbox[x ^ d]) as usize] += 1; } }
    let mut spec = std::collections::BTreeMap::new();
    for d in 1..16usize { for g in 0..16usize { *spec.entry(ddt[d][g]).or_insert(0u32) += 1; } }
    out += &format!("E13a SR(7,4,4,4): S-box = L(x^-1)+6 over GF(16), sbox = {:x?}\n  DDT value spectrum over 15x16 cells: {:?}\n", f.sbox, spec);
    out += &format!("  MC = {:?}, MC^-1 = {:?}\n", f.mc, f.imc);
    // cipher round trip: encrypt then check trace lengths / decryption via reconstruction not needed
    let mut rng = Rng::new(0x1a);
    let key = rng.state16().map(|b| b & 0xf); let key: St = key;
    let rk = expand_key(f, &key);
    // right pairs: build a few, verify plaintext/ciphertext difference shape
    let (mut n_pairs, mut shape_ok, mut tot_attempts) = (0, 0, 0u64);
    for i in 0..3u64 {
        let mut r2 = Rng::new(0x100 + i);
        let (pr, att) = make_right_pair(f, &mut r2, &key);
        tot_attempts += att; n_pairs += 1;
        let dp = xor(&pr.p0, &pr.p1);
        let t0 = enc(f, &rk, &pr.p0); let t1 = enc(f, &rk, &pr.p1);
        let dc = xor(&t0.c, &t1.c);
        let dx5 = xor(&t0.xs[5], &t1.xs[5]);
        let dx1 = xor(&t0.xs[1], &t1.xs[1]);
        let diag = [0usize, 5, 10, 15]; let anti = [0usize, 7, 10, 13];
        let ok = (0..16).all(|j| (dp[j] != 0) == diag.contains(&j))
            && (0..16).all(|j| (dc[j] != 0) == anti.contains(&j))
            && dx5[0] != 0 && (1..16).all(|j| dx5[j] == 0)
            && dx1[0] != 0 && (1..16).all(|j| dx1[j] == 0);
        if ok { shape_ok += 1; }
    }
    out += &format!("  right pairs constructed: {n_pairs}/3, differential shape (dP diagonal / dx1 single / dx5 single / dC anti-diagonal): {shape_ok}/3, mean construction attempts {:.0}\n", tot_attempts as f64 / n_pairs as f64);
    out
}

// ---------------------------------------------------------------- E13b: honest bridge stats at 4 bits
pub fn e13b(f: &Gf16, instances: u64) -> String {
    let nthreads = rayon::current_num_threads() as u64;
    let per = instances / nthreads;
    let (n, s0, raw_ok, app_ok, both, neither, alg_raw, alg_app, alg_either, corr, key_either, undef_frame) =
        (0..nthreads).into_par_iter().map(|t| {
            let mut rng = Rng::new(0xB13 + t);
            let (mut n, mut s0, mut ro, mut ao, mut bo, mut ne, mut ar, mut aa, mut ae, mut co, mut ke, mut uf) = (0u64,0,0,0,0,0,0,0,0,0,0,0);
            for _ in 0..per {
                let key: St = rng.state16().map(|b| b & 0xf);
                let rk = expand_key(f, &key);
                let base = rng.state16().map(|b| b & 0xf);
                let k6 = rk[7]; let pos = [0usize, 13, 10, 7];
                let u5 = imc(f, &rk[6]);
                let mut a = [0u8; 16]; let mut w = [0u8; 16];
                for v in 0..16u8 {
                    let mut x1 = base; x1[0] = base[0] ^ v;
                    let p = plaintext_from_x1(f, &rk, &x1);
                    let tr = enc(f, &rk, &p);
                    a[v as usize] = tr.xs[5][0];
                    let mut x6c0 = [0u8; 4];
                    for j in 0..4 { x6c0[j] = f.inv_sbox[(tr.c[pos[j]] ^ k6[pos[j]]) as usize]; }
                    let mut vv = 0u8; for j in 0..4 { vv ^= f.mul(f.imc[0][j], x6c0[j]); }
                    w[v as usize] = vv;
                    debug_assert_eq!(f.inv_sbox[(vv ^ u5[0]) as usize], a[v as usize]);
                }
                n += 1;
                let s = a[0]; if s == 0 { s0 += 1; }
                let g: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(w[0] ^ w[v]) as usize] as usize]).collect();
                let dinv: Vec<u8> = (1..16).map(|v| f.inv[(a[0] ^ a[v]) as usize]).collect();
                let mut ca = 0; let mut cb = 0;
                for v in 1..16 { if a[v] == 0 { ca += 1 } if a[v] == a[0] { cb += 1 } }
                let b = ca + cb;
                // parity claim on P_m with the two usable exponents 7 and 11 (11 ~ 7^8)
                let pon = power_sums16(f, &g, false); let poff = power_sums16(f, &dinv, false);
                let mut pona = pon; pona[0] ^= 1; let mut poffa = poff; poffa[0] ^= 1;
                let mut r_ok = true; let mut a_ok = true;
                for &m in &[7u8, 11, 13, 14] {
                    let s2m = f.pow(s, 2 * m as u32);
                    if pm16(f, m, &pon) != f.mul(s2m, pm16(f, m, &poff)) { r_ok = false; }
                    if pm16(f, m, &pona) != f.mul(s2m, pm16(f, m, &poffa)) { a_ok = false; }
                }
                if r_ok { ro += 1 } if a_ok { ao += 1 } if r_ok && a_ok { bo += 1 } if !r_ok && !a_ok { ne += 1 }
                let _ = b;
                // Alg 1 exact multiset test
                let mut pred: Vec<u8> = dinv.iter().map(|&d| f.mul(f.mul(s, s), d) ^ s).collect();
                let mut gs = g.clone(); gs.sort_unstable(); pred.sort_unstable();
                let raw_pass = gs == pred;
                let mut ga = g.clone(); ga.push(0); ga.sort_unstable();
                let mut pa = pred.clone(); pa.push(s); pa.sort_unstable();
                let app_pass = ga == pa;
                if raw_pass { ar += 1 } if app_pass { aa += 1 } if raw_pass || app_pass { ae += 1 }
                let gc: Vec<u8> = { let mut x: Vec<u8> = g.iter().cloned().filter(|&x| x != 0 && x != s).collect(); x.sort_unstable(); x };
                let pc: Vec<u8> = { let mut x: Vec<u8> = pred.iter().cloned().filter(|&x| x != 0 && x != s).collect(); x.sort_unstable(); x };
                if gc == pc { co += 1 }
                // chi-star agreement raw or appended
                let con = chi_star16(f, &g, false); let cof = chi_star16(f, &dinv, false);
                let cona = chi_star16(f, &g, true); let cofa = chi_star16(f, &dinv, true);
                if con == cof || cona == cofa { ke += 1 }
                if pm16(f, 7, &pon) == 0 { uf += 1 }
            }
            (n, s0, ro, ao, bo, ne, ar, aa, ae, co, ke, uf)
        }).reduce(|| (0,0,0,0,0,0,0,0,0,0,0,0), |x, y| (x.0+y.0,x.1+y.1,x.2+y.2,x.3+y.3,x.4+y.4,x.5+y.5,x.6+y.6,x.7+y.7,x.8+y.8,x.9+y.9,x.10+y.10,x.11+y.11));
    let pc = |x: u64| 100.0 * x as f64 / n as f64;
    format!("E13b honest bridge on SR(7,4,4,4), {} instances (delta-sets of 16):\n  s = 0: {:.3}% (expect 1/16 = 6.25%)\n  P_m parity claim (m in coset {{7,11,13,14}}): raw matches {:.3}%, appended {:.3}%, both {:.3}%, neither {:.3}%\n  Alg.1 exact-multiset test on the true hit: raw {:.2}%, appended {:.2}%, either {:.2}%; corrected (delete {{0,s}}) {:.3}%\n  chi-star (16-bit canonical vector) online == offline (raw or appended): {:.3}%\n  P_7 = 0 online (frame undefined -> brute-min branch): {:.3}% (expect ~1/16)\n",
        n, pc(s0), pc(raw_ok), pc(app_ok), pc(both), pc(neither), pc(alg_raw), pc(alg_app), pc(alg_either), pc(corr), pc(key_either), pc(undef_frame))
}

// ---------------------------------------------------------------- key schedule inversion
/// Given round key rk[idx] (idx = word-block index), recover rk[0] = master key.
pub fn invert_key_schedule(f: &Gf16, rk_last: &St, idx: usize) -> St {
    let nw = 4 * (idx + 1);
    let mut w = vec![[0u8; 4]; nw];
    for c in 0..4 { for j in 0..4 { w[nw - 4 + c][j] = rk_last[4 * c + j]; } }
    // rcon used at word i (i%4==0) is 2^(i/4-1)
    for i in (4..nw).rev() {
        let mut t = w[i - 1];
        if i % 4 == 0 {
            let rc = f.pow(2, (i / 4 - 1) as u32);
            t = [f.sbox[t[1] as usize] ^ rc, f.sbox[t[2] as usize], f.sbox[t[3] as usize], f.sbox[t[0] as usize]];
        }
        for j in 0..4 { w[i - 4][j] = w[i][j] ^ t[j]; }
    }
    let mut k = [0u8; 16]; for c in 0..4 { for j in 0..4 { k[4 * c + j] = w[c][j]; } } k
}

// ---------------------------------------------------------------- E13c: FULL KEY RECOVERY
/// The attack sees: a chosen-plaintext oracle for 7-round SR(4,4,4) under an
/// unknown key, and one right pair (P0,P1) (planted; a genuine right pair
/// found by the harness with key knowledge, since the 2^57-data collection
/// phase is out of reach). It streams the FULL basic offline table.
/// per-key result of the honest sampled-coverage experiment
pub struct KeyRes {
    pub attempts_until_pair: u64, pub rej_s0: u64, pub rej_bodd: u64, pub rej_frame: u64,
    pub n_kc: usize, pub n_uc: usize, pub n_query: usize, pub chosen_pt: u64,
    pub true_entries: u64, pub true_hits: u64, pub true_surv: u64, pub true_completions: u64, pub true_secs: f64,
    pub samp: Vec<(u64, u64, u64)>, pub samp_secs: f64,
    pub recovered: bool, pub recovered_matches: bool,
    pub alg1_exact: i8, pub alg1_corrected: i8, pub s_true: u8,
}
/// E13e: honest downstream pipeline given an oracle-supplied right pair, with
/// SAMPLED table coverage (true outer task fully enumerated with completion;
/// K random other outer tasks fully enumerated for statistics). The attack
/// path never reads the secret key; only the harness's choice of which outer
/// tasks to stream (and the post-hoc labels) uses it.
pub fn e13e(f: &Gf16, seed: u64, ktasks: usize) -> (String, KeyRes) {
    let mut out = String::new();
    let mut rng = Rng::new(0xC13 ^ seed.wrapping_mul(0x9e37));
    let secret_key: St = rng.state16().map(|b| b & 0xf);
    let rk_secret = expand_key(f, &secret_key);
    let oracle = |p: &St| -> St { enc(f, &rk_secret, p).c };
    // ---- oracle-supplied right pair (disclosed shortcut standing in for the ~2^57-CP harvest);
    //      conditioning (s != 0, even bad-slot count, defined chi* frame) is measured, not hidden
    let mut attempts_total = 0u64; let (mut rej_s0, mut rej_bodd, mut rej_frame) = (0u64, 0u64, 0u64);
    let (pair, s_true, e_true, f_true, anchor) = loop {
        let (pr, att) = make_right_pair(f, &mut rng, &secret_key);
        attempts_total += att;
        let t0 = enc(f, &rk_secret, &pr.p0); let t1 = enc(f, &rk_secret, &pr.p1);
        let mut a = [0u8; 16]; let mut w = [0u8; 16];
        let x1 = t0.xs[1];
        let posad = [0usize, 13, 10, 7];
        for v in 0..16u8 {
            let mut x = x1; x[0] ^= v;
            let tr = enc(f, &rk_secret, &plaintext_from_x1(f, &rk_secret, &x));
            a[v as usize] = tr.xs[5][0];
            let mut vv = 0u8;
            for j in 0..4 { vv ^= f.mul(f.imc[0][j], f.inv_sbox[(tr.c[posad[j]] ^ rk_secret[7][posad[j]]) as usize]); }
            w[v as usize] = vv;
        }
        let s = a[0]; let b = (1..16).filter(|&v| a[v] == 0 || a[v] == a[0]).count();
        if s == 0 { rej_s0 += 1; continue; }
        if b % 2 == 1 { rej_bodd += 1; continue; }
        let gt: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(w[0] ^ w[v]) as usize] as usize]).collect();
        if chi_star16_frame(f, &gt).is_none() { rej_frame += 1; continue; } // need P_7 != 0 (frame variant)
        let mut x2c0 = [0u8; 4]; for j in 0..4 { x2c0[j] = t0.xs[2][j]; }
        let dg = [0usize, 5, 10, 15]; let mut x4d = [0u8; 4]; for j in 0..4 { x4d[j] = t0.xs[4][dg[j]]; }
        let e_true = f.sbox[t0.xs[1][0] as usize] ^ f.sbox[t1.xs[1][0] as usize];
        let f_true = t0.xs[5][0] ^ t1.xs[5][0];
        break (pr, s, e_true, f_true, (x2c0, t0.xs[3], x4d));
    };
    eprintln!("  [milestone] right pair + instance ready");
    let (p0, p1) = (pair.p0, pair.p1);
    let (c0, c1) = (oracle(&p0), oracle(&p1));
    out += "E13e SR(7,4,4,4) honest downstream pipeline (oracle-supplied right pair, sampled table coverage)\n";
    out += &format!("  secret key = {:x?}\n  right pair P0={:x?} P1={:x?}\n  (truth for reference: s = {:x}, outer params e = {:x}, f = {:x})\n", secret_key, p0, p1, s_true, e_true, f_true);
    // ---- online candidate keys from the pair (attacker's view only)
    let dg = [0usize, 5, 10, 15]; let anti = [0usize, 13, 10, 7];
    let dp = xor(&p0, &p1); let dc = xor(&c0, &c1);
    let mut kcands: Vec<[u8; 4]> = vec![];
    for k in 0..(1u32 << 16) {
        let kc = [(k & 0xf) as u8, ((k >> 4) & 0xf) as u8, ((k >> 8) & 0xf) as u8, ((k >> 12) & 0xf) as u8];
        // z0 column 0 = SR(y0) col 0 = (y0[0], y0[5], y0[10], y0[15]) = S(P[dg]^kc)
        let mut d = [0u8; 4];
        for r in 0..4 { d[r] = f.sbox[(p0[dg[r]] ^ kc[r]) as usize] ^ f.sbox[(p1[dg[r]] ^ kc[r]) as usize]; }
        let mut dw = [0u8; 4]; for i in 0..4 { for j in 0..4 { dw[i] ^= f.mul(f.mc[i][j], d[j]); } }
        if dw[0] != 0 && dw[1] == 0 && dw[2] == 0 && dw[3] == 0 { kcands.push(kc); }
    }
    let mut ucands: Vec<[u8; 4]> = vec![];
    for k in 0..(1u32 << 16) {
        let uc = [(k & 0xf) as u8, ((k >> 4) & 0xf) as u8, ((k >> 8) & 0xf) as u8, ((k >> 12) & 0xf) as u8];
        let mut d = [0u8; 4];
        for r in 0..4 { d[r] = f.inv_sbox[(c0[anti[r]] ^ uc[r]) as usize] ^ f.inv_sbox[(c1[anti[r]] ^ uc[r]) as usize]; }
        let mut dz = [0u8; 4]; for i in 0..4 { for j in 0..4 { dz[i] ^= f.mul(f.imc[i][j], d[j]); } }
        if dz[0] != 0 && dz[1] == 0 && dz[2] == 0 && dz[3] == 0 { ucands.push(uc); }
    }
    let true_kc = [rk_secret[0][0], rk_secret[0][5], rk_secret[0][10], rk_secret[0][15]];
    let true_uc = [rk_secret[7][0], rk_secret[7][13], rk_secret[7][10], rk_secret[7][7]];
    out += &format!("  online: k_-1[0,5,10,15] candidates from the input differential = {} (contains true: {}); u_6[0,7,10,13] candidates = {} (contains true: {})\n",
        kcands.len(), kcands.contains(&true_kc), ucands.len(), ucands.contains(&true_uc));
    let _ = dp; let _ = dc;
    // ---- for each (kc, uc): build the delta-set (chosen plaintexts), get its ciphertexts, form G and its raw chi-star
    struct Query { kc: [u8; 4], uc: [u8; 4], w: [u8; 16], g: Vec<u8>, hist: [u8; 16], p7: u8 }
    let mut queries: Vec<Query> = vec![];
    let mut chosen_pt = 2u64;
    for kc in kcands.iter() {
        // delta-set plaintexts around P0 for this candidate k_-1 diagonal
        let mut cts = vec![[0u8; 16]; 16];
        for v in 0..16u8 {
            let mut m = p0;
            if v != 0 {
                // want dx1 col0 = (v,0,0,0): dz0 col0 = MC^-1 (v,0,0,0) -> dy0 at diagonal
                for r in 0..4 {
                    let dz = f.mul(f.imc[r][0], v);
                    let x0 = p0[dg[r]] ^ kc[r];
                    let y0 = f.sbox[x0 as usize] ^ dz;
                    m[dg[r]] = f.inv_sbox[y0 as usize] ^ kc[r];
                }
            }
            cts[v as usize] = oracle(&m);
        }
        chosen_pt += 15;
        for uc in ucands.iter() {
            let mut w = [0u8; 16];
            for v in 0..16usize {
                let mut vv = 0u8;
                for j in 0..4 { vv ^= f.mul(f.imc[0][j], f.inv_sbox[(cts[v][anti[j]] ^ uc[j]) as usize]); }
                w[v] = vv;
            }
            let g: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(w[0] ^ w[v]) as usize] as usize]).collect();
            let mut hist = [0u8; 16]; for &x in &g { hist[x as usize] += 1; }
            let p7q = pm16(f, 7, &power_sums16_fast(f, &g, false));
            queries.push(Query { kc: *kc, uc: *uc, w, g, hist, p7: p7q });
        }
    }
    let qfp: Vec<Option<(u16, u8)>> = queries.iter().map(|q| chi_star16_frame(f, &q.g)).collect();
    let mut qidx: Vec<Vec<u16>> = vec![vec![]; 65536];
    let mut qdead = 0usize;
    for (i, fp) in qfp.iter().enumerate() { match fp { Some((x, _)) => qidx[*x as usize].push(i as u16), None => qdead += 1 } }
    out += &format!("  {} online queries (raw parity, {} with undefined frame skipped), chosen plaintexts used = {} (all inside one 2^16 diagonal structure)\n", queries.len(), qdead, chosen_pt);
    eprintln!("  [milestone] {} queries built; starting stream", queries.len());
    // ---- STREAM the full basic offline table
    let ddt: Vec<Vec<Vec<u8>>> = { let mut d = vec![vec![Vec::new(); 16]; 16]; for din in 1..16usize { for x in 0..16u8 { let dout = f.sbox[x as usize] ^ f.sbox[(x as usize) ^ din]; d[din][dout as usize].push(x); } } d };
    let total_entries = AtomicU64::new(0);
    let dead_outer = AtomicU64::new(0);
    let fp_hits = AtomicU64::new(0);
    let true_seen = AtomicU64::new(0);
    let skipped_entries = AtomicU64::new(0);
    let done_tasks = AtomicU64::new(0);
    let survivors: Mutex<Vec<(usize, [u8; 4], St, [u8; 4], u8)>> = Mutex::new(vec![]); // (query idx, x2c0, x3, x4d, s_hat)
    let hist = Mutex::new(vec![0u64; 65536]);
    let (ax2c0, ax3, ax4d) = anchor;
    // ---- inline key completion for a survivor: DS tests at output rows 1..3 (2^20 each) + key-schedule inversion + trial encryption
    let complete = |qi: usize, x2c0: [u8; 4], x3: St, x4d: [u8; 4], _sh: u8| -> Option<St> {
        let q = &queries[qi];
        let mut cts = vec![[0u8; 16]; 16];
        for v in 0..16u8 {
            let mut m = p0;
            if v != 0 { for r in 0..4 { let dz = f.mul(f.imc[r][0], v); let x0 = p0[dg[r]] ^ q.kc[r]; let y0 = f.sbox[x0 as usize] ^ dz; m[dg[r]] = f.inv_sbox[y0 as usize] ^ q.kc[r]; } }
            cts[v as usize] = oracle(&m);
        }
        let mut k6 = [0u8; 16];
        for r in 0..4 { k6[anti[r]] = q.uc[r]; }
        for i in 1..4usize {
            let mut target: Vec<u8> = (1..16u8).map(|ep| {
                let dx2 = [f.mul(f.mc[0][0], ep), f.mul(f.mc[1][0], ep), f.mul(f.mc[2][0], ep), f.mul(f.mc[3][0], ep)];
                let mut dy2 = [0u8; 16]; for j in 0..4 { dy2[j] = f.sbox[x2c0[j] as usize] ^ f.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
                let dx3 = mc(f, &sr(&dy2));
                let mut dy3 = [0u8; 16]; for j in 0..16 { dy3[j] = f.sbox[x3[j] as usize] ^ f.sbox[(x3[j] ^ dx3[j]) as usize]; }
                let dx4 = mc(f, &sr(&dy3));
                let dgp = [0usize, 5, 10, 15]; let mut d = 0u8;
                for r in 0..4 { let dy = f.sbox[x4d[r] as usize] ^ f.sbox[(x4d[r] ^ dx4[dgp[r]]) as usize]; d ^= f.mul(f.mc[i][r], dy); }
                d
            }).collect();
            target.sort_unstable();
            let c = (4 - i) % 4;
            let poss: Vec<usize> = (0..4).map(|r| 4 * ((c + 4 - r) % 4) + r).collect();
            let mut hit: Option<[u8; 5]> = None; let mut nhit = 0;
            for cand in 0..(1u32 << 20) {
                let kk = [(cand & 0xf) as u8, ((cand >> 4) & 0xf) as u8, ((cand >> 8) & 0xf) as u8, ((cand >> 12) & 0xf) as u8, ((cand >> 16) & 0xf) as u8];
                let mut vals = [0u8; 16];
                for v in 0..16usize {
                    let mut z = 0u8;
                    for r in 0..4 { z ^= f.mul(f.imc[i][r], f.inv_sbox[(cts[v][poss[r]] ^ kk[r]) as usize]); }
                    vals[v] = f.inv_sbox[(z ^ kk[4]) as usize];
                }
                let mut m: Vec<u8> = (1..16).map(|v| vals[0] ^ vals[v]).collect();
                m.sort_unstable();
                if m == target { nhit += 1; if hit.is_none() { hit = Some(kk); } }
            }
            match hit { Some(kk) if nhit == 1 => { for r in 0..4 { k6[poss[r]] = kk[r]; } }, _ => return None }
        }
        let master = invert_key_schedule(f, &k6, 7);
        let rkm = expand_key(f, &master);
        if enc(f, &rkm, &p0).c == c0 && enc(f, &rkm, &p1).c == c1 { Some(master) } else { None }
    };
    let recovered_key: Mutex<Option<St>> = Mutex::new(None);
    let stop = std::sync::atomic::AtomicBool::new(false);
    // units: (e, ff, x2c0, is_true_task). True outer task fully (completion on); ktasks random others (completion off).
    let true_task = (e_true, f_true, ax2c0[0]);
    let all_tasks: Vec<(u8, u8, u8)> = (1..16u8).flat_map(|e| (1..16u8).flat_map(move |ff| (0..16u8).map(move |x0| (e, ff, x0)))).collect();
    let mut others: Vec<(u8, u8, u8)> = all_tasks.iter().copied().filter(|t| *t != true_task).collect();
    for i in 0..others.len() { let j = i + (rng.next_u64() as usize) % (others.len() - i); others.swap(i, j); }
    others.truncate(ktasks);
    let mut units: Vec<(u8, u8, [u8; 4], usize)> = vec![]; // task index 0 = true task, 1..=K = sampled
    for (ti, t) in std::iter::once(true_task).chain(others.iter().copied()).enumerate() {
        for x1n in 0..16u8 { for x2n in 0..16u8 { for x3n in 0..16u8 { units.push((t.0, t.1, [t.2, x1n, x2n, x3n], ti)); } } }
    }
    let ntasks = 1 + others.len();
    let task_entries: Vec<AtomicU64> = (0..ntasks).map(|_| AtomicU64::new(0)).collect();
    let task_hits: Vec<AtomicU64> = (0..ntasks).map(|_| AtomicU64::new(0)).collect();
    let task_surv: Vec<AtomicU64> = (0..ntasks).map(|_| AtomicU64::new(0)).collect();
    let completions = AtomicU64::new(0);
    let alg1_exact: AtomicU64 = AtomicU64::new(2); // 2 = not observed, 0/1 = fail/pass on the true (entry, query) pair
    let alg1_corrected: AtomicU64 = AtomicU64::new(2);
    out += &format!("  table coverage: true outer task (e,f,x0)={:?} fully enumerated with completion; {} random other outer tasks fully enumerated for statistics\n", true_task, others.len());
    let t_stream = std::time::Instant::now();
    units.par_iter().for_each(|&(e, ff, x2c0, ti)| {
        let is_true_task = ti == 0;
        let mut ltot = 0u64; let mut lfp = 0u64; let mut lts = 0u64; let mut lskip = 0u64;
        let mut lhist = vec![0u64; 65536];
        let mut lsurv: Vec<(usize, [u8; 4], St, [u8; 4], u8)> = vec![];
        let dgp = [0usize, 5, 10, 15];
        // task-level backward difference structure for the entry's own output nibble ff
        let dz4 = [f.mul(f.imc[0][0], ff), f.mul(f.imc[1][0], ff), f.mul(f.imc[2][0], ff), f.mul(f.imc[3][0], ff)];
        // y3 position feeding z3[4c+i]
        let zpos = |c: usize, i: usize| 4 * ((c + i) % 4) + i;
        {
            if !(is_true_task && stop.load(Ordering::Relaxed)) {
            // per e' (all 15 delta-set differences): dx3_all
            let mut dx3_all = [[0u8; 16]; 15];
            for ep in 1..16u8 {
                let dx2 = [f.mul(f.mc[0][0], ep), f.mul(f.mc[1][0], ep), f.mul(f.mc[2][0], ep), f.mul(f.mc[3][0], ep)];
                let mut dy2 = [0u8; 16];
                for j in 0..4 { dy2[j] = f.sbox[x2c0[j] as usize] ^ f.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
                dx3_all[(ep - 1) as usize] = mc(f, &sr(&dy2));
            }
            let dx3 = dx3_all[(e - 1) as usize];
            for zc in 0..(1u32 << 16) {
                let z4c0 = [(zc & 0xf) as u8, ((zc >> 4) & 0xf) as u8, ((zc >> 8) & 0xf) as u8, ((zc >> 12) & 0xf) as u8];
                let mut x4d = [0u8; 4]; let mut dx4o = [0u8; 16];
                for r in 0..4 {
                    x4d[r] = f.inv_sbox[z4c0[r] as usize];
                    dx4o[dgp[r]] = x4d[r] ^ f.inv_sbox[(z4c0[r] ^ dz4[r]) as usize];
                }
                let dy3o = isr(&imc(f, &dx4o));
                // round-3 branch sets for the entry's own differential (e, ff)
                let mut sets: [&Vec<u8>; 16] = [&ddt[1][0]; 16]; let mut dead = false;
                for j in 0..16 {
                    if dx3[j] == 0 || dy3o[j] == 0 { dead = true; break; }
                    let sset = &ddt[dx3[j] as usize][dy3o[j] as usize];
                    if sset.is_empty() { dead = true; break; }
                    sets[j] = sset;
                }
                if dead { continue; }
                // ---- enumerate branch choices in mixed-radix Gray order with incremental state
                let radix: [usize; 16] = std::array::from_fn(|j| sets[j].len());
                let mut idx = [0usize; 16]; let mut dir = [1i8; 16];
                let mut x3 = [0u8; 16];
                for j in 0..16 { x3[j] = sets[j][0]; }
                // incremental state for all 15 e'
                let mut dy3 = [[0u8; 16]; 15]; let mut dx4d = [[0u8; 4]; 15]; let mut dy4 = [[0u8; 4]; 15]; let mut dcur = [0u8; 15];
                let mut dinv = [0u8; 15]; let mut chi: u16 = 0; let (mut s1, mut s3, mut s5, mut s7) = (0u8, 0u8, 0u8, 0u8);
                for ep in 0..15 {
                    for j in 0..16 { dy3[ep][j] = f.sbox[x3[j] as usize] ^ f.sbox[(x3[j] ^ dx3_all[ep][j]) as usize]; }
                    for r in 0..4 {
                        let mut acc = 0u8;
                        for i in 0..4 { acc ^= f.mul(f.mc[r][i], dy3[ep][zpos(r, i)]); }
                        dx4d[ep][r] = acc;
                        dy4[ep][r] = f.sbox[x4d[r] as usize] ^ f.sbox[(x4d[r] ^ acc) as usize];
                    }
                    let mut dv = 0u8;
                    for r in 0..4 { dv ^= f.mul(f.mc[0][r], dy4[ep][r]); }
                    dcur[ep] = dv;
                }
                for ep in 0..15 { let iv = f.inv[dcur[ep] as usize]; dinv[ep] = iv; chi ^= 1u16 << iv; let t = f.pow357[iv as usize]; s1 ^= iv; s3 ^= t[0]; s5 ^= t[1]; s7 ^= t[2]; }
                loop {
                    // emit entry
                    ltot += 1;
                    if x2c0 == ax2c0 && x3 == ax3 && x4d == ax4d { lts += 1; }
                    let s2 = f.mul(s1, s1); let s4 = f.mul(s2, s2); let s6 = f.mul(s3, s3);
                    let p7 = f.mul(s1, s6) ^ f.mul(s2, s5) ^ f.mul(s3, s4) ^ s7;
                    let (fp, q7e) = if p7 == 0 { lskip += 1; (0xffff, 0) } else {
                        let alpha = f.pow(f.inv[p7 as usize], 13);
                        (permute16(f, chi, alpha, f.mul(alpha, s1)), p7)
                    };
                    if fp != 0xffff {
                    if ltot & 0xf == 0 { lhist[fp as usize] += 1; }
                    for &qi in &qidx[fp as usize] {
                        lfp += 1;
                        let q = &queries[qi as usize];
                        let sh = f.pow(f.div(q.p7, q7e), 14);
                        if sh == 0 { continue; }
                        let s2 = f.mul(sh, sh);
                        let mut ph = [0u8; 16];
                        for i in 0..15 { ph[(f.mul(s2, dinv[i]) ^ sh) as usize] += 1; }
                        let mut ok = true;
                        for v in 1..16usize { if v != sh as usize && ph[v] != q.hist[v] { ok = false; break; } }
                        // harness-side: instrument the TRUE (entry, query) pair with both tests
                        if x2c0 == ax2c0 && x3 == ax3 && x4d == ax4d && q.kc == true_kc && q.uc == true_uc {
                            let mut exact = true;
                            for v in 0..16usize { if ph[v] != q.hist[v] { exact = false; break; } }
                            alg1_exact.store(exact as u64, Ordering::Relaxed);
                            alg1_corrected.store(ok as u64, Ordering::Relaxed);
                        }
                        if ok {
                            lsurv.push((qi as usize, x2c0, x3, x4d, sh));
                            if is_true_task {
                                completions.fetch_add(1, Ordering::Relaxed);
                                if let Some(mk) = complete(qi as usize, x2c0, x3, x4d, sh) {
                                    *recovered_key.lock().unwrap() = Some(mk);
                                    stop.store(true, Ordering::Relaxed);
                                }
                            }
                        }
                    }
                    }
                    // ---- next Gray step: find lowest movable digit
                    let mut jm = 16usize;
                    for j in 0..16 {
                        if radix[j] < 2 { continue; }
                        let ni = idx[j] as i32 + dir[j] as i32;
                        if ni >= 0 && (ni as usize) < radix[j] { jm = j; break; }
                    }
                    if jm == 16 { break; }
                    for j in 0..jm { dir[j] = -dir[j]; }
                    idx[jm] = (idx[jm] as i32 + dir[jm] as i32) as usize;
                    let newx = sets[jm][idx[jm]];
                    x3[jm] = newx;
                    // incremental update for the flipped nibble jm
                    let jr = jm % 4; let ja = jm / 4;
                    let c = (ja + 4 - jr) % 4; // z3 column receiving y3 position jm
                    for ep in 0..15 {
                        let old = dy3[ep][jm];
                        let new = f.sbox[newx as usize] ^ f.sbox[(newx ^ dx3_all[ep][jm]) as usize];
                        if new == old { continue; }
                        dy3[ep][jm] = new;
                        dx4d[ep][c] ^= f.mul(f.mc[c][jr], old ^ new);
                        let ndy4 = f.sbox[x4d[c] as usize] ^ f.sbox[(x4d[c] ^ dx4d[ep][c]) as usize];
                        dcur[ep] ^= f.mul(f.mc[0][c], dy4[ep][c] ^ ndy4);
                        dy4[ep][c] = ndy4;
                        let (o, nn) = (dinv[ep], f.inv[dcur[ep] as usize]);
                        if o != nn {
                            dinv[ep] = nn; chi ^= (1u16 << o) ^ (1u16 << nn);
                            let (to, tn) = (f.pow357[o as usize], f.pow357[nn as usize]);
                            s1 ^= o ^ nn; s3 ^= to[0] ^ tn[0]; s5 ^= to[1] ^ tn[1]; s7 ^= to[2] ^ tn[2];
                        }
                    }
                }
            }
            }
        }
        task_entries[ti].fetch_add(ltot, Ordering::Relaxed);
        task_hits[ti].fetch_add(lfp, Ordering::Relaxed);
        task_surv[ti].fetch_add(lsurv.len() as u64, Ordering::Relaxed);
        total_entries.fetch_add(ltot, Ordering::Relaxed);
        skipped_entries.fetch_add(lskip, Ordering::Relaxed);
        fp_hits.fetch_add(lfp, Ordering::Relaxed);
        true_seen.fetch_add(lts, Ordering::Relaxed);
        if !lsurv.is_empty() { survivors.lock().unwrap().extend(lsurv); }
        let mut h = hist.lock().unwrap(); for i in 0..65536 { h[i] += lhist[i]; }
        let dt = done_tasks.fetch_add(1, Ordering::Relaxed) + 1;
        if dt % 1024 == 0 { eprintln!("  [progress] {} / {} units done, entries {} , stop={}", dt, units.len(), total_entries.load(Ordering::Relaxed), stop.load(Ordering::Relaxed)); }
    });
    let t_all = t_stream.elapsed().as_secs_f64();
    let tot = total_entries.load(Ordering::Relaxed);
    let h = hist.into_inner().unwrap();
    let coll: f64 = h.iter().map(|&c| (c as f64) * (c as f64)).sum::<f64>() / (tot as f64) / (tot as f64);
    out += &format!("  streamed {} entries (log2 = {:.2}) in {:.1}s; 16-bit chi-star collision entropy = {:.2} bits; fp hits = {} (predicted Q*N*2^-H2 = {:.0})\n",
        tot, (tot as f64).log2(), t_all, -coll.log2(), fp_hits.load(Ordering::Relaxed), queries.len() as f64 * (tot as f64) * coll);
    let te0 = task_entries[0].load(Ordering::Relaxed); let th0 = task_hits[0].load(Ordering::Relaxed); let ts0 = task_surv[0].load(Ordering::Relaxed);
    let ncomp = completions.load(Ordering::Relaxed);
    out += &format!("  TRUE outer task: entries={} (log2 {:.2}) hits={} dispatch survivors={} completion calls={}\n", te0, (te0.max(1) as f64).log2(), th0, ts0, ncomp);
    let mut samp: Vec<(u64, u64, u64)> = vec![];
    for i in 1..ntasks { samp.push((task_entries[i].load(Ordering::Relaxed), task_hits[i].load(Ordering::Relaxed), task_surv[i].load(Ordering::Relaxed))); }
    if !samp.is_empty() {
        let me = samp.iter().map(|s| s.0 as f64).sum::<f64>() / samp.len() as f64;
        let mh = samp.iter().map(|s| s.1 as f64).sum::<f64>() / samp.len() as f64;
        let ms = samp.iter().map(|s| s.2 as f64).sum::<f64>() / samp.len() as f64;
        out += &format!("  SAMPLED tasks (n={}): mean entries={:.0} (log2 {:.2}) mean hits={:.0} mean survivors={:.1} -> extrapolated full table (3600 tasks): entries 2^{:.2}, survivors 2^{:.2}\n",
            samp.len(), me, me.max(1.0).log2(), mh, ms, (me.max(1.0) * 3600.0).log2(), (ms.max(1.0) * 3600.0).log2());
    }
    let recovered = recovered_key.into_inner().unwrap();
    let (rec, ok) = match recovered { Some(mk) => (true, mk == secret_key), None => (false, false) };
    out += &format!("  RESULT: key recovered={} matches secret={} ; true entry seen={} ; Alg1 exact test on true hit: {} ; corrected test: {} (1=pass 0=fail 2=not observed)\n",
        rec, ok, true_seen.load(Ordering::Relaxed), alg1_exact.load(Ordering::Relaxed), alg1_corrected.load(Ordering::Relaxed));
    let kr = KeyRes { attempts_until_pair: attempts_total, rej_s0, rej_bodd, rej_frame,
        n_kc: kcands.len(), n_uc: ucands.len(), n_query: queries.len(), chosen_pt,
        true_entries: te0, true_hits: th0, true_surv: ts0, true_completions: ncomp, true_secs: t_all,
        samp, samp_secs: t_all, recovered: rec, recovered_matches: ok,
        alg1_exact: alg1_exact.load(Ordering::Relaxed) as i8, alg1_corrected: alg1_corrected.load(Ordering::Relaxed) as i8, s_true };
    (out, kr)
}

/// multi-key sweep of the honest sampled-coverage experiment
pub fn e13e_sweep(f: &Gf16, nkeys: u64, ktasks: usize, seed0: u64) -> String {
    let mut out = String::new();
    let mut rows: Vec<KeyRes> = vec![];
    out += &format!("E13e SWEEP: {} keys, K={} sampled outer tasks/key (plus the true outer task, full, with completion)\n", nkeys, ktasks);
    for k in 0..nkeys {
        let t0 = std::time::Instant::now();
        let (log, kr) = e13e(f, seed0.wrapping_add(k.wrapping_mul(0x9e37)), ktasks);
        out += &format!("--- key {k} ({:.1}s) ---\n{log}", t0.elapsed().as_secs_f64());
        rows.push(kr);
    }
    let n = rows.len() as f64;
    let nrec = rows.iter().filter(|r| r.recovered_matches).count();
    out += &format!("\n=== AGGREGATE over {} keys ===\n", rows.len());
    out += &format!("  full-key recovery (trial encryption verified): {}/{}\n", nrec, rows.len());
    out += &format!("  pair-conditioning: mean attempts={:.1}; rejections s=0: {} ; odd bad-slot count: {} ; undefined frame: {} (total over sweep)\n",
        rows.iter().map(|r| r.attempts_until_pair as f64).sum::<f64>() / n,
        rows.iter().map(|r| r.rej_s0).sum::<u64>(), rows.iter().map(|r| r.rej_bodd).sum::<u64>(), rows.iter().map(|r| r.rej_frame).sum::<u64>());
    out += &format!("  mean online: |kc|={:.1} |uc|={:.1} queries={:.1} chosen plaintexts={:.1}\n",
        rows.iter().map(|r| r.n_kc as f64).sum::<f64>() / n, rows.iter().map(|r| r.n_uc as f64).sum::<f64>() / n,
        rows.iter().map(|r| r.n_query as f64).sum::<f64>() / n, rows.iter().map(|r| r.chosen_pt as f64).sum::<f64>() / n);
    let mte = rows.iter().map(|r| r.true_entries as f64).sum::<f64>() / n;
    let mts = rows.iter().map(|r| r.true_surv as f64).sum::<f64>() / n;
    let mtc = rows.iter().map(|r| r.true_completions as f64).sum::<f64>() / n;
    out += &format!("  true task: mean entries 2^{:.2}, mean dispatch survivors {:.1}, mean completions {:.1}\n", mte.max(1.0).log2(), mts, mtc);
    let all_s: Vec<&(u64, u64, u64)> = rows.iter().flat_map(|r| r.samp.iter()).collect();
    if !all_s.is_empty() {
        let me = all_s.iter().map(|s| s.0 as f64).sum::<f64>() / all_s.len() as f64;
        let ms = all_s.iter().map(|s| s.2 as f64).sum::<f64>() / all_s.len() as f64;
        let var = all_s.iter().map(|s| (s.0 as f64 - me).powi(2)).sum::<f64>() / (all_s.len() as f64 - 1.0).max(1.0);
        let se = (var / all_s.len() as f64).sqrt();
        out += &format!("  sampled outer tasks (n={}): mean entries {:.0} +- {:.0} (log2 {:.2}); mean dispatch survivors {:.2}/task\n", all_s.len(), me, se, me.max(1.0).log2(), ms);
        out += &format!("  EXTRAPOLATION (3600 outer tasks): full table 2^{:.2} entries; table-wide dispatch survivors 2^{:.2} per right pair\n",
            (me.max(1.0) * 3600.0).log2(), (ms.max(1e-9) * 3600.0).log2());
    }
    let ae: Vec<i8> = rows.iter().map(|r| r.alg1_exact).collect();
    let ac: Vec<i8> = rows.iter().map(|r| r.alg1_corrected).collect();
    out += &format!("  Alg.1 (exact multiset test) on the true hit: pass {} fail {} unobserved {}; corrected test: pass {} fail {} unobserved {}\n",
        ae.iter().filter(|&&x| x == 1).count(), ae.iter().filter(|&&x| x == 0).count(), ae.iter().filter(|&&x| x == 2).count(),
        ac.iter().filter(|&&x| x == 1).count(), ac.iter().filter(|&&x| x == 0).count(), ac.iter().filter(|&&x| x == 2).count());
    out
}

pub fn e13c(f: &Gf16, seed: u64) -> String {
    let mut out = String::new();
    let mut rng = Rng::new(0xC13 ^ seed.wrapping_mul(0x9e37));
    let secret_key: St = rng.state16().map(|b| b & 0xf);
    let rk_secret = expand_key(f, &secret_key);
    let oracle = |p: &St| -> St { enc(f, &rk_secret, p).c };
    // ---- planted right pair (harness uses key knowledge; also require even bad-slot count b so the raw parity applies)
    eprintln!("  [dbg] entering instance loop");
    let (pair, s_true, e_true, f_true, anchor) = loop {
        let (pr, att) = make_right_pair(f, &mut rng, &secret_key);
        eprintln!("  [dbg] got a right pair after {} attempts", att);
        let t0 = enc(f, &rk_secret, &pr.p0); let t1 = enc(f, &rk_secret, &pr.p1);
        let mut a = [0u8; 16]; let mut w = [0u8; 16];
        let x1 = t0.xs[1];
        let posad = [0usize, 13, 10, 7];
        for v in 0..16u8 {
            let mut x = x1; x[0] ^= v;
            let tr = enc(f, &rk_secret, &plaintext_from_x1(f, &rk_secret, &x));
            a[v as usize] = tr.xs[5][0];
            let mut vv = 0u8;
            for j in 0..4 { vv ^= f.mul(f.imc[0][j], f.inv_sbox[(tr.c[posad[j]] ^ rk_secret[7][posad[j]]) as usize]); }
            w[v as usize] = vv;
        }
        let s = a[0]; let b = (1..16).filter(|&v| a[v] == 0 || a[v] == a[0]).count();
        if s == 0 || b % 2 == 1 { continue; }
        let gt: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(w[0] ^ w[v]) as usize] as usize]).collect();
        if chi_star16_frame(f, &gt).is_none() { continue; } // need P_7 != 0 for the true query (frame variant)
        let mut x2c0 = [0u8; 4]; for j in 0..4 { x2c0[j] = t0.xs[2][j]; }
        let dg = [0usize, 5, 10, 15]; let mut x4d = [0u8; 4]; for j in 0..4 { x4d[j] = t0.xs[4][dg[j]]; }
        let e_true = f.sbox[t0.xs[1][0] as usize] ^ f.sbox[t1.xs[1][0] as usize];
        let f_true = t0.xs[5][0] ^ t1.xs[5][0];
        break (pr, s, e_true, f_true, (x2c0, t0.xs[3], x4d));
    };
    eprintln!("  [milestone] right pair + instance ready");
    let (p0, p1) = (pair.p0, pair.p1);
    let (c0, c1) = (oracle(&p0), oracle(&p1));
    out += "E13c SR(7,4,4,4) FULL KEY RECOVERY (single planted right pair, chosen-plaintext oracle)\n";
    out += &format!("  secret key = {:x?}\n  right pair P0={:x?} P1={:x?}\n  (truth for reference: s = {:x}, outer params e = {:x}, f = {:x})\n", secret_key, p0, p1, s_true, e_true, f_true);
    // ---- online candidate keys from the pair (attacker's view only)
    let dg = [0usize, 5, 10, 15]; let anti = [0usize, 13, 10, 7];
    let dp = xor(&p0, &p1); let dc = xor(&c0, &c1);
    let mut kcands: Vec<[u8; 4]> = vec![];
    for k in 0..(1u32 << 16) {
        let kc = [(k & 0xf) as u8, ((k >> 4) & 0xf) as u8, ((k >> 8) & 0xf) as u8, ((k >> 12) & 0xf) as u8];
        // z0 column 0 = SR(y0) col 0 = (y0[0], y0[5], y0[10], y0[15]) = S(P[dg]^kc)
        let mut d = [0u8; 4];
        for r in 0..4 { d[r] = f.sbox[(p0[dg[r]] ^ kc[r]) as usize] ^ f.sbox[(p1[dg[r]] ^ kc[r]) as usize]; }
        let mut dw = [0u8; 4]; for i in 0..4 { for j in 0..4 { dw[i] ^= f.mul(f.mc[i][j], d[j]); } }
        if dw[0] != 0 && dw[1] == 0 && dw[2] == 0 && dw[3] == 0 { kcands.push(kc); }
    }
    let mut ucands: Vec<[u8; 4]> = vec![];
    for k in 0..(1u32 << 16) {
        let uc = [(k & 0xf) as u8, ((k >> 4) & 0xf) as u8, ((k >> 8) & 0xf) as u8, ((k >> 12) & 0xf) as u8];
        let mut d = [0u8; 4];
        for r in 0..4 { d[r] = f.inv_sbox[(c0[anti[r]] ^ uc[r]) as usize] ^ f.inv_sbox[(c1[anti[r]] ^ uc[r]) as usize]; }
        let mut dz = [0u8; 4]; for i in 0..4 { for j in 0..4 { dz[i] ^= f.mul(f.imc[i][j], d[j]); } }
        if dz[0] != 0 && dz[1] == 0 && dz[2] == 0 && dz[3] == 0 { ucands.push(uc); }
    }
    let true_kc = [rk_secret[0][0], rk_secret[0][5], rk_secret[0][10], rk_secret[0][15]];
    let true_uc = [rk_secret[7][0], rk_secret[7][13], rk_secret[7][10], rk_secret[7][7]];
    out += &format!("  online: k_-1[0,5,10,15] candidates from the input differential = {} (contains true: {}); u_6[0,7,10,13] candidates = {} (contains true: {})\n",
        kcands.len(), kcands.contains(&true_kc), ucands.len(), ucands.contains(&true_uc));
    let _ = dp; let _ = dc;
    // ---- for each (kc, uc): build the delta-set (chosen plaintexts), get its ciphertexts, form G and its raw chi-star
    struct Query { kc: [u8; 4], uc: [u8; 4], w: [u8; 16], g: Vec<u8>, hist: [u8; 16], p7: u8 }
    let mut queries: Vec<Query> = vec![];
    let mut chosen_pt = 2u64;
    for kc in kcands.iter() {
        // delta-set plaintexts around P0 for this candidate k_-1 diagonal
        let mut cts = vec![[0u8; 16]; 16];
        for v in 0..16u8 {
            let mut m = p0;
            if v != 0 {
                // want dx1 col0 = (v,0,0,0): dz0 col0 = MC^-1 (v,0,0,0) -> dy0 at diagonal
                for r in 0..4 {
                    let dz = f.mul(f.imc[r][0], v);
                    let x0 = p0[dg[r]] ^ kc[r];
                    let y0 = f.sbox[x0 as usize] ^ dz;
                    m[dg[r]] = f.inv_sbox[y0 as usize] ^ kc[r];
                }
            }
            cts[v as usize] = oracle(&m);
        }
        chosen_pt += 15;
        for uc in ucands.iter() {
            let mut w = [0u8; 16];
            for v in 0..16usize {
                let mut vv = 0u8;
                for j in 0..4 { vv ^= f.mul(f.imc[0][j], f.inv_sbox[(cts[v][anti[j]] ^ uc[j]) as usize]); }
                w[v] = vv;
            }
            let g: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(w[0] ^ w[v]) as usize] as usize]).collect();
            let mut hist = [0u8; 16]; for &x in &g { hist[x as usize] += 1; }
            let p7q = pm16(f, 7, &power_sums16_fast(f, &g, false));
            queries.push(Query { kc: *kc, uc: *uc, w, g, hist, p7: p7q });
        }
    }
    let qfp: Vec<Option<(u16, u8)>> = queries.iter().map(|q| chi_star16_frame(f, &q.g)).collect();
    let mut qidx: Vec<Vec<u16>> = vec![vec![]; 65536];
    let mut qdead = 0usize;
    for (i, fp) in qfp.iter().enumerate() { match fp { Some((x, _)) => qidx[*x as usize].push(i as u16), None => qdead += 1 } }
    out += &format!("  {} online queries (raw parity, {} with undefined frame skipped), chosen plaintexts used = {} (all inside one 2^16 diagonal structure)\n", queries.len(), qdead, chosen_pt);
    eprintln!("  [milestone] {} queries built; starting stream", queries.len());
    // ---- STREAM the full basic offline table
    let ddt: Vec<Vec<Vec<u8>>> = { let mut d = vec![vec![Vec::new(); 16]; 16]; for din in 1..16usize { for x in 0..16u8 { let dout = f.sbox[x as usize] ^ f.sbox[(x as usize) ^ din]; d[din][dout as usize].push(x); } } d };
    let total_entries = AtomicU64::new(0);
    let dead_outer = AtomicU64::new(0);
    let fp_hits = AtomicU64::new(0);
    let true_seen = AtomicU64::new(0);
    let skipped_entries = AtomicU64::new(0);
    let done_tasks = AtomicU64::new(0);
    let survivors: Mutex<Vec<(usize, [u8; 4], St, [u8; 4], u8)>> = Mutex::new(vec![]); // (query idx, x2c0, x3, x4d, s_hat)
    let hist = Mutex::new(vec![0u64; 65536]);
    let (ax2c0, ax3, ax4d) = anchor;
    // ---- inline key completion for a survivor: DS tests at output rows 1..3 (2^20 each) + key-schedule inversion + trial encryption
    let complete = |qi: usize, x2c0: [u8; 4], x3: St, x4d: [u8; 4], _sh: u8| -> Option<St> {
        let q = &queries[qi];
        let mut cts = vec![[0u8; 16]; 16];
        for v in 0..16u8 {
            let mut m = p0;
            if v != 0 { for r in 0..4 { let dz = f.mul(f.imc[r][0], v); let x0 = p0[dg[r]] ^ q.kc[r]; let y0 = f.sbox[x0 as usize] ^ dz; m[dg[r]] = f.inv_sbox[y0 as usize] ^ q.kc[r]; } }
            cts[v as usize] = oracle(&m);
        }
        let mut k6 = [0u8; 16];
        for r in 0..4 { k6[anti[r]] = q.uc[r]; }
        for i in 1..4usize {
            let mut target: Vec<u8> = (1..16u8).map(|ep| {
                let dx2 = [f.mul(f.mc[0][0], ep), f.mul(f.mc[1][0], ep), f.mul(f.mc[2][0], ep), f.mul(f.mc[3][0], ep)];
                let mut dy2 = [0u8; 16]; for j in 0..4 { dy2[j] = f.sbox[x2c0[j] as usize] ^ f.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
                let dx3 = mc(f, &sr(&dy2));
                let mut dy3 = [0u8; 16]; for j in 0..16 { dy3[j] = f.sbox[x3[j] as usize] ^ f.sbox[(x3[j] ^ dx3[j]) as usize]; }
                let dx4 = mc(f, &sr(&dy3));
                let dgp = [0usize, 5, 10, 15]; let mut d = 0u8;
                for r in 0..4 { let dy = f.sbox[x4d[r] as usize] ^ f.sbox[(x4d[r] ^ dx4[dgp[r]]) as usize]; d ^= f.mul(f.mc[i][r], dy); }
                d
            }).collect();
            target.sort_unstable();
            let c = (4 - i) % 4;
            let poss: Vec<usize> = (0..4).map(|r| 4 * ((c + 4 - r) % 4) + r).collect();
            let mut hit: Option<[u8; 5]> = None; let mut nhit = 0;
            for cand in 0..(1u32 << 20) {
                let kk = [(cand & 0xf) as u8, ((cand >> 4) & 0xf) as u8, ((cand >> 8) & 0xf) as u8, ((cand >> 12) & 0xf) as u8, ((cand >> 16) & 0xf) as u8];
                let mut vals = [0u8; 16];
                for v in 0..16usize {
                    let mut z = 0u8;
                    for r in 0..4 { z ^= f.mul(f.imc[i][r], f.inv_sbox[(cts[v][poss[r]] ^ kk[r]) as usize]); }
                    vals[v] = f.inv_sbox[(z ^ kk[4]) as usize];
                }
                let mut m: Vec<u8> = (1..16).map(|v| vals[0] ^ vals[v]).collect();
                m.sort_unstable();
                if m == target { nhit += 1; if hit.is_none() { hit = Some(kk); } }
            }
            match hit { Some(kk) if nhit == 1 => { for r in 0..4 { k6[poss[r]] = kk[r]; } }, _ => return None }
        }
        let master = invert_key_schedule(f, &k6, 7);
        let rkm = expand_key(f, &master);
        if enc(f, &rkm, &p0).c == c0 && enc(f, &rkm, &p1).c == c1 { Some(master) } else { None }
    };
    let recovered_key: Mutex<Option<St>> = Mutex::new(None);
    let stop = std::sync::atomic::AtomicBool::new(false);
    let mut tasks: Vec<(u8, u8, u8)> = (1..16u8).flat_map(|e| (1..16u8).flat_map(move |ff| (0..16u8).map(move |x0| (e, ff, x0)))).collect();
    if std::env::var("E13C_ONLY_TRUE").is_ok() {
        tasks.retain(|t| *t == (e_true, f_true, ax2c0[0]));
        out += "  [ONLY-TRUE MODE: enumeration restricted to the true (e, f, x2c0[0]) outer task -- cost is meaningless, correctness path is real]\n";
    }
    if let Ok(ti) = std::env::var("E13C_ONETASK") {
        let i: usize = ti.parse().unwrap_or(0);
        let keep = tasks[i % tasks.len()];
        tasks.retain(|t| *t == keep);
        out += &format!("  [ONETASK MODE: task {:?} only]\n", keep);
    }
    if std::env::var("E13C_SMOKE").is_ok() {
        // smoke test: only the true (e, f, x2c0[0]) task plus a few random ones
        let keep = (e_true, f_true, ax2c0[0]);
        tasks.retain(|t| *t == keep || (t.2 == 0 && t.1 == 1));
        out += &format!("  [SMOKE MODE: {} of 3600 outer tasks]
", tasks.len());
    }
    tasks.par_iter().for_each(|&(e, ff, x0)| {
        let mut ltot = 0u64; let mut lfp = 0u64; let mut lts = 0u64; let mut lskip = 0u64;
        let mut lhist = vec![0u64; 65536];
        let mut lsurv: Vec<(usize, [u8; 4], St, [u8; 4], u8)> = vec![];
        let dgp = [0usize, 5, 10, 15];
        // task-level backward difference structure for the entry's own output nibble ff
        let dz4 = [f.mul(f.imc[0][0], ff), f.mul(f.imc[1][0], ff), f.mul(f.imc[2][0], ff), f.mul(f.imc[3][0], ff)];
        // y3 position feeding z3[4c+i]
        let zpos = |c: usize, i: usize| 4 * ((c + i) % 4) + i;
        let t_task = std::time::Instant::now();
        for x1n in 0..16u8 { for x2n in 0..16u8 { for x3n in 0..16u8 {
            if stop.load(Ordering::Relaxed) { break; }
            if std::env::var("E13C_PROG").is_ok() && x3n == 0 && x2n % 4 == 0 { eprintln!("  task ({e},{ff},{x0}) x1n={x1n} x2n={x2n}: entries so far {ltot} in {:.0}s", t_task.elapsed().as_secs_f64()); }
            let x2c0 = [x0, x1n, x2n, x3n];
            if std::env::var("E13C_ONLY_TRUE").is_ok() && x2c0 != ax2c0 { continue; }
            // per e' (all 15 delta-set differences): dx3_all
            let mut dx3_all = [[0u8; 16]; 15];
            for ep in 1..16u8 {
                let dx2 = [f.mul(f.mc[0][0], ep), f.mul(f.mc[1][0], ep), f.mul(f.mc[2][0], ep), f.mul(f.mc[3][0], ep)];
                let mut dy2 = [0u8; 16];
                for j in 0..4 { dy2[j] = f.sbox[x2c0[j] as usize] ^ f.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
                dx3_all[(ep - 1) as usize] = mc(f, &sr(&dy2));
            }
            let dx3 = dx3_all[(e - 1) as usize];
            for zc in 0..(1u32 << 16) {
                let z4c0 = [(zc & 0xf) as u8, ((zc >> 4) & 0xf) as u8, ((zc >> 8) & 0xf) as u8, ((zc >> 12) & 0xf) as u8];
                let mut x4d = [0u8; 4]; let mut dx4o = [0u8; 16];
                for r in 0..4 {
                    x4d[r] = f.inv_sbox[z4c0[r] as usize];
                    dx4o[dgp[r]] = x4d[r] ^ f.inv_sbox[(z4c0[r] ^ dz4[r]) as usize];
                }
                let dy3o = isr(&imc(f, &dx4o));
                // round-3 branch sets for the entry's own differential (e, ff)
                let mut sets: [&Vec<u8>; 16] = [&ddt[1][0]; 16]; let mut dead = false;
                for j in 0..16 {
                    if dx3[j] == 0 || dy3o[j] == 0 { dead = true; break; }
                    let sset = &ddt[dx3[j] as usize][dy3o[j] as usize];
                    if sset.is_empty() { dead = true; break; }
                    sets[j] = sset;
                }
                if dead { continue; }
                // ---- enumerate branch choices in mixed-radix Gray order with incremental state
                let radix: [usize; 16] = std::array::from_fn(|j| sets[j].len());
                let mut idx = [0usize; 16]; let mut dir = [1i8; 16];
                let mut x3 = [0u8; 16];
                for j in 0..16 { x3[j] = sets[j][0]; }
                // incremental state for all 15 e'
                let mut dy3 = [[0u8; 16]; 15]; let mut dx4d = [[0u8; 4]; 15]; let mut dy4 = [[0u8; 4]; 15]; let mut dcur = [0u8; 15];
                let mut dinv = [0u8; 15]; let mut chi: u16 = 0; let (mut s1, mut s3, mut s5, mut s7) = (0u8, 0u8, 0u8, 0u8);
                for ep in 0..15 {
                    for j in 0..16 { dy3[ep][j] = f.sbox[x3[j] as usize] ^ f.sbox[(x3[j] ^ dx3_all[ep][j]) as usize]; }
                    for r in 0..4 {
                        let mut acc = 0u8;
                        for i in 0..4 { acc ^= f.mul(f.mc[r][i], dy3[ep][zpos(r, i)]); }
                        dx4d[ep][r] = acc;
                        dy4[ep][r] = f.sbox[x4d[r] as usize] ^ f.sbox[(x4d[r] ^ acc) as usize];
                    }
                    let mut dv = 0u8;
                    for r in 0..4 { dv ^= f.mul(f.mc[0][r], dy4[ep][r]); }
                    dcur[ep] = dv;
                }
                for ep in 0..15 { let iv = f.inv[dcur[ep] as usize]; dinv[ep] = iv; chi ^= 1u16 << iv; let t = f.pow357[iv as usize]; s1 ^= iv; s3 ^= t[0]; s5 ^= t[1]; s7 ^= t[2]; }
                loop {
                    // emit entry
                    ltot += 1;
                    if x2c0 == ax2c0 && x3 == ax3 && x4d == ax4d { lts += 1; }
                    let s2 = f.mul(s1, s1); let s4 = f.mul(s2, s2); let s6 = f.mul(s3, s3);
                    let p7 = f.mul(s1, s6) ^ f.mul(s2, s5) ^ f.mul(s3, s4) ^ s7;
                    let (fp, q7e) = if p7 == 0 { lskip += 1; (0xffff, 0) } else {
                        let alpha = f.pow(f.inv[p7 as usize], 13);
                        (permute16(f, chi, alpha, f.mul(alpha, s1)), p7)
                    };
                    if fp != 0xffff {
                    if ltot & 0xf == 0 { lhist[fp as usize] += 1; }
                    for &qi in &qidx[fp as usize] {
                        lfp += 1;
                        let q = &queries[qi as usize];
                        let sh = f.pow(f.div(q.p7, q7e), 14);
                        if sh == 0 { continue; }
                        let s2 = f.mul(sh, sh);
                        let mut ph = [0u8; 16];
                        for i in 0..15 { ph[(f.mul(s2, dinv[i]) ^ sh) as usize] += 1; }
                        let mut ok = true;
                        for v in 1..16usize { if v != sh as usize && ph[v] != q.hist[v] { ok = false; break; } }
                        if ok {
                            lsurv.push((qi as usize, x2c0, x3, x4d, sh));
                            if std::env::var("E13C_NOCOMPLETE").is_err() { if let Some(mk) = complete(qi as usize, x2c0, x3, x4d, sh) {
                                *recovered_key.lock().unwrap() = Some(mk);
                                stop.store(true, Ordering::Relaxed);
                            } }
                        }
                    }
                    }
                    // ---- next Gray step: find lowest movable digit
                    let mut jm = 16usize;
                    for j in 0..16 {
                        if radix[j] < 2 { continue; }
                        let ni = idx[j] as i32 + dir[j] as i32;
                        if ni >= 0 && (ni as usize) < radix[j] { jm = j; break; }
                    }
                    if jm == 16 { break; }
                    for j in 0..jm { dir[j] = -dir[j]; }
                    idx[jm] = (idx[jm] as i32 + dir[jm] as i32) as usize;
                    let newx = sets[jm][idx[jm]];
                    x3[jm] = newx;
                    // incremental update for the flipped nibble jm
                    let jr = jm % 4; let ja = jm / 4;
                    let c = (ja + 4 - jr) % 4; // z3 column receiving y3 position jm
                    for ep in 0..15 {
                        let old = dy3[ep][jm];
                        let new = f.sbox[newx as usize] ^ f.sbox[(newx ^ dx3_all[ep][jm]) as usize];
                        if new == old { continue; }
                        dy3[ep][jm] = new;
                        dx4d[ep][c] ^= f.mul(f.mc[c][jr], old ^ new);
                        let ndy4 = f.sbox[x4d[c] as usize] ^ f.sbox[(x4d[c] ^ dx4d[ep][c]) as usize];
                        dcur[ep] ^= f.mul(f.mc[0][c], dy4[ep][c] ^ ndy4);
                        dy4[ep][c] = ndy4;
                        let (o, nn) = (dinv[ep], f.inv[dcur[ep] as usize]);
                        if o != nn {
                            dinv[ep] = nn; chi ^= (1u16 << o) ^ (1u16 << nn);
                            let (to, tn) = (f.pow357[o as usize], f.pow357[nn as usize]);
                            s1 ^= o ^ nn; s3 ^= to[0] ^ tn[0]; s5 ^= to[1] ^ tn[1]; s7 ^= to[2] ^ tn[2];
                        }
                    }
                }
            }
        }}}
        let lsurv_n = lsurv.len();
        total_entries.fetch_add(ltot, Ordering::Relaxed);
        skipped_entries.fetch_add(lskip, Ordering::Relaxed);
        fp_hits.fetch_add(lfp, Ordering::Relaxed);
        true_seen.fetch_add(lts, Ordering::Relaxed);
        if !lsurv.is_empty() { survivors.lock().unwrap().extend(lsurv); }
        let mut h = hist.lock().unwrap(); for i in 0..65536 { h[i] += lhist[i]; }
        let dt = done_tasks.fetch_add(1, Ordering::Relaxed) + 1;
        eprintln!("  [progress] {} tasks done, entries {} , stop={}, this task: entries={} hits={} surv={} in {:.1}s", dt, total_entries.load(Ordering::Relaxed), stop.load(Ordering::Relaxed), ltot, lfp, lsurv_n, t_task.elapsed().as_secs_f64());
    });
    let tot = total_entries.load(Ordering::Relaxed);
    let h = hist.into_inner().unwrap();
    let coll: f64 = h.iter().map(|&c| (c as f64) * (c as f64)).sum::<f64>() / (tot as f64) / (tot as f64);
    out += &format!("  offline table: {} entries streamed (log2 = {:.2}) over {} tasks; {} entries with undefined frame (P_7=0) skipped; true entry present in enumeration: {}\n",
        tot, (tot as f64).log2(), done_tasks.load(Ordering::Relaxed), skipped_entries.load(Ordering::Relaxed), true_seen.load(Ordering::Relaxed));
    out += &format!("  16-bit chi-star collision entropy over the table = {:.2} bits; expected fingerprint hits = queries * N * 2^-H2 = {:.0}; observed = {}\n",
        -coll.log2(), queries.len() as f64 * (tot as f64) * coll, fp_hits.load(Ordering::Relaxed));
    let survivors = survivors.into_inner().unwrap();
    out += &format!("  survivors of the corrected multiset-consistency dispatch = {}\n", survivors.len());
    let recovered = recovered_key.into_inner().unwrap();
    for (si, &(qi, _x2, _x3, _x4, sh)) in survivors.iter().enumerate() {
        let q = &queries[qi];
        out += &format!("     survivor {si}: query kc={:x?} uc={:x?}, s_hat={:x} (kc is true: {}, uc is true: {})\n", q.kc, q.uc, sh, q.kc == true_kc, q.uc == true_uc);
    }
    match recovered {
        Some(mk) => out += &format!("  RESULT: master key recovered = {:x?}; matches the secret key: {}\n", mk, mk == secret_key),
        None => out += "  RESULT: no key recovered\n",
    }
    out
}

// (old streaming diagnostics kept for reference)
#[allow(dead_code)]
pub fn e13c_diag(f: &Gf16, wrong_queries: usize) -> String {
    let mut out = String::new();
    // 1. instance: key + right pair with even bad-slot count b
    let mut rng = Rng::new(0xC13);
    let key: St = rng.state16().map(|b| b & 0xf);
    let rk = expand_key(f, &key);
    let (pair, avals, wvals, x1base, e_true, f_true, anchor, tries) = loop {
        let (pr, _att) = make_right_pair(f, &mut rng, &key);
        let t0 = enc(f, &rk, &pr.p0);
        // online delta-set at x1[0]
        let x1 = t0.xs[1];
        let k6 = rk[7]; let pos = [0usize, 13, 10, 7];
        let mut a = [0u8; 16]; let mut w = [0u8; 16];
        for v in 0..16u8 {
            let mut x = x1; x[0] = x1[0] ^ v;
            let p = plaintext_from_x1(f, &rk, &x);
            let tr = enc(f, &rk, &p);
            a[v as usize] = tr.xs[5][0];
            let mut x6c0 = [0u8; 4];
            for j in 0..4 { x6c0[j] = f.inv_sbox[(tr.c[pos[j]] ^ k6[pos[j]]) as usize]; }
            let mut vv = 0u8; for j in 0..4 { vv ^= f.mul(f.imc[0][j], x6c0[j]); }
            w[v as usize] = vv;
        }
        let s = a[0];
        let b = (1..16).filter(|&v| a[v] == 0 || a[v] == a[0]).count();
        if s == 0 || b % 2 == 1 { continue; }
        // true anchor bytes for the reference message
        let mut x2c0 = [0u8; 4]; for j in 0..4 { x2c0[j] = t0.xs[2][j]; }
        let x3 = t0.xs[3];
        let dg = [0usize, 5, 10, 15]; let mut x4d = [0u8; 4]; for j in 0..4 { x4d[j] = t0.xs[4][dg[j]]; }
        let t1 = enc(f, &rk, &pr.p1);
        // outer param e = dy1[0] = S(x1) ^ S(x1')
        let e_true = f.sbox[t0.xs[1][0] as usize] ^ f.sbox[t1.xs[1][0] as usize];
        let f_true = { let dw4 = xor(&t0.xs[5], &t1.xs[5]); dw4[0] };
        break (pr, a, w, x1, e_true, f_true, (x2c0, x3, x4d), 0u32);
    };
    let _ = (pair.p0, tries, x1base);
    let s = avals[0];
    let g_true: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(wvals[0] ^ wvals[v]) as usize] as usize]).collect();
    let q_true = chi_star16(f, &g_true, false);
    // wrong-key queries: peel with wrong u6 nibbles at the same delta-set ciphertexts
    let mut cts = vec![[0u8; 16]; 16];
    for v in 0..16u8 { let mut x = x1base; x[0] ^= v; let p = plaintext_from_x1(f, &rk, &x); cts[v as usize] = enc(f, &rk, &p).c; }
    let pos = [0usize, 13, 10, 7];
    let mut wrongq: Vec<u16> = vec![];
    let mut wrng = Rng::new(0xBAD);
    while wrongq.len() < wrong_queries {
        let guess = [wrng.byte() & 0xf, wrng.byte() & 0xf, wrng.byte() & 0xf, wrng.byte() & 0xf];
        let true_g = [rk[7][0], rk[7][13], rk[7][10], rk[7][7]];
        if guess == true_g { continue; }
        let mut w = [0u8; 16];
        for v in 0..16usize {
            let mut vv = 0u8;
            for j in 0..4 { vv ^= f.mul(f.imc[0][j], f.inv_sbox[(cts[v][pos[j]] ^ guess[j]) as usize]); }
            w[v] = vv;
        }
        let g: Vec<u8> = (1..16).map(|v| f.inv[f.linv[(w[0] ^ w[v]) as usize] as usize]).collect();
        wrongq.push(chi_star16(f, &g, false));
    }
    // 2. stream the FULL basic offline table
    let ddt: Vec<Vec<Vec<u8>>> = { let mut d = vec![vec![Vec::new(); 16]; 16]; for din in 1..16usize { for x in 0..16u8 { let dout = f.sbox[x as usize] ^ f.sbox[(x as usize) ^ din]; d[din][dout as usize].push(x); } } d };
    let total_entries = AtomicU64::new(0);
    let dead_outer = AtomicU64::new(0);
    let hits_true = AtomicU64::new(0);        // fingerprint == true query
    let hits_true_pass = AtomicU64::new(0);   // ... and pass the corrected multiset consistency
    let true_entry_seen = AtomicU64::new(0);  // entry equal to the true anchor state
    let true_entry_pass = AtomicU64::new(0);
    let hits_wrong: Vec<AtomicU64> = (0..wrongq.len()).map(|_| AtomicU64::new(0)).collect();
    let hist = Mutex::new(vec![0u64; 65536]); // fingerprint distribution over the whole table
    let (ax2c0, ax3, ax4d) = anchor;
    // outer loop: e (15) x fpar (15) x x2c0 (2^16) x z4c0 (2^16); parallelize over (e, fpar, x2c0[0..2])
    let tasks: Vec<(u8, u8, u8)> = (1..16u8).flat_map(|e| (1..16u8).flat_map(move |ff| (0..16u8).map(move |x0| (e, ff, x0)))).collect();
    tasks.par_iter().for_each(|&(e, ff, x0)| {
        let mut ltot = 0u64; let mut ldead = 0u64; let mut lht = 0u64; let mut lhtp = 0u64;
        let mut lts = 0u64; let mut ltp = 0u64;
        let mut lhw = vec![0u64; wrongq.len()];
        let mut lhist = vec![0u64; 65536];
        let dx2 = [f.mul(f.mc[0][0], e), f.mul(f.mc[1][0], e), f.mul(f.mc[2][0], e), f.mul(f.mc[3][0], e)];
        // backward from ff: dz4 col0 = imc * (ff,0,0,0); dy4 diag; needs z4c0 values -> dx4 diag
        let dz4 = [f.mul(f.imc[0][0], ff), f.mul(f.imc[1][0], ff), f.mul(f.imc[2][0], ff), f.mul(f.imc[3][0], ff)];
        let s2 = f.mul(s, s);
        for x1n in 0..16u8 { for x2n in 0..16u8 { for x3n in 0..16u8 {
            let x2c0 = [x0, x1n, x2n, x3n];
            let mut dy2 = [0u8; 16];
            for j in 0..4 { dy2[j] = f.sbox[x2c0[j] as usize] ^ f.sbox[(x2c0[j] ^ dx2[j]) as usize]; }
            let dx3 = mc(f, &sr(&dy2));
            for zc in 0..(1u32 << 16) {
                let z4c0 = [(zc & 0xf) as u8, ((zc >> 4) & 0xf) as u8, ((zc >> 8) & 0xf) as u8, ((zc >> 12) & 0xf) as u8]; // = y4 diagonal values
                let mut x4d = [0u8; 4]; let mut dx4 = [0u8; 16]; let dg = [0usize, 5, 10, 15];
                for r in 0..4 {
                    x4d[r] = f.inv_sbox[z4c0[r] as usize];
                    dx4[dg[r]] = f.inv_sbox[z4c0[r] as usize] ^ f.inv_sbox[(z4c0[r] ^ dz4[r]) as usize];
                }
                if (0..4).any(|r| dz4[r] == 0) { ldead += 1; continue; }
                let dy3 = isr(&imc(f, &dx4));
                // round-3 branch sets
                let mut sets: [&Vec<u8>; 16] = [&ddt[1][0]; 16]; let mut dead = false;
                for j in 0..16 {
                    if dx3[j] == 0 || dy3[j] == 0 { dead = true; break; }
                    let sset = &ddt[dx3[j] as usize][dy3[j] as usize];
                    if sset.is_empty() { dead = true; break; }
                    sets[j] = sset;
                }
                if dead { ldead += 1; continue; }
                // enumerate branch choices (product over 16 nibbles of |set|)
                let mut idx = [0usize; 16];
                loop {
                    let mut x3 = [0u8; 16];
                    for j in 0..16 { x3[j] = sets[j][idx[j]]; }
                    // entry -> multiset of 15 differences d(e') for e' = 1..15
                    let ds: Vec<u8> = (1..16u8).map(|ep| keyless_d16(f, &x2c0, &x3, &x4d, ep)).collect();
                    let dinv: Vec<u8> = ds.iter().map(|&d| f.inv[d as usize]).collect();
                    let fp = chi_star16(f, &dinv, false);
                    ltot += 1;
                    if ltot & 0xf == 0 { lhist[fp as usize] += 1; }
                    let is_true = x2c0 == ax2c0 && x3 == ax3 && x4d == ax4d;
                    if is_true { lts += 1; }
                    if fp == q_true {
                        lht += 1;
                        // corrected multiset consistency: recover s from Q_7 ratio then compare mod {0,s}
                        let poff = power_sums16(f, &dinv, false);
                        let pon = power_sums16(f, &g_true, false);
                        let (p7, q7) = (pm16(f, 7, &pon), pm16(f, 7, &poff));
                        let mut pass = false;
                        if p7 != 0 && q7 != 0 {
                            let sh = f.pow(f.div(p7, q7), 14); // (P/Q)^{1/14}; 14^{-1} mod 15 = 14
                            let mut pred: Vec<u8> = dinv.iter().map(|&d| f.mul(f.mul(sh, sh), d) ^ sh).filter(|&x| x != 0 && x != sh).collect();
                            let mut gg: Vec<u8> = g_true.iter().cloned().filter(|&x| x != 0 && x != sh).collect();
                            pred.sort_unstable(); gg.sort_unstable();
                            pass = pred == gg;
                        }
                        let _ = s2;
                        if pass { lhtp += 1; if is_true { ltp += 1; } }
                    }
                    for (qi, &q) in wrongq.iter().enumerate() { if fp == q { lhw[qi] += 1; } }
                    // next combination
                    let mut j = 0;
                    loop {
                        idx[j] += 1;
                        if idx[j] < sets[j].len() { break; }
                        idx[j] = 0; j += 1;
                        if j == 16 { break; }
                    }
                    if j == 16 { break; }
                }
            }
        }}}
        total_entries.fetch_add(ltot, Ordering::Relaxed);
        dead_outer.fetch_add(ldead, Ordering::Relaxed);
        hits_true.fetch_add(lht, Ordering::Relaxed);
        hits_true_pass.fetch_add(lhtp, Ordering::Relaxed);
        true_entry_seen.fetch_add(lts, Ordering::Relaxed);
        true_entry_pass.fetch_add(ltp, Ordering::Relaxed);
        for qi in 0..wrongq.len() { hits_wrong[qi].fetch_add(lhw[qi], Ordering::Relaxed); }
        let mut h = hist.lock().unwrap(); for i in 0..65536 { h[i] += lhist[i]; }
    });
    let tot = total_entries.load(Ordering::Relaxed);
    let h = hist.into_inner().unwrap();
    let coll: f64 = h.iter().map(|&c| (c as f64) * (c as f64)).sum::<f64>() / (tot as f64) / (tot as f64);
    let h2 = -coll.log2();
    let hw: Vec<u64> = hits_wrong.iter().map(|x| x.load(Ordering::Relaxed)).collect();
    let mean_hw = hw.iter().sum::<u64>() as f64 / hw.len().max(1) as f64;
    out += &format!("E13c SR(7,4,4,4) end-to-end with the FULL basic offline table (streaming):\n");
    out += &format!("  planted key = {:x?}; right pair found; s = a_0 = {:x}; true outer params (e, f) = ({:x}, {:x})\n", key, s, e_true, f_true);
    out += &format!("  offline table: {} entries streamed ({} dead outer tuples); log2(entries) = {:.2}\n", tot, dead_outer.load(Ordering::Relaxed), (tot as f64).log2());
    out += &format!("  chi-star fingerprint = full 16-bit canonical vector; collision entropy over the whole table = {:.2} bits (max 16; AGL-orbit information cap ~19 bits)\n", h2);
    out += &format!("  TRUE query: fingerprint hits = {} (predicted N*2^-H2 = {:.0}); of these, corrected multiset-consistency survivors = {}\n", hits_true.load(Ordering::Relaxed), (tot as f64) * coll, hits_true_pass.load(Ordering::Relaxed));
    out += &format!("  the true entry (exact 24-nibble anchor state) appeared {} time(s) in the enumeration and passed the consistency check {} time(s)\n", true_entry_seen.load(Ordering::Relaxed), true_entry_pass.load(Ordering::Relaxed));
    out += &format!("  {} wrong-key queries: mean fingerprint hits = {:.1} (predicted {:.0}); min {} max {}\n", hw.len(), mean_hw, (tot as f64) * coll, hw.iter().min().cloned().unwrap_or(0), hw.iter().max().cloned().unwrap_or(0));
    out
}


pub fn e13d(f: &Gf16) -> String {
    // stage counters for the right-pair construction
    let mut rng = Rng::new(5);
    let key: St = rng.state16().map(|b| b & 0xf);
    let rk = expand_key(f, &key);
    let mut ddt: Vec<Vec<Vec<u8>>> = vec![vec![Vec::new(); 16]; 16];
    for din in 1..16usize { for x in 0..16u8 { let dout = f.sbox[x as usize] ^ f.sbox[(x as usize) ^ din]; ddt[din][dout as usize].push(x); } }
    let (mut n, mut compat, mut cols_ok, mut chk_ok, mut r4_ok) = (0u64, 0u64, 0u64, 0u64, 0u64);
    while n < 3_000_000 {
        n += 1;
        let e = 1 + (rng.byte() % 15) as u8;
        let mut dvec = [0u8; 4]; for r in 0..4 { dvec[r] = 1 + (rng.byte() % 15) as u8; }
        let dx2 = [f.mul(f.mc[0][0], e), f.mul(f.mc[1][0], e), f.mul(f.mc[2][0], e), f.mul(f.mc[3][0], e)];
        let mut dy2 = [0u8; 16];
        for j in 0..4 { let opts: Vec<u8> = (1..16u8).filter(|&g| !ddt[dx2[j] as usize][g as usize].is_empty()).collect(); dy2[j] = opts[(rng.byte() as usize) % opts.len()]; }
        let dx3 = mc(f, &sr(&dy2));
        let mut dx4 = [0u8; 16]; let diag = [0usize, 5, 10, 15];
        for r in 0..4 { dx4[diag[r]] = dvec[r]; }
        let dy3 = isr(&imc(f, &dx4));
        if (0..16).any(|j| dx3[j] == 0 || dy3[j] == 0 || ddt[dx3[j] as usize][dy3[j] as usize].is_empty()) { continue; }
        compat += 1;
        let mut x2 = [0u8; 16]; let mut x3 = [0u8; 16]; let mut fail = false;
        for c in 0..4 {
            let r0 = (4 - c) % 4;
            let sols_col0 = &ddt[dx2[r0] as usize][dy2[r0] as usize];
            let mut found = None;
            'outer: for &x2v in sols_col0.iter() { let y_r0 = f.sbox[x2v as usize];
                for a in 0..16u8 { for b in 0..16u8 { for cc in 0..16u8 {
                    let mut zcol = [0u8; 4]; let mut it = [a, b, cc].into_iter();
                    for row in 0..4 { zcol[row] = if row == r0 { y_r0 } else { it.next().unwrap() }; }
                    let mut good = true; let mut col = [0u8; 4];
                    for i in 0..4 { let mut v = 0u8; for j in 0..4 { v ^= f.mul(f.mc[i][j], zcol[j]); } v ^= rk[3][4 * c + i];
                        if !ddt[dx3[4*c+i] as usize][dy3[4*c+i] as usize].contains(&v) { good = false; break; } col[i] = v; }
                    if good { found = Some((x2v, zcol, col)); break 'outer; }
                }}} }
            match found { None => { fail = true; break; } Some((x2v, zcol, col)) => {
                x2[r0] = x2v;
                for row in 0..4 { let src = 4 * ((c + row) % 4) + row; if row != r0 { x2[src] = f.inv_sbox[zcol[row] as usize]; } }
                for i in 0..4 { x3[4 * c + i] = col[i]; } } }
        }
        if fail { continue; }
        cols_ok += 1;
        let x3_chk = xor(&mc(f, &sr(&sb(f, &x2))), &rk[3]);
        if x3_chk != x3 { continue; }
        chk_ok += 1;
        let mut x3p = x3; for j in 0..16 { x3p[j] ^= dx3[j]; }
        let x4 = xor(&mc(f, &sr(&sb(f, &x3))), &rk[4]);
        let x4p = xor(&mc(f, &sr(&sb(f, &x3p))), &rk[4]);
        if xor(&x4, &x4p) != dx4 { continue; }
        let dw4 = xor(&mc(f, &sr(&sb(f, &x4))), &mc(f, &sr(&sb(f, &x4p))));
        if dw4[0] == 0 || (1..16).any(|j| dw4[j] != 0) { continue; }
        r4_ok += 1;
    }
    format!("E13d pair-construction stages over {n} attempts: round3-compat {compat}, column-search ok {cols_ok}, x3 consistency ok {chk_ok}, round-4 collapse (right pair) {r4_ok}\n")
}

// ---------------------------------------------------------------- E13s: honest sizing
/// Black-box harvesting sizing: for `nkeys` random keys, query one full diagonal
/// structure (2^16 chosen plaintexts), keep every pair whose ciphertext
/// difference vanishes outside positions {0,13,10,7} (attacker-visible
/// filter), then -- harness side only -- label each surviving pair as a true
/// right pair (dx1 col0 = (*,0,0,0) and dx5 single nibble at 0) or not, and
/// count the DDT-enumerated online key candidates it would generate.
pub fn e13s(f: &Gf16, nkeys: u64) -> String {
    use rayon::prelude::*;
    let dg = [0usize, 5, 10, 15]; let anti = [0usize, 13, 10, 7];
    let res: Vec<(u64, u64, u64, u64, f64)> = (0..nkeys).into_par_iter().map(|ki| {
        let mut rng = Rng::new(0x5123 ^ ki.wrapping_mul(0x9e37_79b9_7f4a_7c15));
        let key: St = rng.state16().map(|b| b & 0xf);
        let rk = expand_key(f, &key);
        let base: St = rng.state16().map(|b| b & 0xf);
        let t0 = std::time::Instant::now();
        // structure: all 2^16 diagonal values
        let mut ps: Vec<St> = Vec::with_capacity(1 << 16);
        let mut cs: Vec<St> = Vec::with_capacity(1 << 16);
        for v in 0..(1u32 << 16) {
            let mut p = base;
            for r in 0..4 { p[dg[r]] = ((v >> (4 * r)) & 0xf) as u8; }
            ps.push(p);
            cs.push(enc(f, &rk, &p).c);
        }
        // bucket by the 12 inactive ciphertext nibbles
        let mut buckets: std::collections::HashMap<u64, Vec<u32>> = std::collections::HashMap::new();
        for (i, c) in cs.iter().enumerate() {
            let mut k = 0u64;
            for j in 0..16 { if !anti.contains(&j) { k = (k << 4) | c[j] as u64; } }
            buckets.entry(k).or_default().push(i as u32);
        }
        let (mut npair, mut ntrue, mut nquery) = (0u64, 0u64, 0u64);
        for idxs in buckets.values() {
            if idxs.len() < 2 { continue; }
            for a in 0..idxs.len() { for b in a + 1..idxs.len() {
                let (i0, i1) = (idxs[a] as usize, idxs[b] as usize);
                let (p0, p1, c0, c1) = (ps[i0], ps[i1], cs[i0], cs[i1]);
                // attacker-visible candidate counts
                let mut nkc = 0u64;
                for k in 0..(1u32 << 16) {
                    let kc = [(k & 0xf) as u8, ((k >> 4) & 0xf) as u8, ((k >> 8) & 0xf) as u8, ((k >> 12) & 0xf) as u8];
                    let mut d = [0u8; 4];
                    for r in 0..4 { d[r] = f.sbox[(p0[dg[r]] ^ kc[r]) as usize] ^ f.sbox[(p1[dg[r]] ^ kc[r]) as usize]; }
                    let mut dw = [0u8; 4]; for i in 0..4 { for j in 0..4 { dw[i] ^= f.mul(f.mc[i][j], d[j]); } }
                    if dw[0] != 0 && dw[1] == 0 && dw[2] == 0 && dw[3] == 0 { nkc += 1; }
                }
                let mut nuc = 0u64;
                for k in 0..(1u32 << 16) {
                    let uc = [(k & 0xf) as u8, ((k >> 4) & 0xf) as u8, ((k >> 8) & 0xf) as u8, ((k >> 12) & 0xf) as u8];
                    let mut d = [0u8; 4];
                    for r in 0..4 { d[r] = f.inv_sbox[(c0[anti[r]] ^ uc[r]) as usize] ^ f.inv_sbox[(c1[anti[r]] ^ uc[r]) as usize]; }
                    let mut dz = [0u8; 4]; for i in 0..4 { for j in 0..4 { dz[i] ^= f.mul(f.imc[i][j], d[j]); } }
                    if dz[0] != 0 && dz[1] == 0 && dz[2] == 0 && dz[3] == 0 { nuc += 1; }
                }
                npair += 1; nquery += nkc * nuc;
                // harness-side label (uses the key; not visible to the attack)
                let t0r = enc(f, &rk, &p0); let t1r = enc(f, &rk, &p1);
                let dx1 = xor(&t0r.xs[1], &t1r.xs[1]); let dx5 = xor(&t0r.xs[5], &t1r.xs[5]);
                if dx1[0] != 0 && (1..16).all(|j| dx1[j] == 0) && dx5[0] != 0 && (1..16).all(|j| dx5[j] == 0) { ntrue += 1; }
            }}
        }
        (npair, ntrue, nquery, buckets.len() as u64, t0.elapsed().as_secs_f64())
    }).collect();
    let tot_pairs: u64 = res.iter().map(|r| r.0).sum();
    let tot_true: u64 = res.iter().map(|r| r.1).sum();
    let tot_q: u64 = res.iter().map(|r| r.2).sum();
    let mut out = format!("E13s honest harvesting sizing over {} keys, one 2^16 diagonal structure each (2^16 chosen plaintexts/key)\n", nkeys);
    out += &format!("  candidate pairs (ciphertext diff confined to anti-diagonal): total {} = {:.3}/key; true right pairs among them: {} = {:.3}/key\n",
        tot_pairs, tot_pairs as f64 / nkeys as f64, tot_true, tot_true as f64 / nkeys as f64);
    out += &format!("  online (kc,uc) queries generated by ALL candidate pairs: total {} = {:.1}/key = 2^{:.2}/key\n",
        tot_q, tot_q as f64 / nkeys as f64, (tot_q as f64 / nkeys as f64).log2());
    for (ki, r) in res.iter().enumerate() {
        out += &format!("    key {ki}: cand_pairs={} true={} queries={} ({:.1}s)\n", r.0, r.1, r.2, r.4);
    }
    out
}
