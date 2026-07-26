// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! The bridge fingerprints: pairwise power sums P_m, the I_{m,n} 12-byte
//! fingerprint, and the chi-star canonicalization.
use crate::gf::Gf;

pub const E15: [u8; 15] = [7, 11, 13, 19, 23, 29, 31, 37, 43, 47, 53, 59, 61, 91, 127];

/// pow_odd[v][j] = v^(2j+1), j = 0..63  (all odd exponents 1..127)
pub struct PowTab { pub row: Vec<[u8; 64]> }

impl PowTab {
    pub fn new(gf: &Gf) -> PowTab {
        let mut row = vec![[0u8; 64]; 256];
        for v in 0..256 {
            for j in 0..64 { row[v][j] = gf.pow(v as u8, 2 * j as u32 + 1); }
        }
        PowTab { row }
    }
}

/// p[r] for r = 0..127 : p[0] = |V| mod 2, p[r] = XOR of v^r.
/// `extra_zeros` counts appended synthetic zeros (they change only p[0]).
pub fn power_sums(gf: &Gf, pt: &PowTab, values: &[u8], extra_zeros: usize) -> [u8; 128] {
    let mut acc = [0u64; 8];
    for &v in values {
        let r = &pt.row[v as usize];
        for w in 0..8 {
            acc[w] ^= u64::from_le_bytes(r[8 * w..8 * w + 8].try_into().unwrap());
        }
    }
    let mut odd = [0u8; 64];
    for w in 0..8 { odd[8 * w..8 * w + 8].copy_from_slice(&acc[w].to_le_bytes()); }
    let mut p = [0u8; 128];
    p[0] = (((values.len() + extra_zeros) & 1) as u8) & 1;
    for j in 0..64 { p[2 * j + 1] = odd[j]; }
    for r in (2..128).step_by(2) { p[r] = gf.mul(p[r / 2], p[r / 2]); }
    p
}

/// P_m from single-index power sums: XOR over interior submask pairs
/// {k, m-k} of p_k p_{m-k}, plus [|V| odd] p_m.  m must be odd, m <= 127.
pub fn pm_from_p(gf: &Gf, m: u8, p: &[u8; 128]) -> u8 {
    let m32 = m as u32;
    let mut acc = 0u8;
    // iterate submasks k of m
    let mut k: u32 = m32;
    loop {
        // next submask (Gosper-free trick): k = (k-1) & m
        k = (k.wrapping_sub(1)) & m32;
        let comp = m32 ^ k;
        if k != 0 && k < comp {
            acc ^= gf.mul(p[k as usize], p[comp as usize]);
        }
        if k == 0 { break; }
    }
    if p[0] & 1 == 1 { acc ^= p[m as usize]; }
    // the {0,m} pair together with the lone-term parity: an extra p_m only when wt(m)=1
    if m32.count_ones() == 1 { acc ^= p[m as usize]; }
    acc
}

/// Direct O(n^2) evaluation of P_m for cross-checking.
pub fn pm_direct(gf: &Gf, m: u8, values: &[u8]) -> u8 {
    let mut acc = 0u8;
    for i in 0..values.len() {
        for j in i + 1..values.len() {
            acc ^= gf.pow(values[i] ^ values[j], m as u32);
        }
    }
    acc
}

pub struct Imn {
    pub bytes: [u8; 12],
    pub m0: u8,
    pub sel: [u8; 12],
}

/// The paper's 12-byte I_{m,n} fingerprint: first nonzero P_{m0} along E,
/// then the next twelve nonzero exponents n_j; byte j = P_{m0}^{n_j}/P_{n_j}^{m0}.
pub fn imn_fingerprint(gf: &Gf, pms: &[u8; 15]) -> Option<Imn> {
    let mut idx = 0usize;
    while idx < 15 && pms[idx] == 0 { idx += 1; }
    if idx >= 15 { return None; }
    let m0 = E15[idx];
    let pm0 = pms[idx];
    let mut bytes = [0u8; 12];
    let mut sel = [0u8; 12];
    let mut cnt = 0usize;
    let mut j = idx + 1;
    while j < 15 && cnt < 12 {
        if pms[j] != 0 {
            let n = E15[j];
            let num = gf.pow(pm0, n as u32);
            let den = gf.pow(pms[j], m0 as u32);
            bytes[cnt] = gf.div(num, den);
            sel[cnt] = n;
            cnt += 1;
        }
        j += 1;
    }
    if cnt < 12 { return None; }
    Some(Imn { bytes, m0, sel })
}

pub fn all_pms(gf: &Gf, p: &[u8; 128]) -> [u8; 15] {
    let mut out = [0u8; 15];
    for (i, &m) in E15.iter().enumerate() { out[i] = pm_from_p(gf, m, p); }
    out
}

// ---------------------------------------------------------------------
// chi-star canonicalization

pub const CHI_CHAIN: [u8; 6] = [7, 11, 13, 23, 31, 127];

pub struct ChiOut {
    pub vec: [u8; 32],       // canonical 256-bit vector (bit d at byte d/8, bit d%8)
    pub chain_level: usize,  // which P_m in CHI_CHAIN gave alpha*
    pub branch: u8,          // 0=n odd, 1=even normal, 2=even & S1=0 (min over 256 beta)
    pub two_roots: bool,
}

pub fn chi_vector(values: &[u8], extra_zeros: usize) -> [u64; 4] {
    let mut chi = [0u64; 4];
    for &v in values { chi[(v >> 6) as usize] ^= 1u64 << (v & 63); }
    if extra_zeros % 2 == 1 { chi[0] ^= 1; }
    chi
}

#[inline(always)]
fn get_bit(chi: &[u64; 4], d: u8) -> u64 { (chi[(d >> 6) as usize] >> (d & 63)) & 1 }

/// apply pi_{a,b}: out[a d + b] = chi[d]
fn permute(gf: &Gf, chi: &[u64; 4], a: u8, b: u8) -> [u8; 32] {
    let mut out = [0u8; 32];
    // out bit t = chi[a^{-1}(t ^ b)]
    let ainv = gf.inv[a as usize];
    for t in 0..256usize {
        let d = gf.mul(ainv, (t as u8) ^ b);
        if get_bit(chi, d) == 1 { out[t >> 3] |= 1 << (t & 7); }
    }
    out
}

pub fn chi_star(gf: &Gf, pt: &PowTab, values: &[u8], extra_zeros: usize) -> Option<ChiOut> {
    let p = power_sums(gf, pt, values, extra_zeros);
    let n_odd = p[0] & 1 == 1;
    // alpha* from the chain
    let mut alpha_star = 0u8;
    let mut level = usize::MAX;
    for (li, &m) in CHI_CHAIN.iter().enumerate() {
        let pmv = pm_from_p(gf, m, &p);
        if pmv != 0 {
            // alpha* = P_m^{-1/m}
            alpha_star = gf.root(gf.inv[pmv as usize], m as u32);
            level = li;
            break;
        }
    }
    if level == usize::MAX { return None; }
    let s1 = p[1];
    let s3 = p[3];
    let chi = chi_vector(values, extra_zeros);
    let mut branch = 0u8;
    let mut two_roots = false;
    let candidates: Vec<u8> = if n_odd {
        vec![gf.mul(alpha_star, s1)]
    } else if s1 != 0 {
        branch = 1;
        // c = S3 / S1^3, solve u^2 + u = c (or c ^ tau if Tr(c)=1)
        let mut c = gf.div(s3, gf.pow(s1, 3));
        if gf.trace(c) == 1 {
            // tau: fixed element with Tr(tau) = 1
            let mut tau = 1u8;
            for x in 1..=255u8 { if gf.trace(x) == 1 { tau = x; break; } }
            c ^= tau;
        }
        let mut u = 0u8;
        let mut found = false;
        for cand in 0..=255u8 {
            if gf.mul(cand, cand) ^ cand == c { u = cand; found = true; break; }
        }
        if !found { return None; }
        two_roots = true;
        let s1p = gf.mul(alpha_star, s1);
        vec![gf.mul(u, s1p), gf.mul(u ^ 1, s1p)]
    } else {
        branch = 2;
        (0..=255u8).collect()
    };
    let mut best: Option<[u8; 32]> = None;
    for &b in &candidates {
        let v = permute(gf, &chi, alpha_star, b);
        best = Some(match best { None => v, Some(prev) => if v < prev { v } else { prev } });
    }
    Some(ChiOut { vec: best.unwrap(), chain_level: level, branch, two_roots })
}
