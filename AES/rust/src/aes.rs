// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! AES-128 with round-by-round state tracing. Byte j sits at row j%4,
//! column j/4 (column-major), matching the paper's indexing.
use crate::gf::Gf;

pub type State = [u8; 16];

pub fn sub_bytes(gf: &Gf, s: &State) -> State { let mut r = [0u8; 16]; for i in 0..16 { r[i] = gf.sbox[s[i] as usize] } r }
pub fn inv_sub_bytes(gf: &Gf, s: &State) -> State { let mut r = [0u8; 16]; for i in 0..16 { r[i] = gf.inv_sbox[s[i] as usize] } r }

/// row r rotated left by r: z[4c+r] = y[4*((c+r)%4)+r]
pub fn shift_rows(s: &State) -> State {
    let mut r = [0u8; 16];
    for c in 0..4 { for row in 0..4 { r[4 * c + row] = s[4 * ((c + row) % 4) + row] } }
    r
}
pub fn inv_shift_rows(s: &State) -> State {
    let mut r = [0u8; 16];
    for c in 0..4 { for row in 0..4 { r[4 * ((c + row) % 4) + row] = s[4 * c + row] } }
    r
}
pub fn mix_columns(gf: &Gf, s: &State) -> State {
    let mut r = [0u8; 16];
    for c in 0..4 {
        let a = &s[4 * c..4 * c + 4];
        r[4 * c] = gf.mul2[a[0] as usize] ^ gf.mul3[a[1] as usize] ^ a[2] ^ a[3];
        r[4 * c + 1] = a[0] ^ gf.mul2[a[1] as usize] ^ gf.mul3[a[2] as usize] ^ a[3];
        r[4 * c + 2] = a[0] ^ a[1] ^ gf.mul2[a[2] as usize] ^ gf.mul3[a[3] as usize];
        r[4 * c + 3] = gf.mul3[a[0] as usize] ^ a[1] ^ a[2] ^ gf.mul2[a[3] as usize];
    }
    r
}
pub fn inv_mix_columns(gf: &Gf, s: &State) -> State {
    let mut r = [0u8; 16];
    for c in 0..4 {
        let a = &s[4 * c..4 * c + 4];
        r[4 * c] = gf.mul14[a[0] as usize] ^ gf.mul11[a[1] as usize] ^ gf.mul13[a[2] as usize] ^ gf.mul9[a[3] as usize];
        r[4 * c + 1] = gf.mul9[a[0] as usize] ^ gf.mul14[a[1] as usize] ^ gf.mul11[a[2] as usize] ^ gf.mul13[a[3] as usize];
        r[4 * c + 2] = gf.mul13[a[0] as usize] ^ gf.mul9[a[1] as usize] ^ gf.mul14[a[2] as usize] ^ gf.mul11[a[3] as usize];
        r[4 * c + 3] = gf.mul11[a[0] as usize] ^ gf.mul13[a[1] as usize] ^ gf.mul9[a[2] as usize] ^ gf.mul14[a[3] as usize];
    }
    r
}
pub fn xor(a: &State, b: &State) -> State { let mut r = [0u8; 16]; for i in 0..16 { r[i] = a[i] ^ b[i] } r }

/// AES-128 key expansion: returns 11 round keys rk[0..=10]; rk[0] = key.
/// In the paper's notation k_{-1} = rk[0], k_i = rk[i+1].
pub fn expand_key(gf: &Gf, key: &State) -> Vec<State> {
    let mut w = [[0u8; 4]; 44];
    for i in 0..4 { for j in 0..4 { w[i][j] = key[4 * i + j]; } }
    let rcon = [0x01u8, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36];
    for i in 4..44 {
        let mut t = w[i - 1];
        if i % 4 == 0 {
            t = [gf.sbox[t[1] as usize] ^ rcon[i / 4 - 1], gf.sbox[t[2] as usize], gf.sbox[t[3] as usize], gf.sbox[t[0] as usize]];
        }
        for j in 0..4 { w[i][j] = w[i - 4][j] ^ t[j]; }
    }
    let mut rks = Vec::with_capacity(11);
    for r in 0..11 {
        let mut k = [0u8; 16];
        for c in 0..4 { for j in 0..4 { k[4 * c + j] = w[4 * r + c][j]; } }
        rks.push(k);
    }
    rks
}

/// Trace of an R-round encryption: xs[i] = x_i (round input, after the
/// preceding key add), for i = 0..R; c = ciphertext. Last round omits MC.
pub struct Trace { pub xs: Vec<State>, pub c: State }

pub fn encrypt_trace(gf: &Gf, rk: &[State], p: &State, rounds: usize) -> Trace {
    let mut xs = Vec::with_capacity(rounds + 1);
    let mut x = xor(p, &rk[0]);
    let mut c = [0u8; 16];
    for i in 0..rounds {
        xs.push(x);
        let y = sub_bytes(gf, &x);
        let z = shift_rows(&y);
        if i + 1 < rounds {
            let w = mix_columns(gf, &z);
            x = xor(&w, &rk[i + 1]);
        } else {
            c = xor(&z, &rk[rounds]);
        }
    }
    Trace { xs, c }
}

pub fn encrypt(gf: &Gf, rk: &[State], p: &State, rounds: usize) -> State { encrypt_trace(gf, rk, p, rounds).c }

/// Given the state x_1 (input to round 1), invert round 0 to a plaintext.
pub fn plaintext_from_x1(gf: &Gf, rk: &[State], x1: &State) -> State {
    let w0 = xor(x1, &rk[1]);
    let z0 = inv_mix_columns(gf, &w0);
    let y0 = inv_shift_rows(&z0);
    let x0 = inv_sub_bytes(gf, &y0);
    xor(&x0, &rk[0])
}

pub fn hex(s: &State) -> String { s.iter().map(|b| format!("{:02x}", b)).collect() }
