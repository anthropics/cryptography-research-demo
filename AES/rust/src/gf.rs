// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! GF(2^8) arithmetic (AES polynomial x^8+x^4+x^3+x+1), the AES S-box,
//! its linear part L, and the MixColumns matrices — all as tables.

pub struct Gf {
    pub log: [u8; 256],
    pub exp: [u8; 512],  // exp[i] = 3^i, doubled length to avoid a mod
    pub inv: [u8; 256],  // 0^{-1} := 0 (AES convention)
    pub sbox: [u8; 256],
    pub inv_sbox: [u8; 256],
    pub ltab: [u8; 256],   // L(y): the GF(2)-linear part of the S-box
    pub linv: [u8; 256],   // L^{-1}
    pub mul2: [u8; 256], pub mul3: [u8; 256], pub mul9: [u8; 256],
    pub mul11: [u8; 256], pub mul13: [u8; 256], pub mul14: [u8; 256],
}

fn xtime(a: u8) -> u8 { let hi = a & 0x80; let mut r = a << 1; if hi != 0 { r ^= 0x1b } r }

pub fn mul_slow(mut a: u8, mut b: u8) -> u8 {
    let mut p = 0u8;
    while b != 0 {
        if b & 1 != 0 { p ^= a }
        a = xtime(a);
        b >>= 1;
    }
    p
}

impl Gf {
    pub fn new() -> Gf {
        let mut log = [0u8; 256];
        let mut exp = [0u8; 512];
        let mut x = 1u8;
        for i in 0..255 {
            exp[i] = x;
            log[x as usize] = i as u8;
            x = mul_slow(x, 3);
        }
        for i in 255..512 { exp[i] = exp[i - 255]; }
        let mut inv = [0u8; 256];
        for a in 1..256 { inv[a] = exp[255 - log[a] as usize]; }
        // affine part: L(y) then xor 0x63
        let mut ltab = [0u8; 256];
        for y in 0..256usize {
            let mut r = 0u8;
            for i in 0..8 {
                let bit = ((y >> i) & 1) ^ ((y >> ((i + 4) % 8)) & 1)
                    ^ ((y >> ((i + 5) % 8)) & 1) ^ ((y >> ((i + 6) % 8)) & 1)
                    ^ ((y >> ((i + 7) % 8)) & 1);
                r |= (bit as u8) << i;
            }
            ltab[y] = r;
        }
        let mut linv = [0u8; 256];
        for y in 0..256 { linv[ltab[y] as usize] = y as u8; }
        let mut sbox = [0u8; 256];
        let mut inv_sbox = [0u8; 256];
        for a in 0..256 { sbox[a] = ltab[inv[a] as usize] ^ 0x63; }
        for a in 0..256 { inv_sbox[sbox[a] as usize] = a as u8; }
        assert_eq!(sbox[0], 0x63); assert_eq!(sbox[1], 0x7c); assert_eq!(sbox[0x53], 0xed);
        let mut mul2 = [0u8; 256]; let mut mul3 = [0u8; 256]; let mut mul9 = [0u8; 256];
        let mut mul11 = [0u8; 256]; let mut mul13 = [0u8; 256]; let mut mul14 = [0u8; 256];
        for a in 0..256usize {
            let a8 = a as u8;
            mul2[a] = mul_slow(a8, 2); mul3[a] = mul_slow(a8, 3); mul9[a] = mul_slow(a8, 9);
            mul11[a] = mul_slow(a8, 11); mul13[a] = mul_slow(a8, 13); mul14[a] = mul_slow(a8, 14);
        }
        Gf { log, exp, inv, sbox, inv_sbox, ltab, linv, mul2, mul3, mul9, mul11, mul13, mul14 }
    }
    #[inline(always)]
    pub fn mul(&self, a: u8, b: u8) -> u8 {
        if a == 0 || b == 0 { 0 } else {
            self.exp[self.log[a as usize] as usize + self.log[b as usize] as usize]
        }
    }
    #[inline(always)]
    pub fn div(&self, a: u8, b: u8) -> u8 {
        // b != 0 required
        if a == 0 { 0 } else {
            let e = (self.log[a as usize] as i32 - self.log[b as usize] as i32).rem_euclid(255);
            self.exp[e as usize]
        }
    }
    /// a^e with integer exponent e >= 0, and 0^e = 0 for e>0, x^0 = 1.
    #[inline(always)]
    pub fn pow(&self, a: u8, e: u32) -> u8 {
        if e == 0 { return 1 }
        if a == 0 { return 0 }
        let l = self.log[a as usize] as u32;
        self.exp[((l * (e % 255)) % 255) as usize]
    }
    /// a^{1/e} := a^{e^{-1} mod 255}, requires gcd(e,255)=1
    pub fn root(&self, a: u8, e: u32) -> u8 {
        if a == 0 { return 0 }
        let einv = modinv(e as i64, 255) as u32;
        self.pow(a, einv)
    }
    #[inline(always)]
    pub fn trace(&self, x: u8) -> u8 {
        // Tr(x) = x + x^2 + x^4 + ... + x^128, result in {0,1}
        let mut t = 0u8; let mut y = x;
        for _ in 0..8 { t ^= y; y = self.mul(y, y); }
        t
    }
}

pub fn modinv(a: i64, m: i64) -> i64 {
    // extended Euclid
    let (mut old_r, mut r) = (a.rem_euclid(m), m);
    let (mut old_s, mut s) = (1i64, 0i64);
    while r != 0 {
        let q = old_r / r;
        let t = old_r - q * r; old_r = r; r = t;
        let t = old_s - q * s; old_s = s; s = t;
    }
    assert_eq!(old_r, 1, "not invertible");
    old_s.rem_euclid(m)
}
