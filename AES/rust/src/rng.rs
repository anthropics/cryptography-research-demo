// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
//! Small fast PRNG (xoshiro256**), seeded via splitmix64.
#[derive(Clone)]
pub struct Rng { s: [u64; 4] }

fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

impl Rng {
    pub fn new(seed: u64) -> Rng {
        let mut x = seed;
        Rng { s: [splitmix64(&mut x), splitmix64(&mut x), splitmix64(&mut x), splitmix64(&mut x)] }
    }
    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0]; self.s[3] ^= self.s[1]; self.s[1] ^= self.s[2]; self.s[0] ^= self.s[3];
        self.s[2] ^= t; self.s[3] = self.s[3].rotate_left(45);
        result
    }
    #[inline(always)]
    pub fn next_u32(&mut self) -> u32 { (self.next_u64() >> 32) as u32 }
    #[inline(always)]
    pub fn byte(&mut self) -> u8 { (self.next_u64() >> 56) as u8 }
    pub fn fill(&mut self, buf: &mut [u8]) {
        let mut i = 0;
        while i + 8 <= buf.len() { buf[i..i + 8].copy_from_slice(&self.next_u64().to_le_bytes()); i += 8; }
        while i < buf.len() { buf[i] = self.byte(); i += 1; }
    }
    pub fn state16(&mut self) -> [u8; 16] { let mut s = [0u8; 16]; self.fill(&mut s); s }
    /// uniform in [0, n)
    pub fn below(&mut self, n: u32) -> u32 {
        loop { let x = self.next_u32(); let m = (x as u64) * (n as u64); let l = m as u32;
            if l >= (u32::MAX - n + 1) % n { return (m >> 32) as u32; } }
    }
}
