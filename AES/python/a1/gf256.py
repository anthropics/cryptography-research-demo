# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""GF(256) arithmetic with the AES reduction polynomial x^8+x^4+x^3+x+1.
This is the only algebraic primitive the attack needs.  Runnable, tested."""

ALOG = [0] * 256
LOG  = [0] * 256

def _init():
    x = 1
    for i in range(255):
        ALOG[i] = x
        LOG[x]  = i
        # x *= 3  in GF(256):  x ^= x*2,  where x*2 = (x<<1) ^ (0x11b if msb set)
        x ^= ((x << 1) ^ (0x11b if x & 0x80 else 0)) & 0xff
    ALOG[255] = 1
_init()

def mul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return ALOG[(LOG[a] + LOG[b]) % 255]

def inv(a: int) -> int:
    return 0 if a == 0 else ALOG[255 - LOG[a]]

def gpow(a: int, m: int) -> int:
    return 0 if a == 0 else ALOG[(LOG[a] * m) % 255]

def _rol8(x, n): return ((x << n) | (x >> (8 - n))) & 0xff

def _affine(x: int) -> int:
    """The fixed GF(2)-affine part of the AES S-box."""
    return x ^ _rol8(x,1) ^ _rol8(x,2) ^ _rol8(x,3) ^ _rol8(x,4) ^ 0x63

SBOX  = bytes(_affine(inv(x)) for x in range(256))
ISBOX = bytes(SBOX.index(x)   for x in range(256))

# Inverse of the affine map (needed for the Möbius bridge: strip M to get inv).
# _affine is a bijection on 256 values → build inverse by table.
_AFF_INV = [0]*256
for _x in range(256): _AFF_INV[_affine(_x)] = _x
def m_inv_byte(b: int) -> int:
    """Apply (affine map)^{-1} to b, i.e. undo the linear+0x63 part of SBOX."""
    return _AFF_INV[b]

# ─── self-test ───
assert mul(0x57, 0x83) == 0xc1
assert inv(0x53) == 0xca and mul(0x53, 0xca) == 0x01
assert SBOX[0x00] == 0x63 and SBOX[0x53] == 0xed and SBOX[0x01] == 0x7c
assert all(ISBOX[SBOX[x]] == x for x in range(256))
assert all(m_inv_byte(_affine(x)) == x for x in range(256))  # M^{-1}∘M = id

if __name__ == "__main__":
    print("gf256 self-test OK")
