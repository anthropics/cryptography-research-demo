# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""7-round AES-128 (the cipher under attack).  Column-major state[16].
Matches FIPS-197 for the first 7 rounds; last round omits MixColumns."""

from gf256 import SBOX, ISBOX, mul

NROUNDS = 7
RCON = (0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1b,0x36)

def key_schedule(K0: bytes) -> list[bytes]:
    """Expand 16-byte master key → NROUNDS+1 round keys."""
    w = list(K0)
    for i in range(4, 4*(NROUNDS+1)):
        t = w[4*(i-1):4*i]
        if i % 4 == 0:
            t = [SBOX[t[1]] ^ RCON[i//4 - 1], SBOX[t[2]], SBOX[t[3]], SBOX[t[0]]]
        w += [w[4*(i-4)+j] ^ t[j] for j in range(4)]
    return [bytes(w[16*r:16*r+16]) for r in range(NROUNDS+1)]

def sub_bytes(s):  return bytes(SBOX[b]  for b in s)
def isub_bytes(s): return bytes(ISBOX[b] for b in s)

# ShiftRows on column-major s[row + 4*col]: row r shifts left by r.
_SR  = [((c + r) % 4)*4 + r for c in range(4) for r in range(4)]
_ISR = [((c - r) % 4)*4 + r for c in range(4) for r in range(4)]
def shift_rows(s):  return bytes(s[_SR[i]]  for i in range(16))
def ishift_rows(s): return bytes(s[_ISR[i]] for i in range(16))

_MC  = ((2,3,1,1),(1,2,3,1),(1,1,2,3),(3,1,1,2))
_IMC = ((14,11,13,9),(9,14,11,13),(13,9,14,11),(11,13,9,14))
def _mix(s, M):
    out = bytearray(16)
    for c in range(4):
        col = s[4*c:4*c+4]
        for r in range(4):
            out[4*c+r] = mul(M[r][0],col[0])^mul(M[r][1],col[1])^ \
                         mul(M[r][2],col[2])^mul(M[r][3],col[3])
    return bytes(out)
def mix_columns(s):  return _mix(s, _MC)
def imix_columns(s): return _mix(s, _IMC)

def ark(s, k): return bytes(a^b for a,b in zip(s,k))

def encrypt(P: bytes, rk: list[bytes]) -> bytes:
    s = ark(P, rk[0])
    for r in range(1, NROUNDS):
        s = ark(mix_columns(shift_rows(sub_bytes(s))), rk[r])
    return ark(shift_rows(sub_bytes(s)), rk[NROUNDS])   # no MC in last round

def trace(P: bytes, rk: list[bytes]) -> dict:
    """Return every intermediate state x_r (input to round-r SB), for attack use."""
    x = {1: ark(P, rk[0])}
    s = x[1]
    for r in range(1, NROUNDS):
        s = ark(mix_columns(shift_rows(sub_bytes(s))), rk[r]); x[r+1] = s
    return x

# ─── self-test vs FIPS-197 Appendix B (first round) ───
_K = bytes.fromhex("2b7e151628aed2a6abf7158809cf4f3c")
_P = bytes.fromhex("3243f6a8885a308d313198a2e0370734")
_rk = key_schedule(_K)
assert ark(_P,_rk[0]).hex() == "193de3bea0f4e22b9ac68d2ae9f84808"
assert sub_bytes(ark(_P,_rk[0])).hex() == "d42711aee0bf98f1b8b45de51e415230"

if __name__ == "__main__":
    print("aes7r self-test OK; C =", encrypt(_P,_rk).hex())
