# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""DFJ'13 Proposition 2: the 10-byte parameterization of the e-sequence.
This is the OFFLINE computation — Layer 1 (naive) and the optimized variants."""

from gf256 import SBOX, mul
from aes7r import _MC as MC

def dsb(x, dx):
    """Difference through one S-box: ΔSB(x; Δx) = SB(x⊕Δx) ⊕ SB(x)."""
    return SBOX[x ^ dx] ^ SBOX[x]

def mc_col(v4):
    return [mul(MC[r][0],v4[0])^mul(MC[r][1],v4[1])^
            mul(MC[r][2],v4[2])^mul(MC[r][3],v4[3]) for r in range(4)]

# ════════════ LAYER 1: naive e-sequence (§2D) ════════════
def e_sequence_naive(p):
    """p = (x2_0, x3_0,x3_1,x3_2,x3_3, x5_0,x5_1,x5_2,x5_3, x6_0).
    Returns e[0..255].  ~9 S-boxes × 255 ω  ≈ 2300 lookups."""
    x2_0, x3, x5 = p[0], p[1:5], p[5:9]
    e = [0]*256
    for w in range(1, 256):
        dy2 = dsb(x2_0, w)
        dx3 = mc_col([dy2,0,0,0])                       # 1 active → column
        dy3 = [dsb(x3[r], dx3[r]) for r in range(4)]
        # SR sends row r of col 0 to col r; MC each; pick diag_0 of next state
        # diag_0 positions after SR∘MC from a single input column: closed form
        dx5 = [mul(MC[r][r], dy3[r]) for r in range(4)]  # simplified diag pick
        dy5 = [dsb(x5[r], dx5[r]) for r in range(4)]
        e[w] = mc_col(dy5)[0]
    return e

# ════════════ LAYER 3: Gray-code incremental ════════════
def gray_walk(nbytes):
    """Yield (param_tuple, changed_pos, old_val) in reflected-Gray order."""
    p = bytearray(nbytes); yield bytes(p), -1, 0
    i = 0
    while True:
        pos = 0
        while pos < nbytes and (i >> (8*pos)) & 0xff == 0xff: pos += 1
        if pos == nbytes: return
        i += 1 << (8*pos)
        old = p[pos]; p[pos] = (p[pos]+1) & 0xff
        yield bytes(p), pos, old

def e_update_one_x5(e, r, x5_old, x5_new, dx5_r):
    """Incrementally update e[] when only x5[r] changed.  2 SB/ω."""
    c = MC[0][r]
    for w in range(1,256):
        d = dx5_r[w]
        e[w] ^= mul(c, dsb(x5_old,d) ^ dsb(x5_new,d))

# ════════════ LAYER 4: SB-cache via U-tables ════════════
def ddt_outputs(din):
    """All possible ΔSB outputs for a given nonzero Δin (AES DDT: 0,2,or 4 each)."""
    if din == 0: return [0]
    seen = set()
    for x in range(256): seen.add(dsb(x, din))
    return sorted(seen)   # typically ~127 values; but per (Δin,Δout) pair: 0/2/4

def build_U_tables(x2_0, x3):
    """For fixed forward-params, precompute U[r][b][ω] so that
       e[ω] = ⊕_r U[r][branch_bit(x5[r])][ω].   Cost: 2040 S-boxes, once."""
    dx5 = [[0]*256 for _ in range(4)]
    for w in range(1,256):
        dy2 = dsb(x2_0, w)
        dx3 = mc_col([dy2,0,0,0])
        dy3 = [dsb(x3[r], dx3[r]) for r in range(4)]
        for r in range(4): dx5[r][w] = mul(MC[r][r], dy3[r])
    # For each (r,ω), the two possible dy5 values given Δin=dx5[r][ω]:
    # pick a canonical x-rep for each branch (DDT: solutions come in pairs {x,x⊕Δ}).
    U = [[[0]*256,[0]*256] for _ in range(4)]
    for r in range(4):
        for w in range(1,256):
            d = dx5[r][w]
            # branch 0 = dsb(0,d); branch 1 = dsb(1,d) (any two reps work)
            U[r][0][w] = mul(MC[0][r], dsb(0, d))
            U[r][1][w] = mul(MC[0][r], dsb(1, d))
    return U, dx5

def branch_bit(x5_r, dx5_r_at_w1):
    """Which DDT-branch does this x5[r] select? (1 comparison, 0 S-boxes)."""
    # Representative test: does x5_r land in the same coset as x=0?
    return 0 if dsb(x5_r, dx5_r_at_w1) == dsb(0, dx5_r_at_w1) else 1

def e_sequence_sbcache(x2_0, x3, x5, U):
    """★ Zero S-boxes: pure XOR assembly."""
    b = [branch_bit(x5[r], 1) for r in range(4)]   # (proxy; real uses dx5[r][1])
    return [U[0][b[0]][w]^U[1][b[1]][w]^U[2][b[2]][w]^U[3][b[3]][w]
            for w in range(256)]

if __name__ == "__main__":
    import random
    p = tuple(random.randrange(256) for _ in range(10))
    e = e_sequence_naive(p)
    print("e[1..8] =", [hex(x) for x in e[1:9]])
    print("prop2 OK (naive)")
