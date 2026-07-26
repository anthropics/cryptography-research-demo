# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Exact e-sequence from the internal state bytes, using the full index maps
(no simplified diagonal wiring).  Used by test_vectors.py as the golden
reference; prop2.py holds the naive / Gray / U-cache variants."""

from gf256 import SBOX, mul
from index_maps import X5_POS, Y3_TO_X4, MC_ROW0
from aes7r import shift_rows, mix_columns


def e_sequence_exact(x2_0, x3, x4_full, x5):
    """e_w for w=1..255 (e_0=0), from x_2[0], x_3[col0], the full x_4 state,
    and x_5 at the diag-0 positions.  This is the 24-reference-byte version;
    the 10-byte version drops x4_full by computing Δx₄ from x3 + Δy₃ (linear)
    — see prop2.py."""
    e = [0] * 256
    for w in range(1, 256):
        dy2 = SBOX[x2_0 ^ w] ^ SBOX[x2_0]
        dx3 = [mul(c, dy2) for c in (2, 1, 1, 3)]
        dy3 = [SBOX[x3[r] ^ dx3[r]] ^ SBOX[x3[r]] for r in range(4)]
        dx4 = [0] * 16
        for p in range(16):
            for r3 in range(4):
                dx4[p] ^= mul(Y3_TO_X4[p][r3], dy3[r3])
        dy4 = [SBOX[x4_full[p] ^ dx4[p]] ^ SBOX[x4_full[p]] for p in range(16)]
        # Δx₅[X5_POS] = Δw₄[X5_POS] = MC(SR(Δy₄))[X5_POS]
        dw4 = mix_columns(shift_rows(bytes(dy4)))
        dx5 = [dw4[X5_POS[r]] for r in range(4)]
        dy5 = [SBOX[x5[r] ^ dx5[r]] ^ SBOX[x5[r]] for r in range(4)]
        e[w] = 0
        for r in range(4):
            e[w] ^= mul(MC_ROW0[r], dy5[r])
    return e
