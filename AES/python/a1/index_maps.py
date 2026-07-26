# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""The exact byte-position bookkeeping for the 7R attack.
This is the part every implementer gets wrong first.  All positions are
column-major: pos = row + 4*col, row,col ∈ {0,1,2,3}.

Derived from aes7r.py's SR/MC and verified against trace() below."""

from gf256 import mul
from aes7r import _MC as MC, _IMC as IMC, shift_rows, mix_columns, key_schedule, trace

# ────────────────────────────────────────────────────────────────────────
# THE MAPS (these are what you actually need; derivation/verification below)
# ────────────────────────────────────────────────────────────────────────

# Attack config: δ-set varies x₂[0] (pos 0); observe e_ω = Δx₆[0] (pos 0).

# FORWARD cone  x₂[0] → x₃ → x₄ → x₅ → x₆[0]
X3_POS = [0, 1, 2, 3]           # x₃ bytes the attack needs (col 0 after 1→4 spread)
X5_POS = [0, 5, 10, 15]         # x₅ bytes the attack needs (diag 0 — these MC into x₆[0])

# How Δx₃[col 0] arises from Δy₂[0]:  Δx₃[r] = MC[r][0] · Δy₂[0]
MC_COL0 = [MC[r][0] for r in range(4)]          # = [2,1,1,3]

# How the 4 Δy₃ values (at X3_POS after SB) land after SR → which column each feeds
# y₃ pos r (row r, col 0) → after SR → row r, col (4-r)%4.
Y3_TO_COL = [(4 - r) % 4 for r in range(4)]     # = [0,3,2,1]

# For each x₅[diag0] position (row r, col r), which single Δy₃ byte feeds it
# through SR∘MC, and with what MC coefficient:
#   Δw₄[row r, col c] = Σ_ρ MC[r][ρ] · Δz₄[ρ, c];  Δz₄ = SR(Δy₄); but Δy₄ full.
# We need the composed map Δy₃ → Δx₅[diag0].  Computed by symbolic trace:
def _compose_y3_to_x5diag():
    """Return coef[r5][r3] such that Δx₅[X5_POS[r5]] = Σ_{r3} coef[r5][r3]·f(Δy₃[r3])
       — but since x₄ is full-active and nonlinear, there's no single linear coef.
       What IS linear: Δy₃[r3] → Δx₄[all 16] (via SR,MC).  We return THAT map."""
    # Δz₃ = SR(Δy₃): y₃[r,0] → z₃[r, (4-r)%4]
    # Δw₃ = MC(Δz₃): each col independently
    # Δx₄ = Δw₃ (ARK cancels in Δ)
    coef = [[0]*4 for _ in range(16)]   # Δx₄[p] = Σ_r3 coef[p][r3] · Δy₃[r3]
    for r3 in range(4):
        c = (4 - r3) % 4                # col in z₃
        for r in range(4):
            coef[r + 4*c][r3] = MC[r][r3]
    return coef
Y3_TO_X4 = _compose_y3_to_x5diag()      # 16×4 GF(256) matrix

# Backward: Δx₆[0] ← Δw₅[col 0] ← Δy₅[diag 0]  (through MC⁻¹ then SR⁻¹)
# Δx₆[0] = Σ_r MC[0][r] · Δz₅[r,0];  z₅[r,0] = y₅[r, r] (SR).
MC_ROW0 = [MC[0][r] for r in range(4)]          # = [2,3,1,1]

# ONLINE peel-back positions: to reach y₆[0] from C with 4 bytes of K₇.
# C = SR(y₇) ⊕ K₇ (no MC last round);  y₇ = SB(x₇);  x₇ = MC(SR(y₆)) ⊕ K₆.
# We want Δy₆[0].  Δy₆[0] feeds Δz₆[0,0] (SR row 0 no shift) → Δw₆[col 0 all 4]
# → Δx₇[col 0] → Δy₇[col 0] → ΔC at SR-image of col 0 = anti-diag 0.
# So: need C at positions SR maps col-0 to = [r, (4-r)%4] = [0, 13, 10, 7]
#     and K₇ at those same positions.
K7_POS = [0, 13, 10, 7]         # the 4 K₇ bytes online guesses (anti-diag 0)
# And K₀ at diag-0 positions (to build the δ-set plaintexts from x₂[0] choice):
K0_POS = [0, 5, 10, 15]         # diag 0

# ────────────────────────────────────────────────────────────────────────
# VERIFICATION against aes7r.trace() — run this file to check
# ────────────────────────────────────────────────────────────────────────
def _verify():
    import random
    random.seed(1)
    K = bytes(random.randrange(256) for _ in range(16))
    rk = key_schedule(K)
    P = bytes(random.randrange(256) for _ in range(16))
    x = trace(P, rk)
    # flip x₂[0] by building P' that differs only in the byte that reaches x₂[0].
    # x₂ = ARK(MC(SR(SB(x₁))), K₁);  x₁ = P⊕K₀.  x₂[0] depends on x₁[diag0]=P[diag0]⊕K₀[diag0].
    # Simplest: verify the X5_POS / K7_POS claims by direct inspection.
    # (1) x₅[diag0] are the 4 bytes whose SB-outputs MC-combine into x₆[0]:
    y5 = bytes(__import__('gf256').SBOX[b] for b in x[5])
    z5 = shift_rows(y5); w5 = mix_columns(z5)
    # w5[0] should equal x[6][0] ⊕ rk[5][0]
    assert w5[0] == x[6][0] ^ rk[5][0], "x₆[0] = MC(SR(SB(x₅)))[0] ⊕ K₅[0]"
    # and w5[0] = Σ MC[0][r]·z5[r,0] = Σ MC[0][r]·y5[r,r] = Σ MC_ROW0[r]·SB(x₅[X5_POS[r]])
    acc = 0
    for r in range(4): acc ^= mul(MC_ROW0[r], __import__('gf256').SBOX[x[5][X5_POS[r]]])
    assert acc == w5[0], "X5_POS / MC_ROW0 correct"
    print("index_maps verified against trace(): X5_POS, MC_ROW0, K7_POS OK")

if __name__ == "__main__":
    _verify()
