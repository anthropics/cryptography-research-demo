#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""project.py -- combine the MEASURED per-entry / per-candidate constants (cycles,
same machine) with the COUNTED exponents of the attack to project the full-scale
7-round AES-128 complexities of the three variants at D = 2^105.

NOTHING here runs 2^80..2^96 work: the per-unit constants are measured on
2^24-entry / 2^16-2^20-candidate real-AES samples; the exponents are counted.

Usage: python3 project.py            (constants are filled in from the logs)
"""
import math

# ---------------------------------------------------------------------------
# MEASURED CONSTANTS (core cycles; medians of 5 pinned reps, Intel Xeon Platinum 8481C).
# Filled in from logs/offline3.log, logs/online3.log, logs/denom.log.
# ---------------------------------------------------------------------------
M = dict(
    # OFFLINE cycles per table entry (full 255-element sequence + reduction)
    off_A = 375.7,       # (A) DDT-Gray walk + brute-free chi* fingerprint pair (vector kernel)  [347.35 ticks]
    off_B = 286.5,       # (B) same walk + multiset hash of the raw sequence                          [264.91 ticks]
    off_C = 11493.3,     # (C) cold per-entry trace + multiset hash                                  [10554.6 ticks]
    # ONLINE cycles per (k_-1, u_6) candidate
    on_A_plain  = 2150.3,    # (A) scalar fused-table peel + L^-1 + vector chi* kernel + 2 probes
    on_A_allvec = 648.6,     # (A) GFNI peel + GFNI L^-1 + vector chi* + 2 probes
    on_B_plain  = 199816.9,  # (B) shared peel + 256 x (iSBOX pass + hash + probe)
    on_B_allvec = 117141.0,  # (B) GFNI peel + per-s GFNI iSBOX pass + gather hash + probes
    on_C_plain  = 371236.0,  # (C) published: peel recomputed inside every s-iteration, scalar
    # DENOMINATORS (cycles per AES-128 encryption, same machine, 10^8 blocks)
    den_ttable = 195.7,      # software T-table encryption  [182.16 ticks/block]
    den_aesni  = 5.63,       # AES-NI, 8-wide pipelined throughput  [5.241 ticks/block]
)

# ---------------------------------------------------------------------------
# COUNTED EXPONENTS (assumptions; from the attack geometry, D = 2^105)
# ---------------------------------------------------------------------------
lg = math.log2
D = 105.0
A_Noff, A_Non = 88.0, 88.0          # (A) table entries (x2 parity words stored), online candidates
B_Noff, B_Non = 80.0, 88.0          # (B),(C): DFJ 10-byte table; 2^88 (k_-1,u_6) candidates
C_Noff, C_Non = 80.0, 88.0          # each candidate runs the 256-value s-loop => 2^96 inner iterations
A_entry_bytes = 26.0                # ~12-byte fp per parity x2 + alignment
B_entry_bytes = 16.0                # 8-byte multiset hash + 8-byte payload
C_entry_bytes = 64.0                # DFJ multiset ~512 bits

# ---------------------------------------------------------------------------
# COUNTED LOOKUP LEDGER (160-lookup units), per offline entry / online candidate
# ---------------------------------------------------------------------------
CAES = 160.0
L = dict(
    # offline per entry
    off_A_impl  = 255*(1.0664+1.0) + 8,    # walk SB (measured counter) + 1 accumulate + C_post(chi*')
    off_A_paper = 263.0,                   # paper's final row with UU-cache (c_off = 1.00) -> P = 2^88.72
    off_B       = 255*(1.0664+1.0),        # same walk + 1 multiset-hash read per element
    off_C       = 255*(24.094+1.0),        # cold: 24 SB/element (measured counter) + 1 hash read
    off_C_paper = 255*(20.0+1.0),          # the paper's "~20 S-boxes/element" cold accounting
    # online per (k_-1,u_6) candidate
    on_A_impl   = 255*(4+1+1) + 8,         # as implemented here: 4-iSBOX peel + L^-1 + accumulate, chi*' post
    on_A_paper  = 388.0,                   # paper's final row (DDT-Gray peel + W-cache + chi*')
    on_B        = 255*4 + 256*(255*(1+1) + 1),   # shared peel + s-loop (1 iSBOX + 1 hash per element) + probe
    on_C        = 256*(255*(4+1+1) + 1),         # per-s full peel (4) + iSBOX + hash per element + probe
    on_DFJ_pub  = 256*1280,                # DFJ'13's own accounting: 5 SB/element per s, hash free
)

def fmt2(x): return "2^%.2f" % lg(x) if x > 0 else "-"

def proj(noff, cent, non, ccand):
    """total cycles = 2^noff * cent + 2^non * ccand (in log2 to avoid overflow)"""
    off = noff + lg(cent)
    on  = non  + lg(ccand)
    tot = lg(2**(off-200) + 2**(on-200)) + 200
    return off, on, tot

def row(name, noff, cent, non, ccand):
    off, on, tot = proj(noff, cent, non, ccand)
    return name, off, on, tot

print("== per-unit constants (measured) ==")
for k in ["off_A","off_B","off_C","on_A_plain","on_A_allvec","on_B_plain","on_B_allvec","on_C_plain","den_ttable","den_aesni"]:
    print("  %-12s %10.1f cyc" % (k, M[k]))

rows = []
# natural implementation level per variant: A = ours best (allvec), B = DFJ-opt best (allvec), C = published plain
rows.append(row("(A) OURS         [off: walk+chi*, on: A-allvec]", A_Noff, M["off_A"], A_Non, M["on_A_allvec"]))
rows.append(row("(B) DFJ-OPTIMIZED[off: walk+mset,  on: B-allvec]", B_Noff, M["off_B"], B_Non, M["on_B_allvec"]))
rows.append(row("(C) DFJ-PUBLISHED[off: cold+mset,  on: C-plain ]", C_Noff, M["off_C"], C_Non, M["on_C_plain"]))
rowsP = []
rowsP.append(row("(A) OURS plain   [on: A-plain ]", A_Noff, M["off_A"], A_Non, M["on_A_plain"]))
rowsP.append(row("(B) DFJ-OPT plain[on: B-plain ]", B_Noff, M["off_B"], B_Non, M["on_B_plain"]))
rowsP.append(row("(C) DFJ-PUB plain[on: C-plain ]", C_Noff, M["off_C"], C_Non, M["on_C_plain"]))

def emit(rows, title):
    print("\n== %s ==" % title)
    print("%-46s %12s %12s %12s %12s %12s %12s" % ("variant", "offline cyc", "online cyc", "TOTAL cyc", "in T-table", "in AES-NI", "(M bytes)"))
    for (name, off, on, tot), eb, noff in zip(rows, [A_entry_bytes, B_entry_bytes, C_entry_bytes], [A_Noff, B_Noff, C_Noff]):
        mb = noff + lg(eb)
        print("%-46s %12s %12s %12s %12s %12s %12s" % (name, "2^%.2f"%off, "2^%.2f"%on, "2^%.2f"%tot,
              "2^%.2f"%(tot-lg(M["den_ttable"])), "2^%.2f"%(tot-lg(M["den_aesni"])), "2^%.1f"%mb))
    tA = rows[0][3]
    print("ratios (cycles): total(B)/total(A) = 2^%.2f = %.1fx ; total(C)/total(A) = 2^%.2f = %.1fx" %
          (rows[1][3]-tA, 2**(rows[1][3]-tA), rows[2][3]-tA, 2**(rows[2][3]-tA)))

emit(rows, "PROJECTION at D=2^105, each variant at its natural implementation level")
emit(rowsP, "PROJECTION, like-for-like scalar/plain-peel online rows")

print("\n== counted-lookup ledger (160-lookup units) ==")
ledger = [
    ("(A) OURS as implemented (DDT-Gray only)", A_Noff, L["off_A_impl"], A_Non, L["on_A_impl"]),
    ("(A) OURS paper-final (UU-cache off / W-cache+chi*' on)", A_Noff, L["off_A_paper"], A_Non, L["on_A_paper"]),
    ("(B) DFJ-OPTIMIZED (shared peel, charged hash)", B_Noff, L["off_B"], B_Non, L["on_B"]),
    ("(C) DFJ-AS-PUBLISHED (per-s peel, charged hash)", C_Noff, L["off_C"], C_Non, L["on_C"]),
    ("(C) DFJ'13 own accounting (1280/s-iter, hash free, ~20SB cold)", C_Noff, L["off_C_paper"], C_Non, L["on_DFJ_pub"]),
]
print("%-58s %9s %9s %10s %10s %10s" % ("variant", "l_off", "l_on", "T_off", "T_on", "TOTAL"))
Lrows = []
for name, noff, loff, non, lon in ledger:
    off = noff + lg(loff/CAES); on = non + lg(lon/CAES)
    tot = lg(2**(off-200) + 2**(on-200)) + 200
    Lrows.append((name, off, on, tot))
    print("%-58s %9.0f %9.0f %10s %10s %10s" % (name, loff, lon, "2^%.2f"%off, "2^%.2f"%on, "2^%.2f"%tot))
print("ledger ratios: total(C-published)/total(A paper-final) = 2^%.2f   |  structural N ratio 2^96/2^88 = 2^8.00" %
      (Lrows[4][3] - Lrows[1][3]))
print("               total(B-optimized)/total(A paper-final)   = 2^%.2f" % (Lrows[2][3] - Lrows[1][3]))
print("               total(C-published)/total(A as-implemented) = 2^%.2f" % (Lrows[4][3] - Lrows[0][3]))
print("\nD = 2^%.0f for all variants; online candidate counts: A 2^%d, B/C 2^%d (k_-1,u_6) x 256 s = 2^96 inner." % (D, int(A_Non), int(B_Non)))
