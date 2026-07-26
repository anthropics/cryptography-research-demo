#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
a1_floor.py — A1 FLOOR PROOF @ D=2^{105}, |K_on|=8: cost model & floor derivation

PROVEN FLOOR: T >= 2^{88} * 254 / 160 = 2^{88.667}  (LIN-STRUCT model, per-entry)

Derivation:
  N = 2^{88}              [offline-tuple count, paper cost section]
  |V| = 255               [MULTISET-FORCED lemma, Jordan 1871]
  fp: S_{255}-sym + Borel-inv + b>=88  [SEPARABILITY + BOREL]
  dim(LS(fp)) in {0,8}    [LS-DICHOTOMY lemma]
  ell >= 255 - dim(LS)/8 >= 254  [LINEAR-STRUCTURE, floor-81.1]
  T = N * ell / C_AES >= 2^{88} * 254 / 160
"""

import math

def lg(x): return math.log2(x)

N        = 2**88
V        = 255
C_AES    = 160

print("="*72)
print("A1 FLOOR PROOF @ D=2^{105}, |K_on|=8")
print("="*72)
print()
print("PROVEN CONSTRAINTS:")
print(f"  N = 2^{lg(N):.0f}                  [offline-tuple count]")
print(f"  |V| = {V}                  [MULTISET-FORCED lemma]")
print(f"  fp: Borel-inv, S_255-sym   [SEPARABILITY+BOREL lemmas]")
print(f"  dim(LS(fp)) in {{0, 8}}      [LS-DICHOTOMY lemma]")
print(f"  C_AES = {C_AES}                [AES cost convention: 160 S-box lookups = 1 enc]")
print()

# ---- THE FLOOR ----
ell_floor_8  = 255 - 8//8   # dim(LS)=8 (PGL-inv branch)
ell_floor_0  = 255 - 0      # dim(LS)=0 (AGL-inv branch)
T_floor_8    = N * ell_floor_8 / C_AES
T_floor_0    = N * ell_floor_0 / C_AES

print("THE FLOOR (LIN-STRUCT model, per-entry):")
print(f"  ell >= 255 - dim(LS)/8")
print(f"  dim(LS)=8 (PGL-inv): ell >= {ell_floor_8}, T >= 2^{lg(T_floor_8):.3f}  [PROVEN]")
print(f"  dim(LS)=0 (AGL-inv): ell >= {ell_floor_0}, T >= 2^{lg(T_floor_0):.3f}  [all known fp]")
print()
print(f"  ==> T >= 2^{lg(T_floor_8):.3f}   (PROVEN FLOOR)")
print()

# ---- ACHIEVED (full stack comparison) ----
print("ACHIEVED at |K_on|=8, D=2^{105} (LIN-STRUCT model accounting):")
print(f"  {'technique stack':<40} {'ell/entry':>10} {'T':>12} {'vs floor':>10}")
print(f"  {'-'*40} {'-'*10} {'-'*12} {'-'*10}")

stacks = [
    ("POWINV_wide + chi-CANON' + U7+UU",      255+8,          "best-known"),
    ("POWINV_wide + chi-CANON (C_post=15)",   255+15,         "baseline"),
    ("HIST+CHUNK + chi-CANON' (LIN-STRUCT)",  255+32+8,       "strict accounting"),
    ("POWINV_wide + LOG-LUCAS (C_post=87)",   255+87,         "online-cut"),
    ("DDT-GRAY (no cache, c_w=2.18)",         int(255*2.18),  "earlier floor"),
]
for name, ell, note in stacks:
    T = N*ell/C_AES
    gap = lg(T) - lg(T_floor_8)
    print(f"  {name:<40} {ell:>10} 2^{lg(T):>8.3f} {'+'+format(gap,'.3f'):>10}")

print()
print(f"  {'A1-FLOOR (PROVEN)':<40} {ell_floor_8:>10} 2^{lg(T_floor_8):>8.3f} {'---':>10}")
print()

# ---- PRIOR "FLOORS" (all superseded) ----
print("PRIOR 'FLOORS' at |K_on|=8, D=2^{105} (all superseded by achieved):")
priors = [
    ("earlier floor (c_w>=2)",          2**88 * 255*2 / 160,        "89.67", "SB-CACHE broke c_SB"),
    ("refined floor (c_w>=2.182)",  2**88 * 255*2.182 / 160,    "89.80", "caching bypasses"),
    ("cache-aware floor (c>=1.5,Cp>=87)", 2**88 * (255*1.5+87) / 160, "89.55", "U7-SHIFT+chi-CANON"),
    ("online-cut (c>=1,Cp>=87)",    2**88 * (255*1+87) / 160,   "89.10", "chi-CANON Cp->8"),
]
for name, T, claimed, why_fell in priors:
    print(f"  {name:<32} 2^{lg(T):.2f} (claimed {claimed})  fell: {why_fell}")
print()
print(f"  {'FLOOR (this derivation)':<32} 2^{lg(T_floor_8):.3f}  -- survives all known tricks")
print()

# ---- MODEL DEPENDENCE (honest) ----
print("MODEL DEPENDENCE:")
print(f"  {'model':<35} {'what counts':<30} {'T floor':>12}")
print(f"  {'-'*35} {'-'*30} {'-'*12}")
models = [
    ("LIN-STRUCT (used here)",      "nonlinear byte->wide map", lg(T_floor_8)),
    ("SB-eq strict",           "256B indexed load",        lg(T_floor_8)),
    ("SB-eq permissive",       "load only; varshift free", None),
    ("wall-clock (measured)",  "cycles/enc-cycles",        None),
]
for name, what, flr in models:
    fstr = f"2^{flr:.3f}" if flr else "no proven floor"
    print(f"  {name:<35} {what:<30} {fstr:>12}")
print()
print("  LIN-STRUCT is the natural model: 1 AES-enc = 160 S-boxes = 160")
print("  nonlinear-byte-ops. c_w>=254/255 PROVEN via dim(LS)<=8 (LS-DICHOTOMY).")
print()

# ---- THE PROOF CHAIN ----
print("PROOF CHAIN (LS-DICHOTOMY):")
chain = [
    "fp valid match-fp ==> Borel-inv on E (= AGL-inv on D)  [SEPARABILITY+BOREL]",
    "LS(fp) is S_255-invariant GF(2)-subspace              [A1]",
    "u in LS, u_i != u_j ==> (a,a,0,..)@{i,j} in LS, c=0   [A2]",
    "==> G_a subset LS (dim 254)                           [A3, S_255-closure]",
    "G_a-inv + Borel-inv ==> <=1 bit  [G_a-BOREL-INCOMPAT] [A4, NEW]",
    "  (scale conj G_a->G_{a/lam} ==> all G_b ==> fp=h(sum e); scale ==> 1 bit)",
    "CONTRADICTION b>=2 ==> LS subset diagonal             [A5]",
    "(c,..c) in LS + scale-inv ==> all-diag in LS          [A6]",
    "  ==> E-trans-inv + Borel-inv ==> PGL-inv",
    "==> dim(LS) in {0, 8}",
    "LIN-STRUCT: ell >= 255 - dim(LS)/8 >= 254",
    "T = N*ell/C_AES >= 2^{88}*254/160 = 2^{88.667}",
]
for step in chain:
    print(f"   {step}")
print()
print("="*72)
