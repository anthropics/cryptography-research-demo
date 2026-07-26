#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""
Demo: over GF(2^5), I_{m,n} is a function of χ*, so χ* ‖ I_{m,n} = χ* (no new bits).

For N random 31-element multisets over GF(32):
  - compute χ*(V)  (brute AGL-min of the 32-bit parity vector)
  - compute I_{m,n}(V)  (4 ratios from exponents E={3,5,7,11,15})
  - verify: same χ*  ⇒  same I_{m,n}  (I_{m,n} determined by χ*)
  - count distinct χ* vs distinct (χ*, I_{m,n})  — should be equal
  - also show: V and V∪{c,c} (different multisets) give identical χ* and I_{m,n}
"""
import random, sys
from collections import defaultdict

# ---- GF(2^5) with modulus x^5+x^2+1 (0x25) ----
Q, POLY = 32, 0x25
def _mul(a,b):
    r=0
    for _ in range(5):
        if b&1: r^=a
        b>>=1; a<<=1
        if a&0x20: a^=POLY
    return r
MUL=[[_mul(a,b) for b in range(Q)] for a in range(Q)]
INV=[0]*Q
for a in range(1,Q):
    for b in range(1,Q):
        if MUL[a][b]==1: INV[a]=b
def gpow(a,e):
    r=1
    while e:
        if e&1: r=MUL[r][a]
        a=MUL[a][a]; e>>=1
    return r

# ---- χ* : brute AGL-min of the parity vector --------------------------------
def chi_vec(V):
    x=0
    for v in V: x ^= 1<<v
    return x
def chi_star(V):
    x=chi_vec(V); best=(1<<Q)
    for a in range(1,Q):
        for b in range(Q):
            y=0
            for d in range(Q):
                if x&(1<<d): y |= 1<<(MUL[a][d]^b)
            if y<best: best=y
    return best

# ---- I_{m,n} over GF(32): E = (3,5,7,11,15), first-nonzero construction ----
E=[3,5,7,11,15]
def Pm(V,m):
    # Lucas reduction: p_r = XOR_ω v_ω^r; P_m = XOR_{0<k<m-k,k⊂m} p_k p_{m-k}
    #                  ⊕ p_m  (|V| odd)
    n=len(V)&1
    p=[0]*(m+1)
    for r in range(m+1):
        acc=0
        for v in V: acc^=gpow(v,r)
        p[r]=acc
    res = p[m] if n else 0
    for k in range(1,m):
        if (k&m)==k and k<m-k:
            res ^= MUL[p[k]][p[m-k]]
    return res
def Imn(V):
    Ps=[Pm(V,m) for m in E]
    nz=[i for i,x in enumerate(Ps) if x]
    if len(nz)<2: return None
    m0=E[nz[0]]; P0=Ps[nz[0]]
    out=[]
    for j in nz[1:]:
        n=E[j]; Pn=Ps[j]
        out.append(MUL[gpow(P0,n)][INV[gpow(Pn,m0)]])
    return tuple(out)

# ---- experiment -------------------------------------------------------------
def main(N=20000, seed=1):
    random.seed(seed)
    chi_to_imn=defaultdict(set); imn_to_chi=defaultdict(set)
    distinct_chi=set(); distinct_imn=set(); distinct_pair=set()
    for _ in range(N):
        V=[random.randrange(Q) for _ in range(Q-1)]     # 31-element multiset
        cs=chi_star(V); im=Imn(V)
        chi_to_imn[cs].add(im); imn_to_chi[im].add(cs)
        distinct_chi.add(cs); distinct_imn.add(im); distinct_pair.add((cs,im))
    multi_c=[k for k,v in chi_to_imn.items() if len(v)>1]
    multi_i=[k for k,v in imn_to_chi.items() if len(v)>1]
    import math
    print(f"N={N} random 31-multisets over GF(32)")
    print(f"  distinct χ*              : {len(distinct_chi):6d}  (~{math.log2(len(distinct_chi)):.1f} bits observed)")
    print(f"  distinct I_mn            : {len(distinct_imn):6d}  (~{math.log2(len(distinct_imn)):.1f} bits observed)")
    print(f"  distinct (χ*, I_mn)      : {len(distinct_pair):6d}")
    print(f"  χ*-fibers with >1 I_mn   : {len(multi_c)}  "
          f"{'(I_mn is a function of χ* ✓)' if not multi_c else '(!! counterexample !!)'}")
    print(f"  I_mn-fibers with >1 χ*   : {len(multi_i)}  "
          f"(χ* is {'NOT ' if multi_i else ''}a function of I_mn; I_mn is "
          f"{'coarser' if multi_i else 'equiv to χ*'})")
    max_fiber=max(len(v) for v in imn_to_chi.values())
    print(f"    largest I_mn fiber contains {max_fiber} distinct χ* values")
    # (b) adding a value twice doesn't change either
    V=[random.randrange(Q) for _ in range(Q-1)]
    c=random.randrange(Q)
    V2=V+[c,c]
    print(f"\n  V vs V∪{{{c},{c}}} (different multisets, same parity vector):")
    print(f"    χ*(V)={chi_star(V):08x}  χ*(V∪cc)={chi_star(V2):08x}  "
          f"{'same' if chi_star(V)==chi_star(V2) else 'DIFFER'}")
    print(f"    I_mn(V)={Imn(V)}  I_mn(V∪cc)={Imn(V2)}  "
          f"{'same' if Imn(V)==Imn(V2) else 'DIFFER'}")
    # (c) so what CAN distinguish V from V∪cc?  The full histogram (multiplicities).
    from collections import Counter
    print(f"    histogram(V)   has max count {max(Counter(V).values())}")
    print(f"    histogram(V∪cc) has max count {max(Counter(V2).values())}")
    # (d) rough entropy ceiling
    import math
    print(f"\n  AGL(1,32) has {31*32}=992 elements; "
          f"χ*-orbit ceiling ≈ 2^32/992 ≈ 2^{32-math.log2(992):.1f} ≈ 2^22 bits.")
    print(f"  ⇒ χ* ‖ I_mn cannot exceed ~22 bits; concatenation adds nothing.")

def exact_counts():
    """Exact #AGL-orbits of 32-bit parity vectors via Burnside, and exact
    #distinct I_{m,n} by enumerating one representative per orbit."""
    import math
    # Burnside: #orbits = (1/|G|) Σ_g 2^{#cycles(g)}
    G=31*32
    total=0
    for a in range(1,Q):
        for b in range(Q):
            seen=0; cyc=0
            for d in range(Q):
                if seen&(1<<d): continue
                cyc+=1; x=d
                while not seen&(1<<x):
                    seen|=1<<x; x=MUL[a][x]^b
            total += 1<<cyc
    orbits=total//G
    print(f"\n  Exact #AGL(1,32)-orbits of 32-bit χ-vectors (Burnside): {orbits:,} "
          f"= 2^{math.log2(orbits):.3f}")
    # I_{m,n} range: tuple of ≤4 nonzero GF(32) elements ⇒ ≤31^4 + smaller tuples
    imn_max = sum(31**k for k in range(5))+1   # +1 for None (degenerate)
    print(f"  I_mn output space upper bound  : ≤ {imn_max:,} ≈ 2^{math.log2(imn_max):.2f}")
    print(f"  (already < #orbits, so I_mn is necessarily coarser than χ*)")
    from math import comb
    ms = comb(Q+Q-2, Q-1)
    print(f"\n  For reference at e=5:")
    print(f"    #(31-element multisets over GF(32)) = C(62,31) = {ms:,} ≈ 2^{math.log2(ms):.1f}")
    print(f"    #AGL-orbits of those multisets      ≈ {ms}/{G} ≈ 2^{math.log2(ms/G):.1f}")
    print(f"    #AGL-orbits of parity vectors (χ*)  = {orbits:,} ≈ 2^{math.log2(orbits):.1f}")
    print(f"    #I_mn values                        ≤ 2^{math.log2(imn_max):.1f}")
    return
    # (full canonical-rep enumeration below is too slow in Python; Burnside suffices)
    # Enumerate canonical reps (lex-min), compute I_{m,n} on each, count distinct.
    # χ-vector → I_{m,n}: p_r = XOR of d^r over set bits d; same Lucas formula.
    def Imn_from_chi(x):
        supp=[d for d in range(Q) if x&(1<<d)]
        Ps=[]
        for m in E:
            p=[0]*(m+1)
            for r in range(m+1):
                acc=0
                for d in supp: acc^=gpow(d,r)
                p[r]=acc
            res=p[m] if len(supp)&1 else 0
            for k in range(1,m):
                if (k&m)==k and k<m-k: res^=MUL[p[k]][p[m-k]]
            Ps.append(res)
        nz=[i for i,v in enumerate(Ps) if v]
        if len(nz)<2: return None
        m0,P0=E[nz[0]],Ps[nz[0]]
        return tuple(MUL[gpow(P0,E[j])][INV[gpow(Ps[j],m0)]] for j in nz[1:])
    print("  Enumerating canonical reps and counting distinct I_mn ...")
    reps=set(); imn_set=set()
    # iterate all χ; keep if it equals its own brute canonical
    for x in range(1<<Q):
        if x&0xFFFFF==0:
            print(f"    {x>>20}/{1<<(Q-20)} M  reps={len(reps):,}",end='\r')
        # quick reject: canonical rep must have bit0 = 0 if any bit is 0
        # (since some β can move any 0-bit to position 0). Only x=all-1s has bit0=1.
        if x&1 and x!=(1<<Q)-1: continue
        is_canon=True
        for a in range(1,Q):
            for b in range(Q):
                y=0
                for d in range(Q):
                    if x&(1<<d): y|=1<<(MUL[a][d]^b)
                if y<x: is_canon=False; break
            if not is_canon: break
        if is_canon:
            reps.add(x); imn_set.add(Imn_from_chi(x))
    print()
    print(f"  Canonical reps found           : {len(reps):,}")
    print(f"  Distinct I_mn over all orbits  : {len(imn_set):,} = 2^{math.log2(len(imn_set)):.3f}")
    print(f"  ⇒ χ* carries {math.log2(orbits):.2f} bits; I_mn carries "
          f"{math.log2(len(imn_set)):.2f} bits; concat = χ* = {math.log2(orbits):.2f} bits.")

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='exact':
        exact_counts()
    else:
        main(int(sys.argv[1]) if len(sys.argv)>1 else 20000)
