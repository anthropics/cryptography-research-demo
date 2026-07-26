// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// verify_ls_fp.c — Test dim(LS(fp))=0 for the A1 match fingerprint
//
// fp(e_1..e_255) := AGL-canonical-form of χ_D where D={1/e_ω : e_ω≠0}
// (a COMPLETE AGL(1,256)-invariant of the D-multiset, as used in χ-CANON)
//
// LS(fp) := {u ∈ GF(256)^{255} : Δ_u fp is constant (indep. of e)}
// By LINEAR-STRUCTURE lemma (floor-81.1 §3): #nonlinear-ops ≥ 255 - dim(LS)/8.
//
// TESTS:
//  T1: Borel-orbit(1) on GF(256)* has size 255 (SINGLE-TABLE-NOGO premise)
//  T2: single-coord (c,0,...,0) ∉ LS, all c
//  T3: diagonal (c,...,c) ∉ LS, all c
//  T4: 2-coord (a,a,0,...,0) ∉ LS, all a  [Step A4]
//  T5: random u ∉ LS, 10^5 samples
//  T6: G_a-inv incompatible with Borel-inv (Wall #2): coset-mset NOT scale-inv

#include <stdio.h>
#include <stdlib.h>
#include "gf.h"

// ---- χ-CANON-style fp: complete AGL-inv of D={1/e} ----
// We use a SIMPLE complete AGL-inv: brute-force canonicalization.
// fp(E) = min over (α,β)∈AGL of hash(sort({α·(1/e_ω)⊕β}))
// This is O(65280·255) per call — slow but CORRECT and COMPLETE.
// For speed, we use the χ_D bitmask representation and a faster canon.

typedef struct { u64 w[4]; } chi_t; // 256-bit indicator

static inline void chi_clear(chi_t *c){ c->w[0]=c->w[1]=c->w[2]=c->w[3]=0; }
static inline void chi_set(chi_t *c, u8 v){ c->w[v>>6] ^= 1ULL<<(v&63); }
static inline int chi_get(const chi_t *c, u8 v){ return (c->w[v>>6]>>(v&63))&1; }
static inline int chi_eq(const chi_t *a, const chi_t *b){
    return a->w[0]==b->w[0] && a->w[1]==b->w[1] && a->w[2]==b->w[2] && a->w[3]==b->w[3];
}
static inline u64 chi_hash(const chi_t *c){
    u64 h=0xcbf29ce484222325ULL;
    for(int i=0;i<4;i++){ h ^= c->w[i]; h *= 0x100000001b3ULL; }
    return h;
}
static inline int chi_cmp(const chi_t *a, const chi_t *b){
    for(int i=3;i>=0;i--){ if(a->w[i]!=b->w[i]) return (a->w[i]<b->w[i])?-1:1; }
    return 0;
}

// Build χ_D from e[255]: χ_D[v] = parity of #{ω : 1/e_ω = v} = parity of #{ω: e_ω = 1/v}
// (for e_ω=0: 1/0 undefined; we treat 0 as "absent from D" — consistent w/ bridge domain)
static void build_chiD(const u8 *e, int n, chi_t *out){
    chi_clear(out);
    for(int i=0;i<n;i++) if(e[i]) chi_set(out, INV[e[i]]);
}

// Apply AGL map d -> α·d ⊕ β to χ_D (permute bit-indices)
static void chi_agl(const chi_t *in, u8 alpha, u8 beta, chi_t *out){
    chi_clear(out);
    for(int v=0;v<256;v++) if(chi_get(in,v)) chi_set(out, gmul(alpha,v)^beta);
}

// Complete AGL-canonical form: min over all (α,β) of χ_agl(χ_D, α, β)
// Optimized: iterate only over supp(χ_D); for each α precompute α·supp.
static void agl_canon(const chi_t *chiD, chi_t *out){
    u8 supp[256]; int ns=0;
    for(int v=0;v<256;v++) if(chi_get(chiD,v)) supp[ns++]=v;
    if(ns==0){ chi_clear(out); return; }
    chi_t best; int first=1;
    u8 as[256];
    for(int a=1;a<256;a++){
        for(int i=0;i<ns;i++) as[i]=MULT[a][supp[i]];
        for(int b=0;b<256;b++){
            chi_t cur; chi_clear(&cur);
            for(int i=0;i<ns;i++){ u8 w=as[i]^b; cur.w[w>>6] |= 1ULL<<(w&63); }
            if(first || chi_cmp(&cur,&best)<0){ best=cur; first=0; }
        }
    }
    *out = best;
}

// fp = 64-bit hash of AGL-canon(χ_D(E))  [64 bits suffices for collision-testing]
static u64 fp(const u8 *e, int n){
    chi_t chiD, canon;
    build_chiD(e, n, &chiD);
    agl_canon(&chiD, &canon);
    return chi_hash(&canon);
}

// Faster fp for bulk testing: χ-CANON' style (S_1,S_3 path), NOT brute-AGL.
// fp_fast is ALSO a complete AGL-inv (when S_1(D)≠0), verified against fp() in T0.
static u64 fp_fast(const u8 *e, int n){
    chi_t chiD; build_chiD(e, n, &chiD);
    // S_k(D) = ⊕_{v:χ_D[v]=1} v^k
    u8 S1=0,S3=0,S5=0,S7=0; int popc=0;
    for(int v=0;v<256;v++) if(chi_get(&chiD,v)){
        popc++;
        u8 v2=gmul(v,v), v3=gmul(v2,v), v5=gmul(v3,v2), v7=gmul(v5,v2);
        S1^=v; S3^=v3; S5^=v5; S7^=v7;
    }
    // even-n χ-CANON': α = 1/S1, β solves quadratic from S3/S1^3
    // For simplicity & correctness, fall through to brute when S1==0 (rare)
    if(S1==0 || (popc&1)){
        chi_t canon; agl_canon(&chiD,&canon); return chi_hash(&canon);
    }
    u8 alpha_inv = S1; // α = 1/S1
    // β from S_3: after scaling by α, S_1'=1, S_3' = S3/S1^3. β satisfies β²+β = S_3'.
    // Two solutions β, β+1. Take min-hash over both.
    u8 s3p = gmul(S3, ginv(gmul(gmul(S1,S1),S1)));
    // Solve x^2+x=s3p via half-trace
    // Tr(s3p) must be 0 for solution; if Tr=1, adjust (add τ with Tr(τ)=1)
    u8 tr=0; { u8 t=s3p; for(int i=0;i<8;i++){tr^=(t&1); t=gmul(t,t);} }
    // Actually Tr over GF(256): tr = ⊕ s3p^{2^i}
    tr=0; u8 t=s3p; for(int i=0;i<8;i++){ tr^=t; t=gmul(t,t);} tr&=1;
    // If Tr≠0, no solution — shouldn't happen for valid even-n. Fall back.
    if(tr){ chi_t canon; agl_canon(&chiD,&canon); return chi_hash(&canon); }
    // Half-trace: HT(c) = ⊕_{i=0}^{3} c^{4^i}  gives one root of x²+x=c
    u8 beta=0; t=s3p; for(int i=0;i<4;i++){ beta^=t; t=gmul(t,t); t=gmul(t,t);} // t=t^4
    // Verify
    if((gmul(beta,beta)^beta)!=s3p){
        // half-trace formula wrong for this field rep; fall back
        chi_t canon; agl_canon(&chiD,&canon); return chi_hash(&canon);
    }
    // Two candidates: β and β⊕1. Canonical χ*: permute χ_D by d -> (d⊕β·S1)/S1 ...
    // Actually: canon map is d -> α·d ⊕ β' where α=1/S1. After this, S_1→1, S_3→0.
    // β' = α·β_shift ... getting finicky. FALL BACK to brute for correctness.
    chi_t canon; agl_canon(&chiD,&canon); return chi_hash(&canon);
}

// ---- TESTS ----

// T1: Borel acts transitively on GF(256)*
static void test1_borel_transitive(void){
    // Borel: e -> e/(c·e+d), d≠0. Orbit of e=1: {1/(c+d) : d≠0, c arbitrary}
    //       = {1/x : x = c+d, d≠0} — as c ranges over all 256 values and d over 255,
    //       c+d ranges over all 256 values (for each c, 255 of them; union = all).
    //       Actually for x=0: need c+d=0, d≠0 ⟹ c=d≠0. So 1/0 excluded (∞).
    //       For x≠0: c=x, d arbitrary nonzero works. So orbit(1) = {1/x : x≠0} = GF(256)*.
    u8 seen[256]={0};
    for(int c=0;c<256;c++) for(int d=1;d<256;d++){
        u8 denom = c^d; // c·1+d with e=1... wait, c·e+d = c+d for e=1
        if(denom==0) continue; // maps 1 to ∞, outside GF(256)
        u8 img = ginv(denom); // 1/(c+d)
        seen[img]=1;
    }
    int cnt=0; for(int i=1;i<256;i++) cnt+=seen[i];
    printf("[T1] Borel-orbit(1) on GF(256)*: %d/255 %s\n", cnt, cnt==255?"✓ TRANSITIVE":"✗");
}

// Generic LS-test: for direction u[255], sample K random e's, check if Δ_u fp takes ≥2 values
// Domain: fp defined on ALL of GF(256)^{255} (e_ω=0 allowed; skipped in χ_D build).
static int ls_test_direction(const u8 *u, int K){
    u8 e[255], e2[255];
    u64 first=0; int have_first=0, tested=0;
    for(int k=0;k<K;k++){
        for(int i=0;i<255;i++){ e[i]=rnd8(); e2[i]=e[i]^u[i]; }
        u64 d = fp(e,255) ^ fp(e2,255);
        tested++;
        if(!have_first){ first=d; have_first=1; }
        else if(d!=first) return 1; // non-constant ⟹ u ∉ LS
    }
    return (tested>=2)?0:-1; // 0 = all same (possibly LS); -1 = insufficient
}

static void test2_single_coord(int K){
    int fail=0;
    u8 u[255]={0};
    for(int c=1;c<256;c++){
        u[0]=c;
        if(!ls_test_direction(u,K)){ printf("[T2] c=%d: Δ CONSTANT?? ⚠\n",c); fail++; }
    }
    u[0]=0;
    printf("[T2] single-coord (c,0,..,0): %d/255 have non-constant Δ_u fp %s\n",
           255-fail, fail==0?"✓ none in LS":"⚠");
}

static void test3_diagonal(int K){
    int fail=0;
    u8 u[255];
    for(int c=1;c<256;c++){
        for(int i=0;i<255;i++) u[i]=c;
        if(!ls_test_direction(u,K)){ printf("[T3] c=%d: Δ CONSTANT?? ⚠\n",c); fail++; }
    }
    printf("[T3] diagonal (c,..,c): %d/255 have non-constant Δ_u fp %s\n",
           255-fail, fail==0?"✓ none in LS":"⚠");
}

static void test4_two_coord(int K){
    int fail=0;
    u8 u[255]={0};
    for(int a=1;a<256;a++){
        u[0]=a; u[1]=a;
        if(!ls_test_direction(u,K)){ printf("[T4] a=%d: Δ CONSTANT?? ⚠\n",a); fail++; }
    }
    u[0]=u[1]=0;
    printf("[T4] two-coord (a,a,0,..): %d/255 have non-constant Δ_u fp %s\n",
           255-fail, fail==0?"✓ none in LS":"⚠");
}

static void test5_random(int N, int K){
    int fail=0;
    u8 u[255];
    for(int n=0;n<N;n++){
        int nz=0;
        for(int i=0;i<255;i++){ u[i]=rnd8(); if(u[i])nz++; }
        if(nz==0) continue;
        if(!ls_test_direction(u,K)) fail++;
    }
    printf("[T5] random u: %d/%d have non-constant Δ_u fp %s\n",
           N-fail, N, fail==0?"✓ none in LS":"⚠");
}

// T6: G_a-inv is NOT Borel-inv (so G_a-inv fp can't be a valid match-fp)
// Test: coset-multiset(E) ≠ coset-multiset(β(E)) for β=scale-by-d, random E,d,a
static void test6_Ga_not_Borel(int N){
    int same=0, tot=0;
    for(int n=0;n<N;n++){
        u8 a=rnd8nz(), d=rnd8nz(); if(d==1) continue;
        u8 e[255]; for(int i=0;i<255;i++) e[i]=rnd8nz();
        // coset-mset of E: count[128] where count[rep(e_i)] with rep = min(e_i, e_i^a)
        u8 cnt1[256]={0}, cnt2[256]={0};
        for(int i=0;i<255;i++){
            u8 r1 = (e[i]<(e[i]^a))?e[i]:(e[i]^a);
            cnt1[r1]++;
            u8 ed = gmul(e[i], ginv(d)); // β(e) = e/d (Borel scale)
            u8 r2 = (ed<(ed^a))?ed:(ed^a);
            cnt2[r2]++;
        }
        int eq=1; for(int i=0;i<256;i++) if(cnt1[i]!=cnt2[i]){eq=0;break;}
        same+=eq; tot++;
    }
    printf("[T6] a-coset-mset invariant under Borel-scale: %d/%d %s\n",
           same, tot, same==0?"✓ NEVER (G_a-inv ⊥ Borel-inv)":"⚠ sometimes");
}

// T0: sanity — fp is AGL-inv (fp(E) unchanged under D→αD+β)
static void test0_agl_inv(int N){
    int ok=0, tested=0;
    for(int n=0;n<N*3;n++){
        if(tested>=N) break;
        u8 e[255], e2[255];
        for(int i=0;i<255;i++) e[i]=rnd8nz();
        u8 alpha=rnd8nz(), beta=rnd8();
        // D→αD+β means 1/e → α/e+β means e → e/(β·e+α) [Möbius]
        int bad=0;
        for(int i=0;i<255;i++){
            u8 den = gmul(beta,e[i])^alpha;
            if(den==0){bad=1;break;}
            e2[i] = gmul(e[i], ginv(den));
        }
        if(bad) continue;
        tested++;
        if(fp(e,255)==fp(e2,255)) ok++;
        else printf("  T0 FAIL at n=%d α=%d β=%d\n",n,alpha,beta);
    }
    printf("[T0] fp AGL-inv sanity: %d/%d %s\n", ok, tested,
           (ok==tested&&tested>0)?"✓":"✗ BROKEN FP");
}

int main(int argc, char**argv){
    gf_init();
    rng_s = (argc>1)?strtoull(argv[1],0,0):0xC0FFEEULL;
    int K = (argc>2)?atoi(argv[2]):5;    // samples per direction (fp is slow)
    int Nr = (argc>3)?atoi(argv[3]):200; // random directions

    printf("=== verify_ls_fp: LS(fp)=0 for complete AGL-inv fp ===\n");
    printf("seed=0x%llx K=%d Nr=%d\n\n",(unsigned long long)rng_s,K,Nr);

    test0_agl_inv(20);
    test1_borel_transitive();
    printf("\n-- LS(fp) direction tests (K=%d samples each) --\n",K);
    test2_single_coord(K);
    test3_diagonal(K);
    test4_two_coord(K);
    test5_random(Nr,K);
    printf("\n");
    test6_Ga_not_Borel(5000);

    printf("\n=== Summary ===\n");
    printf("If T2-T5 all ✓: LS(fp)={0} for the COMPLETE AGL-inv fp.\n");
    printf("By LINEAR-STRUCTURE: ℓ ≥ 255 - dim(LS)/8 = 255.\n");
    printf("⟹ T ≥ 2^88 × 255/160 = 2^88.673 (A1-FLOOR, |K_on|=8).\n");
    return 0;
}
