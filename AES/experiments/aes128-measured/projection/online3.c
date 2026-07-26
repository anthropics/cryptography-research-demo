// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* online_bench.c -- REAL-SCALE microbenchmark of the ONLINE inner loop of the
 * 7-round AES-128 attack, for ONE candidate-pair context.
 *
 *   (A) OURS  : per candidate u_6[antidiag] guess: peel 255 ciphertext diffs
 *              (4 iSBOX + iMC-row-0 combine per omega), Linv lookup + field
 *              inversion -> 255-element sequence g_omega, complete chi-star
 *              fingerprint pair (fp0,fp1) with the v2 (brute-free) kernel,
 *              one hash-table probe per fp.
 *   (B) DFJ   : same peel, then an inner loop over all 256 values of
 *              s = u_5[0]: per s, iSBOX(v^s) to get the x_5[0] difference
 *              sequence, order-free multiset hash (sum of 64-bit random
 *              table entries), probe.  Early exit disabled for timing
 *              (both variants measured).
 *   Both with plain and DDT-Gray-amortized peel (4-bit BRGC over the u_6
 *   branch bits: cold 4 iSBOX/omega every 16 candidates, 1 iSBOX/omega on
 *   the 15 flips, as in ddt_gray_online.c).
 *
 * NOTE ON THE DELTA-SET: here we are COSTING the loop, not attacking.  The
 * delta-set of 256 chosen plaintexts is built "honestly" from a random base
 * plaintext and the TRUE k_{-1}[0,5,10,15] (diagonal) of a random key
 * (same construction as build_delta_set() in a1_e2e.c); the 256 ciphertexts
 * are REAL 7-round AES-128 encryptions (aes_core.h, NR=7).
 *
 * Round/byte naming (a1_e2e.c convention, rounds r=1..7, x_r/y_r/z_r/w_r):
 *   paper k_{-1}[0,5,10,15]  = CTX.rk[0][DIAG0]          (delta-set side)
 *   paper u_6[0,7,10,13]     = CTX.rk[7][ADIAG={0,13,10,7}] (last-round peel)
 *   paper u_5[0]  (kappa)    = (MC^{-1} rk[6])[0], absorbed by the Mobius bridge
 *   paper x_5[0]  (a_omega)  = internal x_6[0] in this indexing
 *   peeled value  v_omega    = MC^{-1}(x_7)[0] = S(a_omega) ^ kappa
 *   online seq    g_omega    = 1 / L^{-1}(v_omega ^ v_0)   (omega=1..255)
 *
 * Build: gcc -O3 -march=native -o online_bench online_bench.c
 * Usage: taskset -c <core> ./online_bench [seed] [logN_A] [logN_B] [reps]
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
#include "aes_core.h"
#include "kernel_v2_vec.h"   /* kernel_v2_vec + kernel_full_v2 (scalar) + tables */

#define NR 7

static aesctx CTX;
static u8 KEY[16];

/* ---------------- rng ---------------- */
static u64 rng_s[2];
static u64 rng(void){u64 a=rng_s[0],b=rng_s[1];rng_s[0]=b;a^=a<<23;rng_s[1]=a^b^(a>>18)^(b>>5);return rng_s[1]+b;}
static void rng_seed(u64 s){rng_s[0]=s^0x9e3779b97f4a7c15ULL;rng_s[1]=s^0xbf58476d1ce4e5b9ULL;for(int i=0;i<10;i++)rng();}

/* ---------------- timing ---------------- */
static inline u64 rdtsc_(void){unsigned a;__asm__ volatile("":::"memory");u64 t=__rdtscp(&a);__asm__ volatile("":::"memory");return t;}
static double calib_ratio(void){   /* core cycles per TSC tick via imul latency 3 */
    u64 best=~0ULL;
    for(int rep=0;rep<200;rep++){
        u64 y=3,t0=rdtsc_();
        for(long i=0;i<20000;i++){
            __asm__ volatile("imulq %1,%0\n\timulq %1,%0\n\timulq %1,%0\n\timulq %1,%0\n\timulq %1,%0":"+r"(y):"r"(5UL));
        }
        u64 t1=rdtsc_();
        __asm__ volatile(""::"r"(y));
        if(t1-t0<best)best=t1-t0;
    }
    return (3.0*100000.0)/(double)best;
}
/* force a value/array to be materialized (anti dead-code) */
#define ESCAPE(p) __asm__ volatile("" : : "g"(p) : "memory")
static volatile u64 SINK;

/* ---------------- tables ---------------- */
static u8  TT[4][256];        /* TT[i][x] = iMC[0][i] * iSBOX[x]  (fused peel table) */
static u64 R64[256];          /* random 64-bit per value, for the DFJ multiset hash */
#define HTBITS 20
static u64*HT;                /* dummy 2^20-entry open hash set (8 MB) */
static u8 ZERO256[256] __attribute__((aligned(64)));

static void tables_init(void){
    for(int i=0;i<4;i++) for(int x=0;x<256;x++) TT[i][x]=MUL[iMCc[0][i]][iSBOX[x]];
    for(int i=0;i<256;i++) R64[i]=rng();
    HT=aligned_alloc(64,sizeof(u64)<<HTBITS);
    for(u64 i=0;i<(1ULL<<HTBITS);i++) HT[i]=rng()|1;   /* random occupants: probes miss */
    memset(ZERO256,0,sizeof ZERO256);
}
static inline u64 ht_probe(u64 fp){
    return HT[(fp*0x9E3779B97F4A7C15ULL)>>(64-HTBITS)] == fp;
}

/* ---------------- GFNI helpers for the fully-vectorized variants ----------------
 * iSBOX(y) = inv( iAFF(y) ),  iAFF(y) = Linv(y) ^ Linv(0x63)   (affine, then inverse)
 *   -> vgf2p8affineqb (Linv matrix, c) followed by vgf2p8affineinvqb (identity).
 * Linv is GF(2)-linear -> one vgf2p8affineqb for the L^{-1} sequence step.
 * Both matrices are self-tested against the scalar tables at init.          */
#define IAFF_C 0x05          /* AES inverse-affine constant = Linv(0x63), asserted below */
static __m512i VM_LINV, VIDENT;
static void gfni_init(void){
    u8 col[8]; for(int j=0;j<8;j++) col[j]=Linv[1<<j];
    VM_LINV=_mm512_set1_epi64((long long)kv_mat_from_cols(col));
    VIDENT =_mm512_set1_epi64((long long)KV_IDENT);
    if(Linv[0x63]!=IAFF_C){fprintf(stderr,"IAFF_C mismatch\n");abort();}
    /* self-test: vectorized iSBOX and Linv equal the scalar tables on all 256 bytes */
    u8 in[256] __attribute__((aligned(64))), o1[256] __attribute__((aligned(64))), o2[256] __attribute__((aligned(64)));
    for(int x=0;x<256;x++) in[x]=(u8)x;
    for(int c=0;c<4;c++){
        __m512i y=_mm512_load_si512(in+64*c);
        __m512i t=_mm512_gf2p8affine_epi64_epi8(y,VM_LINV,0);
        _mm512_store_si512(o2+64*c,t);
        t=_mm512_gf2p8affine_epi64_epi8(y,VM_LINV,IAFF_C);
        t=_mm512_gf2p8affineinv_epi64_epi8(t,VIDENT,0);
        _mm512_store_si512(o1+64*c,t);
    }
    for(int x=0;x<256;x++){
        if(o1[x]!=iSBOX[x]||o2[x]!=Linv[x]){
            fprintf(stderr,"GFNI self-test FAIL at %d: isb %02x/%02x linv %02x/%02x\n",x,o1[x],iSBOX[x],o2[x],Linv[x]);
            abort();
        }
    }
}
static inline __m512i visbox512(__m512i y){
    __m512i t=_mm512_gf2p8affine_epi64_epi8(y,VM_LINV,IAFF_C);   /* iAFF(y) */
    return _mm512_gf2p8affineinv_epi64_epi8(t,VIDENT,0);        /* inv -> iSBOX(y) */
}

/* ---------------- online context: one candidate-pair delta-set ---------------- */
typedef struct {
    u8 P[16];            /* base plaintext */
    u8 Cd[256][16];      /* real ciphertexts, indexed by dv (difference at w_1[0]) */
    u8 Cs[4][256] __attribute__((aligned(64)));  /* SoA: Cs[i][dv]=Cd[dv][ADIAG[i]] */
    u8 a_tr[256];        /* internal x_6[0] (paper: x_5[0]) per element */
    u8 e_tr[256];        /* true offline sequence Delta x_6[0] (e_tr[0]=0) */
    u8 klast[4];         /* true last-round anti-diagonal bytes (paper: u_6[0,7,10,13]) */
    u8 kappa;            /* true absorbed byte (paper: u_5[0] = kappa) */
    u64 fD0,fD1;         /* chi-star fp of the TRUE offline sequence (vec kernel) */
    u64 H_tr;            /* DFJ multiset hash of the TRUE Delta x_6[0] sequence */
} octx_t;

static void octx_build(octx_t*O){
    for(int i=0;i<16;i++) O->P[i]=rng()&0xff;
    u8 km1[4]; for(int i=0;i<4;i++) km1[i]=CTX.rk[0][DIAG0[i]];
    /* w_1[col0] of the base plaintext, via the TRUE k_{-1} (we are costing,
     * not attacking: the delta-set is built honestly from P + k_{-1}) */
    u8 z1r[4]; for(int i=0;i<4;i++) z1r[i]=SBOX[O->P[DIAG0[i]]^km1[i]];
    u8 w1r[4]; for(int r=0;r<4;r++) w1r[r]=gmul(MCc[r][0],z1r[0])^gmul(MCc[r][1],z1r[1])^gmul(MCc[r][2],z1r[2])^gmul(MCc[r][3],z1r[3]);
    for(int dv=0;dv<256;dv++){
        u8 w1[4]={(u8)(w1r[0]^dv),w1r[1],w1r[2],w1r[3]};
        u8 z1[4]; for(int r=0;r<4;r++) z1[r]=gmul(iMCc[r][0],w1[0])^gmul(iMCc[r][1],w1[1])^gmul(iMCc[r][2],w1[2])^gmul(iMCc[r][3],w1[3]);
        u8 P[16]; memcpy(P,O->P,16);
        for(int i=0;i<4;i++) P[DIAG0[i]]=iSBOX[z1[i]]^km1[i];
        trace_t T; aes_enc_trace(&CTX,P,O->Cd[dv],&T);   /* REAL 7-round AES-128 */
        O->a_tr[dv]=T.x[6][0];
        for(int i=0;i<4;i++) O->Cs[i][dv]=O->Cd[dv][ADIAG[i]];
    }
    for(int dv=0;dv<256;dv++) O->e_tr[dv]=O->a_tr[dv]^O->a_tr[0];
    for(int i=0;i<4;i++) O->klast[i]=CTX.rk[7][ADIAG[i]];
    O->kappa=gmul(iMCc[0][0],CTX.rk[6][0])^gmul(iMCc[0][1],CTX.rk[6][1])
            ^gmul(iMCc[0][2],CTX.rk[6][2])^gmul(iMCc[0][3],CTX.rk[6][3]);
    /* truth fingerprints */
    kernel_v2_vec(O->e_tr,ZERO256,ZERO256,ZERO256,&O->fD0,&O->fD1);
    O->H_tr=0; for(int dv=1;dv<256;dv++) O->H_tr+=R64[(u8)(O->a_tr[dv]^O->a_tr[0])];
}

/* ---------------- the online peel ---------------- */
/* one element of the peel: v_w = MC^{-1}(x_7)[0] of ciphertext w under the
 * 4-byte last-round guess k (4 fused iSBOX*iMC lookups + 3 XOR) */
static inline u8 peel_elem(const octx_t*O,const u8 k[4],int w){
    return (u8)(TT[0][(u8)(O->Cs[0][w]^k[0])]^TT[1][(u8)(O->Cs[1][w]^k[1])]
               ^TT[2][(u8)(O->Cs[2][w]^k[2])]^TT[3][(u8)(O->Cs[3][w]^k[3])]);
}
/* plain: 4 fused (iSBOX * iMC-coef) lookups + 3 XOR per omega */
static inline void peel_plain(const octx_t*O,const u8 k[4],u8*v){
    const u8*C0=O->Cs[0],*C1=O->Cs[1],*C2=O->Cs[2],*C3=O->Cs[3];
    const u8 k0=k[0],k1=k[1],k2=k[2],k3=k[3];
    const u8*T0=TT[0],*T1=TT[1],*T2=TT[2],*T3=TT[3];
    for(int w=0;w<256;w++)
        v[w]=(u8)(T0[(u8)(C0[w]^k0)]^T1[(u8)(C1[w]^k1)]^T2[(u8)(C2[w]^k2)]^T3[(u8)(C3[w]^k3)]);
}
/* DDT-Gray amortized: cache the 4 per-position peeled columns; a candidate
 * differing in ONE branch byte costs 1 lookup/omega (ddt_gray_online.c). */
typedef struct { u8 pf[4][256]; u8 v[256]; u8 k[4]; } gray_t;
static inline void peel_gray_cold(const octx_t*O,gray_t*G,const u8 k[4]){
    for(int i=0;i<4;i++){
        const u8*C=O->Cs[i]; const u8*T=TT[i]; const u8 ki=k[i]; u8*p=G->pf[i];
        for(int w=0;w<256;w++) p[w]=T[(u8)(C[w]^ki)];
        G->k[i]=ki;
    }
    for(int w=0;w<256;w++) G->v[w]=(u8)(G->pf[0][w]^G->pf[1][w]^G->pf[2][w]^G->pf[3][w]);
}
static inline void peel_gray_flip(const octx_t*O,gray_t*G,int r,u8 knew){
    const u8*C=O->Cs[r]; const u8*T=TT[r]; u8*p=G->pf[r]; u8*v=G->v;
    for(int w=0;w<256;w++){ u8 np=T[(u8)(C[w]^knew)]; v[w]^=(u8)(np^p[w]); p[w]=np; }
    G->k[r]=knew;
}
/* sequence: e_omega = L^{-1}(v_omega ^ v_0)  (the kernel inverts it -> g) */
static inline void seq_e(const u8*v,u8*eH){
    const u8 v0=v[0]; eH[0]=0;
    for(int w=1;w<256;w++) eH[w]=Linv[(u8)(v[w]^v0)];
}
/* (B) DFJ per-s step: 256 iSBOX + 255 random-table adds */
static inline u64 dfj_hash_s(const u8*v,u8 s){
    u8 a0=iSBOX[(u8)(v[0]^s)];
    u64 h=0;
    for(int w=1;w<256;w++) h+=R64[(u8)(iSBOX[(u8)(v[w]^s)]^a0)];
    return h;
}

/* (B_published) DFJ AS PUBLISHED: for EVERY s the FULL peel is recomputed
 * per element (4 iSBOX + MC^-1 combine) plus the s-dependent iSBOX plus the
 * multiset insert: 255 x (4+1) S-box-class lookups + combine per s, no
 * sharing of the peel across the s-loop, no Gray, plain scalar, no early
 * exit.  This is the published accounting ell = 256*5 = 1280 lookups per
 * 9-byte (k_-1,u_6,u_5[0]) candidate. */
static inline u64 dfj_published_s(const octx_t*O,const u8 k[4],u8 s){
    const u8 a0=iSBOX[(u8)(peel_elem(O,k,0)^s)];
    u64 h=0;
    for(int w=1;w<256;w++) h+=R64[(u8)(iSBOX[(u8)(peel_elem(O,k,w)^s)]^a0)];
    return h;
}

/* -------- fully-vectorized (GFNI/AVX-512) variants of the online steps -------- */
/* peel, plain, fresh 4-byte guess: 4 positions x 4 ZMM: xor key, vec-iSBOX,
 * vec-mult by the iMC row-0 coefficient, xor-accumulate  (~0.2 cyc/omega) */
static inline void peel_vec(const octx_t*O,const u8 k[4],u8*v){
    __m512i acc[4];
    for(int i=0;i<4;i++){
        const __m512i kb=_mm512_set1_epi8((char)k[i]), mc=_mm512_set1_epi8((char)iMCc[0][i]);
        for(int c=0;c<4;c++){
            __m512i x=_mm512_xor_si512(_mm512_load_si512((const void*)(O->Cs[i]+64*c)),kb);
            __m512i t=_mm512_gf2p8mul_epi8(visbox512(x),mc);
            acc[c]=(i==0)?t:_mm512_xor_si512(acc[c],t);
        }
    }
    for(int c=0;c<4;c++) _mm512_store_si512((void*)(v+64*c),acc[c]);
}
/* sequence: e = L^{-1}(v ^ v_0) for all 256 lanes with one affine per ZMM */
static inline void seq_vec(const u8*v,u8*eH){
    const __m512i v0=_mm512_set1_epi8((char)v[0]);
    for(int c=0;c<4;c++){
        __m512i d=_mm512_xor_si512(_mm512_load_si512((const void*)(v+64*c)),v0);
        _mm512_store_si512((void*)(eH+64*c),_mm512_gf2p8affine_epi64_epi8(d,VM_LINV,0));
    }
}
/* (B) DFJ per-s step, vectorized iSBOX pass: a = iSBOX(v^s) in 8 vector ops,
 * then the order-free multiset hash sum of R64[Delta a] (inherently a gather) */
static inline u64 dfj_hash_s_vec(const u8*v,u8 s,u8*abuf){
    const __m512i sb=_mm512_set1_epi8((char)s);
    for(int c=0;c<4;c++){
        __m512i a=visbox512(_mm512_xor_si512(_mm512_load_si512((const void*)(v+64*c)),sb));
        _mm512_store_si512((void*)(abuf+64*c),a);
    }
    const u8 a0=abuf[0];
    u64 h=0;
    for(int w=1;w<256;w++) h+=R64[(u8)(abuf[w]^a0)];
    return h;
}

/* ---------------- per-candidate key stream ---------------- */
typedef struct {
    u64 N;
    u8 (*keys)[4];     /* plain: fresh 4 bytes per candidate; gray: group base keys */
    u8 (*dlt)[4];      /* gray: per-group nonzero branch deltas */
} kstream_t;
static void kstream_make(kstream_t*K,u64 N){
    K->N=N;
    K->keys=malloc(N*4); K->dlt=malloc(N*4);
    for(u64 i=0;i<N;i++){
        u64 r=rng();
        for(int j=0;j<4;j++){K->keys[i][j]=(u8)(r>>(8*j)); u8 d=(u8)(r>>(32+8*j)); K->dlt[i][j]=d?d:1;}
    }
}

/* ======================================================================
 * Benchmark bodies.  mode = cumulative stage cut:
 *   A: 1 peel | 2 +seq(L^-1) | 3 +chi* kernel (incl. field inversion) | 4 +2 probes
 *   B: 1 peel | 2 +s-loop iSBOX-only | 3 +multiset hash | 4 +probes (full)
 * ====================================================================== */
/* gray: 0 = plain peel (table), 1 = DDT-Gray amortized peel, 2 = GFNI vector peel+seq */
static u64 bench_A(const octx_t*O,const kstream_t*K,u64 N,int gray,int scalar,int mode){
    u8 v[256] __attribute__((aligned(64)));
    u8 eH[256] __attribute__((aligned(64)));
    gray_t G; u64 acc=0;
    for(u64 n=0;n<N;n++){
        const u8*vp;
        if(gray==0){ peel_plain(O,K->keys[n],v); vp=v; }
        else if(gray==2){ peel_vec(O,K->keys[n],v); vp=v; }
        else{
            u64 g=n>>4; unsigned st=(unsigned)(n&15);
            if(st==0) peel_gray_cold(O,&G,K->keys[g]);
            else{ int r=__builtin_ctz(st); peel_gray_flip(O,&G,r,(u8)(G.k[r]^K->dlt[g][r])); }
            vp=G.v;
        }
        ESCAPE(vp);
        if(mode>=2){ if(gray==2) seq_vec(vp,eH); else seq_e(vp,eH); ESCAPE(eH); }
        if(mode>=3){
            u64 f0,f1;
            if(!scalar) kernel_v2_vec(eH,ZERO256,ZERO256,ZERO256,&f0,&f1);
            else        kernel_full_v2(eH,ZERO256,ZERO256,ZERO256,&f0,&f1);
            acc^=f0^f1;
            if(mode>=4) acc+=ht_probe(f0)+(ht_probe(f1)<<1);
        }
    }
    return acc;
}
/* early: 0=none(all 256 s), 1=exit at a uniformly random planted s (models
 * the best case where a hit always occurs at a random position; in the real
 * attack hits are vanishingly rare, so the no-exit number is the honest one) */
static u64 bench_B(const octx_t*O,const kstream_t*K,u64 N,int gray,int early,int mode,double*iters_out){
    u8 v[256] __attribute__((aligned(64)));
    u8 abuf[256] __attribute__((aligned(64)));
    gray_t G; u64 acc=0; u64 totiter=0;
    for(u64 n=0;n<N;n++){
        const u8*vp;
        if(gray==0){ peel_plain(O,K->keys[n],v); vp=v; }
        else if(gray==2){ peel_vec(O,K->keys[n],v); vp=v; }
        else{
            u64 g=n>>4; unsigned st=(unsigned)(n&15);
            if(st==0) peel_gray_cold(O,&G,K->keys[g]);
            else{ int r=__builtin_ctz(st); peel_gray_flip(O,&G,r,(u8)(G.k[r]^K->dlt[g][r])); }
            vp=G.v;
        }
        ESCAPE(vp);
        if(mode<2) continue;
        int s_hit=early?(int)K->keys[n][0]:-1;
        int s;
        for(s=0;s<256;s++){
            if(mode==2){   /* sequence formation only: 256 iSBOX, byte checksum */
                u8 su=(u8)s, a0=iSBOX[(u8)(vp[0]^su)], c=0;
                for(int w=1;w<256;w++) c^=(u8)(iSBOX[(u8)(vp[w]^su)]^a0);
                acc+=c;
                continue;
            }
            u64 h=(gray==2)?dfj_hash_s_vec(vp,(u8)s,abuf):dfj_hash_s(vp,(u8)s);
            acc^=h;
            u64 hit=0;
            if(mode>=4){ hit=ht_probe(h); acc+=hit; }
            if(early && (hit || s==s_hit)){ s++; break; }
        }
        totiter+=(u64)s;
    }
    if(iters_out) *iters_out=(double)totiter/(double)N;
    return acc;
}

/* DFJ AS PUBLISHED: per candidate (k_-1,u_6) = 256 x [full 255-element peel
 * redone + s-step + multiset hash + probe], no sharing, no early exit.
 * mode (cumulative stage cut, for the phase breakdown):
 *   1: 256 x (255-element 4-iSBOX peel recomputed, byte checksum)
 *   2: + the s-dependent iSBOX per element (sequence a_omega), checksum
 *   3: + multiset hash (R64 sum)            4: + probe (full row)          */
static u64 bench_Bpub(const octx_t*O,const kstream_t*K,u64 N,int mode){
    u64 acc=0;
    for(u64 n=0;n<N;n++){
        const u8*k=K->keys[n];
        for(int s=0;s<256;s++){
            if(mode==1){ u8 c=0; for(int w=0;w<256;w++) c^=peel_elem(O,k,w); acc+=c; continue; }
            if(mode==2){
                u8 su=(u8)s, a0=iSBOX[(u8)(peel_elem(O,k,0)^su)], c=0;
                for(int w=1;w<256;w++) c^=(u8)(iSBOX[(u8)(peel_elem(O,k,w)^su)]^a0);
                acc+=c; continue;
            }
            u64 h=dfj_published_s(O,k,(u8)s);
            acc^=h;
            if(mode>=4) acc+=ht_probe(h);
        }
    }
    return acc;
}

/* ======================================================================
 * correctness gates
 * ====================================================================== */
static int gates(octx_t*O){
    int allok=1;
    printf("=== correctness gates ===\n");
    /* (i) known-answer with the TRUE u_6 bytes (and true s for DFJ) */
    {
        u8 v[256],eH[256];
        peel_plain(O,O->klast,v);
        seq_e(v,eH);
        u64 g0,g1,s0,s1;
        kernel_v2_vec(eH,ZERO256,ZERO256,ZERO256,&g0,&g1);
        kernel_full_v2(eH,ZERO256,ZERO256,ZERO256,&s0,&s1);
        int vs_eq=(g0==s0)&&(g1==s1);
        int either=(g0==O->fD0)||(g1==O->fD1);
        printf("[gate i-a] ours: online fp (vec)    fp0=%016llx fp1=%016llx\n",(unsigned long long)g0,(unsigned long long)g1);
        printf("[gate i-a] ours: online fp (scalar) fp0=%016llx fp1=%016llx  vec==scalar %s\n",(unsigned long long)s0,(unsigned long long)s1,vs_eq?"PASS":"FAIL");
        printf("[gate i-a] true offline seq fp      fD0=%016llx fD1=%016llx  EITHER-match %s\n",(unsigned long long)O->fD0,(unsigned long long)O->fD1,either?"PASS":"FAIL");
        if(!vs_eq||!either) allok=0;
        /* DFJ: true s = kappa -> a_omega must equal internal x_6[0] */
        int aok=1; for(int w=0;w<256;w++) if(iSBOX[(u8)(v[w]^O->kappa)]!=O->a_tr[w]){aok=0;break;}
        u64 H=dfj_hash_s(v,O->kappa);
        int hok=(H==O->H_tr);
        printf("[gate i-b] DFJ: iSBOX(v^kappa_true)==trace x_6[0] all 256: %s; multiset hash at true s == truth: %s (H=%016llx)\n",
               aok?"PASS":"FAIL",hok?"PASS":"FAIL",(unsigned long long)H);
        if(!aok||!hok) allok=0;
        /* DFJ-published (peel redone per s) at the true key and true s must
         * equal the shared-peel hash and the trace-derived truth */
        u64 Hpub=dfj_published_s(O,O->klast,O->kappa);
        int pok=(Hpub==H)&&(Hpub==O->H_tr);
        printf("[gate i-c] DFJ-published (per-s peel): hash at true s == shared-peel hash == truth: %s (Hpub=%016llx)\n",
               pok?"PASS":"FAIL",(unsigned long long)Hpub);
        if(!pok) allok=0;
        /* spot-check 1000 random (candidate,s): published per-s intermediate
         * peel == shared peel (elementwise), and per-s hash == shared-peel hash */
        int pv_ok=1;
        for(int t=0;t<1000&&pv_ok;t++){
            u8 kw[4]; for(int j=0;j<4;j++) kw[j]=rng()&0xff;
            u8 s=(u8)rng();
            u8 vs[256]; peel_plain(O,kw,vs);
            for(int w=0;w<256;w++) if(peel_elem(O,kw,w)!=vs[w]){pv_ok=0;break;}
            if(dfj_published_s(O,kw,s)!=dfj_hash_s(vs,s)) pv_ok=0;
        }
        printf("[gate i-d] DFJ-published per-s peel == shared peel and per-s hash == shared-hash, 1000 random (cand,s): %s\n",
               pv_ok?"PASS":"FAIL");
        if(!pv_ok) allok=0;
    }
    /* (ii) sequence elements vs internal trace, 20 random keys */
    {
        int seq_ok=0, a_ok=0, gray_ok=0, either_ok=0, vec_ok=0, NK=20;
        octx_t*T=aligned_alloc(64,(sizeof(octx_t)+63)&~63ULL);
        u8 savekey[16]; memcpy(savekey,KEY,16);
        for(int k=0;k<NK;k++){
            for(int i=0;i<16;i++) KEY[i]=rng()&0xff;
            aes_set_key(&CTX,KEY,NR);
            octx_build(T);
            u8 v[256],eH[256];
            peel_plain(T,T->klast,v); seq_e(v,eH);
            /* e_omega must equal inv(a_omega)^inv(a_0)  (a=x_6[0] from trace) */
            int ok=1; for(int w=1;w<256;w++) if(eH[w]!=(u8)(GF_inv[T->a_tr[w]]^GF_inv[T->a_tr[0]])){ok=0;break;}
            /* and v_omega must equal SBOX(a_omega)^kappa */
            for(int w=0;w<256;w++) if(v[w]!=(u8)(SBOX[T->a_tr[w]]^T->kappa)){ok=0;break;}
            seq_ok+=ok;
            int aok=1; for(int w=0;w<256;w++) if(iSBOX[(u8)(v[w]^T->kappa)]!=T->a_tr[w]){aok=0;break;}
            a_ok+=aok;
            /* gray walk ends at the true key: must reproduce the plain peel */
            gray_t G; u8 kk[4];
            for(int i=0;i<4;i++) kk[i]=(u8)(T->klast[i]^(u8)(rng()|1));
            peel_gray_cold(T,&G,kk);
            for(int r=0;r<4;r++) peel_gray_flip(T,&G,r,T->klast[r]);
            gray_ok+=(memcmp(G.v,v,256)==0);
            u64 g0,g1; kernel_v2_vec(eH,ZERO256,ZERO256,ZERO256,&g0,&g1);
            either_ok+=((g0==T->fD0)||(g1==T->fD1));
            /* GFNI-vectorized peel / seq / DFJ-hash must equal the scalar ones */
            u8 vv[256] __attribute__((aligned(64))), ev[256] __attribute__((aligned(64))), ab[256] __attribute__((aligned(64)));
            u8 va[256] __attribute__((aligned(64))), eHa[256] __attribute__((aligned(64)));
            memcpy(va,v,256); memcpy(eHa,eH,256);
            peel_vec(T,T->klast,vv); seq_vec(va,ev);
            int vok=(memcmp(vv,v,256)==0)&&(memcmp(ev,eHa,256)==0);
            for(int s=0;s<256&&vok;s+=51) if(dfj_hash_s(va,(u8)s)!=dfj_hash_s_vec(va,(u8)s,ab)) vok=0;
            vec_ok+=vok;
        }
        free(T); memcpy(KEY,savekey,16); aes_set_key(&CTX,KEY,NR);
        printf("[gate ii] 20 keys: e-seq==trace(inv a) %d/20, DFJ a-seq==trace %d/20, gray-peel==plain %d/20, GFNI-vec==scalar %d/20, chi* EITHER %d/20 -> %s\n",
               seq_ok,a_ok,gray_ok,vec_ok,either_ok,(seq_ok==NK&&a_ok==NK&&gray_ok==NK&&vec_ok==NK)?"PASS":"FAIL");
        /* (EITHER has a known ~0.4% degenerate-miss rate (a=0 / parity edge),
         *  reported but not a harness failure unless grossly off) */
        if(seq_ok!=NK||a_ok!=NK||gray_ok!=NK||vec_ok!=NK) allok=0;
        if(either_ok<18){ allok=0; printf("[gate ii] WARNING: EITHER rate too low\n"); }
    }
    /* (iii) fingerprints vary under WRONG u_6 guesses (no constant output) */
    {
        u8 v[256],eH[256];
        peel_plain(O,O->klast,v); seq_e(v,eH);
        u64 t0,t1; kernel_v2_vec(eH,ZERO256,ZERO256,ZERO256,&t0,&t1);
        u64 Htrue=dfj_hash_s(v,O->kappa);
        int NW=64, eq_true=0, h_eq=0;
        u64 seen[64]; int nd=0;
        for(int i=0;i<NW;i++){
            u8 kw[4]; for(int j=0;j<4;j++) kw[j]=rng()&0xff;
            if(!memcmp(kw,O->klast,4)) kw[0]^=1;
            peel_plain(O,kw,v); seq_e(v,eH);
            u64 f0,f1; kernel_v2_vec(eH,ZERO256,ZERO256,ZERO256,&f0,&f1);
            if(f0==t0||f1==t1) eq_true++;
            int dup=0; for(int j=0;j<nd;j++) if(seen[j]==f0){dup=1;break;}
            if(!dup) seen[nd++]=f0;
            if(dfj_hash_s(v,O->kappa)==Htrue) h_eq++;
        }
        int ok=(eq_true==0)&&(nd>=60)&&(h_eq==0);
        printf("[gate iii] 64 wrong u_6: collide-with-true %d, distinct fp0 %d/64, DFJ hash collide %d -> %s\n",
               eq_true,nd,h_eq,ok?"PASS":"FAIL");
        if(!ok) allok=0;
    }
    printf("=== gates: %s ===\n\n",allok?"ALL PASS":"FAIL");
    return allok;
}

/* ======================================================================
 * driver
 * ====================================================================== */
static int cmpd(const void*a,const void*b){double x=*(const double*)a,y=*(const double*)b;return x<y?-1:x>y;}
static inline double cpu_ns(void){struct timespec ts;clock_gettime(CLOCK_THREAD_CPUTIME_ID,&ts);return (double)ts.tv_sec*1e9+(double)ts.tv_nsec;}
/* stat: med/mn/mx in TSC ticks per candidate; cyc = median core cycles per
 * candidate via the core/TSC ratio calibrated immediately around THIS
 * configuration; cpuns = thread-CPU-time ns per candidate over all reps
 * (EXCLUDES time the thread is descheduled by other processes), cpucyc = cpuns *
 * 2.7 GHz * ratio = core cycles from thread CPU time. */
typedef struct { double med,mn,mx,cyc,ratio,cpuns,cpucyc; } stat_t;

static stat_t runA(const octx_t*O,const kstream_t*K,u64 N,int gray,int scalar,int mode,int reps,const char*name){
    double t[16];
    (void)bench_A(O,K,N>>4,gray,scalar,mode);           /* warm-up (untimed) */
    double r0=calib_ratio();
    double c0=cpu_ns();
    for(int r=0;r<reps;r++){
        u64 t0=rdtsc_();
        u64 acc=bench_A(O,K,N,gray,scalar,mode);
        u64 t1=rdtsc_();
        SINK^=acc;
        t[r]=(double)(t1-t0)/(double)N;
    }
    double c1=cpu_ns();
    double r1=calib_ratio();
    qsort(t,reps,sizeof(double),cmpd);
    stat_t S={t[reps/2],t[0],t[reps-1],0,(r0+r1)/2.0,0,0}; S.cyc=S.med*S.ratio; S.cpuns=(c1-c0)/((double)reps*(double)N); S.cpucyc=S.cpuns*2.7*S.ratio;
    printf("%-34s N=2^%d reps=%d  ticks/cand median %9.1f [min %9.1f max %9.1f] | core-cyc/cand median %9.1f (ratio %.2f/%.2f) | cpu-time %9.1f cyc (%.0f ns) | cyc/omega %7.3f\n",
           name,(int)__builtin_ctzll(N),reps,S.med,S.mn,S.mx,S.cyc,r0,r1,S.cpucyc,S.cpuns,S.cyc/255.0);
    fflush(stdout);
    return S;
}
static stat_t runB(const octx_t*O,const kstream_t*K,u64 N,int gray,int early,int mode,int reps,const char*name){
    double t[16],iters=0;
    (void)bench_B(O,K,N>>4,gray,early,mode,&iters);     /* warm-up (untimed) */
    double r0=calib_ratio();
    double c0=cpu_ns();
    for(int r=0;r<reps;r++){
        u64 t0=rdtsc_();
        u64 acc=bench_B(O,K,N,gray,early,mode,&iters);
        u64 t1=rdtsc_();
        SINK^=acc;
        t[r]=(double)(t1-t0)/(double)N;
    }
    double c1=cpu_ns();
    double r1=calib_ratio();
    qsort(t,reps,sizeof(double),cmpd);
    stat_t S={t[reps/2],t[0],t[reps-1],0,(r0+r1)/2.0,0,0}; S.cyc=S.med*S.ratio; S.cpuns=(c1-c0)/((double)reps*(double)N); S.cpucyc=S.cpuns*2.7*S.ratio;
    printf("%-34s N=2^%d reps=%d  ticks/cand median %9.1f [min %9.1f max %9.1f] | core-cyc/cand median %9.1f (ratio %.2f/%.2f) | cpu-time %9.1f cyc (%.0f ns) | s-iters/cand %.1f (cyc/s-iter %.1f)\n",
           name,(int)__builtin_ctzll(N),reps,S.med,S.mn,S.mx,S.cyc,r0,r1,S.cpucyc,S.cpuns,iters,iters>0?S.cyc/iters:0.0);
    fflush(stdout);
    return S;
}
static stat_t runBpub(const octx_t*O,const kstream_t*K,u64 N,int mode,int reps,const char*name){
    double t[16];
    (void)bench_Bpub(O,K,N>>6,mode);                      /* warm-up (untimed) */
    double r0=calib_ratio();
    double c0=cpu_ns();
    for(int r=0;r<reps;r++){
        u64 t0=rdtsc_();
        u64 acc=bench_Bpub(O,K,N,mode);
        u64 t1=rdtsc_();
        SINK^=acc;
        t[r]=(double)(t1-t0)/(double)N;
    }
    double c1=cpu_ns();
    double r1=calib_ratio();
    qsort(t,reps,sizeof(double),cmpd);
    stat_t S={t[reps/2],t[0],t[reps-1],0,(r0+r1)/2.0,0,0}; S.cyc=S.med*S.ratio; S.cpuns=(c1-c0)/((double)reps*(double)N); S.cpucyc=S.cpuns*2.7*S.ratio;
    printf("%-34s N=2^%d reps=%d  ticks/cand median %9.1f [min %9.1f max %9.1f] | core-cyc/cand median %9.1f (ratio %.2f/%.2f) | cpu-time %9.1f cyc (%.0f ns) | cyc/s-iter %.1f\n",
           name,(int)__builtin_ctzll(N),reps,S.med,S.mn,S.mx,S.cyc,r0,r1,S.cpucyc,S.cpuns,S.cyc/256.0);
    fflush(stdout);
    return S;
}
static void loadavg(const char*tag){FILE*f=fopen("/proc/loadavg","r");if(f){char b[128];if(fgets(b,sizeof b,f))printf("loadavg(%s): %s",tag,b);fclose(f);}}

int main(int argc,char**argv){
    u64 seed=argc>1?strtoull(argv[1],0,0):0xA15BE2C4ULL;
    int logNA=argc>2?atoi(argv[2]):20;
    /* B: each candidate contains a 256-iteration s-loop, so 2^14 candidates
     * = 2^22 inner (u_6,u_5[0]) iterations = 2^30 sequence elements */
    int logNB=argc>3?atoi(argv[3]):14;
    int reps =argc>4?atoi(argv[4]):5;
    int logNP=argc>5?atoi(argv[5]):18;   /* DFJ-published candidates */
    rng_seed(seed);
    aes_core_init(); chifast_init(); kernel_v2_vec_init(); tables_init(); gfni_init();

    for(int i=0;i<16;i++) KEY[i]=rng()&0xff;
    aes_set_key(&CTX,KEY,NR);
    printf("online_bench: seed=0x%llx  key=",(unsigned long long)seed);
    for(int i=0;i<16;i++) printf("%02x",KEY[i]);
    printf("  NR=%d  HT=2^%d entries (%zu MB)\n",NR,HTBITS,(size_t)(sizeof(u64)<<HTBITS)>>20);

    octx_t*O=aligned_alloc(64,(sizeof(octx_t)+63)&~63ULL);
    octx_build(O);
    printf("context: a_ref=x_6[0]=%02x kappa(u_5[0])=%02x klast(u_6 adiag)=%02x%02x%02x%02x\n\n",
           O->a_tr[0],O->kappa,O->klast[0],O->klast[1],O->klast[2],O->klast[3]);

    int gates_ok=gates(O);

    loadavg("start");
    printf("core/TSC ratio now (imul calibration): %.3f  [ticks = rdtscp TSC ticks; core-cyc = ticks*ratio, ratio re-calibrated around every configuration]\n\n",calib_ratio());

    u64 NA=1ULL<<logNA, NB=1ULL<<logNB, NP=1ULL<<logNP;
    kstream_t KA,KB,KP; kstream_make(&KA,NA); kstream_make(&KB,NB); kstream_make(&KP,NP);

    printf("=== (A) OURS: per-candidate = peel + L^-1 + chi-star fp (v2, brute-free) + 2 probes ===\n");
    stat_t Avp=runA(O,&KA,NA,0,0,4,reps,"A-plain-vec   (kernel_v2_vec)");
    stat_t Avg=runA(O,&KA,NA,1,0,4,reps,"A-gray-vec    (kernel_v2_vec)");
    stat_t Asp=runA(O,&KA,NA,0,1,4,reps,"A-plain-scalar(kernel_full_v2)");
    stat_t Asg=runA(O,&KA,NA,1,1,4,reps,"A-gray-scalar (kernel_full_v2)");
    printf("-- A breakdown (cumulative stages, plain peel) --\n");
    stat_t A1=runA(O,&KA,NA,0,0,1,reps,"A.peel only");
    stat_t A2=runA(O,&KA,NA,0,0,2,reps,"A.peel+seq(L^-1)");
    stat_t A3=runA(O,&KA,NA,0,0,3,reps,"A.peel+seq+kernel_v2_vec");
    stat_t A3s=runA(O,&KA,NA,0,1,3,reps,"A.peel+seq+kernel_full_v2(scal)");
    printf("-- A breakdown (gray peel) --\n");
    stat_t AG1=runA(O,&KA,NA,1,0,1,reps,"A.gray peel only");
    printf("-- A fully vectorized (GFNI peel + GFNI L^-1 + kernel_v2_vec + probes) --\n");
    stat_t Avv=runA(O,&KA,NA,2,0,4,reps,"A-allvec (GFNI peel)");
    stat_t AV1=runA(O,&KA,NA,2,0,1,reps,"A.vec peel only");
    stat_t AV3=runA(O,&KA,NA,2,0,3,reps,"A.vec peel+seq+kernel (no probe)");

    printf("\n=== (B) DFJ-style: per-candidate = peel + 256 x (iSBOX pass + multiset hash + probe) ===\n");
    stat_t Bp =runB(O,&KB,NB,0,0,4,reps,"B-plain no-early-exit");
    stat_t Bg =runB(O,&KB,NB,1,0,4,reps,"B-gray  no-early-exit");
    stat_t Bpe=runB(O,&KB,NB,0,1,4,reps,"B-plain EARLY-exit@planted-s");
    printf("-- B breakdown (cumulative stages, plain peel, no early exit) --\n");
    stat_t B1=runB(O,&KB,NB,0,0,1,reps,"B.peel only");
    stat_t B2=runB(O,&KB,NB,0,0,2,reps,"B.peel+256x iSBOX-pass");
    stat_t B3=runB(O,&KB,NB,0,0,3,reps,"B.peel+256x(iSBOX+hash)");
    printf("-- B vectorized (GFNI peel + per-s GFNI iSBOX pass; R64 multiset sum + probe stay gathers) --\n");
    stat_t Bv =runB(O,&KB,NB,2,0,4,reps,"B-allvec no-early-exit");
    stat_t Bv3=runB(O,&KB,NB,2,0,3,reps,"B.vec peel+256x(vecISBOX+hash)");

    printf("\n=== (B_published) DFJ AS PUBLISHED: per (k_-1,u_6) candidate = 256 x (full 255-element peel redone + iSBOX(v^s) + multiset hash + probe) ===\n");
    stat_t Bpub=runBpub(O,&KP,NP,4,reps,"DFJ-published plain no-exit");
    printf("-- B_published breakdown (cumulative stages) --\n");
    stat_t P1=runBpub(O,&KP,NP,1,reps,"Bpub.256x(255-elt peel redone)");
    stat_t P2=runBpub(O,&KP,NP,2,reps,"Bpub.256x(peel+iSBOX seq)");
    stat_t P3=runBpub(O,&KP,NP,3,reps,"Bpub.256x(peel+seq+hash)");
    loadavg("end");

    /* ---------------- THREE-WAY table ---------------- */
    printf("\n=== THREE-WAY COMPARISON (per (k_{-1},u_6) candidate; cyc = calibrated core cycles, median of 5 reps; [min,max] ticks) ===\n");
    printf("%-38s %12s %12s  %-24s %16s %16s\n","kernel","core-cyc(TSC)","cyc(cputime)","[min,max] ticks","vs ours-allvec","vs ours-plain");
    #define ROW(nm,S) printf("%-38s %12.0f %12.0f  [%9.0f,%9.0f] %14.2fx %15.2fx\n",nm,(S).cyc,(S).cpucyc,(S).mn,(S).mx,(S).cpucyc/Avv.cpucyc,(S).cpucyc/Avp.cpucyc)
    ROW("DFJ-published (plain)",Bpub);
    ROW("DFJ-optimized plain (shared peel)",Bp);
    ROW("DFJ-optimized gray  (shared+Gray)",Bg);
    ROW("DFJ-optimized allvec (GFNI peel+s)",Bv);
    ROW("ours plain  (vec chi* kernel)",Avp);
    ROW("ours gray   (vec chi* kernel)",Avg);
    ROW("ours allvec (GFNI peel+vec chi*)",Avv);
    #undef ROW
    printf("(ratios from thread-CPU-time cycles; TSC-based in brackets)\n");
    printf("HEADLINE  ours-allvec vs DFJ-optimized-allvec (both optimized):  %.2fx (TSC %.2fx)  [%.0f / %.0f cyc]\n",
           Bv.cpucyc/Avv.cpucyc,Bv.cyc/Avv.cyc,Bv.cpucyc,Avv.cpucyc);
    printf("HEADLINE  ours-allvec vs DFJ-published-plain (ours vs as-published): %.2fx (TSC %.2fx)  [%.0f / %.0f cyc]\n",
           Bpub.cpucyc/Avv.cpucyc,Bpub.cyc/Avv.cyc,Bpub.cpucyc,Avv.cpucyc);
    printf("          ours-plain  vs DFJ-published-plain (same impl level):     %.2fx (TSC %.2fx)  [%.0f / %.0f cyc]\n",
           Bpub.cpucyc/Avp.cpucyc,Bpub.cyc/Avp.cyc,Bpub.cpucyc,Avp.cpucyc);
    printf("          ours-allvec vs DFJ-optimized-plain:                       %.2fx (TSC %.2fx)  [%.0f / %.0f cyc]\n",
           Bp.cpucyc/Avv.cpucyc,Bp.cyc/Avv.cyc,Bp.cpucyc,Avv.cpucyc);

    printf("\n=== SUMMARY (per candidate; [ticks = TSC 2.7GHz reference ticks | cyc = calibrated core cycles]) ===\n");
    printf("A-plain vec %9.1f t %9.1f cyc | A-gray vec %9.1f t %9.1f cyc | A-plain scalar %9.1f t %9.1f cyc | A-gray scalar %9.1f t %9.1f cyc | A-allvec %9.1f t %9.1f cyc\n",
           Avp.med,Avp.cyc,Avg.med,Avg.cyc,Asp.med,Asp.cyc,Asg.med,Asg.cyc,Avv.med,Avv.cyc);
    printf("B-plain     %9.1f t %9.1f cyc | B-gray     %9.1f t %9.1f cyc | B-plain early-exit %9.1f t %9.1f cyc | B-allvec %9.1f t %9.1f cyc\n",
           Bp.med,Bp.cyc,Bg.med,Bg.cyc,Bpe.med,Bpe.cyc,Bv.med,Bv.cyc);
    printf("A breakdown (core-cyc): peel %.1f | +seq(L^-1) %.1f | +chi* vec %.1f (scalar chi* %.1f) | +2 probes %.1f | (GFNI peel %.1f, GFNI peel+seq+kernel %.1f)\n",
           A1.cyc,A2.cyc-A1.cyc,A3.cyc-A2.cyc,A3s.cyc-A2.cyc,Avp.cyc-A3.cyc,AV1.cyc,AV3.cyc);
    printf("B breakdown (core-cyc): peel %.1f | 256x(iSBOX pass+hash) %.1f [= %.1f /s; iSBOX-only stage %.1f/s, +hash %.1f/s] | 256x probe %.1f [= %.1f /s]\n",
           B1.cyc,B3.cyc-B1.cyc,(B3.cyc-B1.cyc)/256.0,(B2.cyc-B1.cyc)/256.0,(B3.cyc-B2.cyc)/256.0,Bp.cyc-B3.cyc,(Bp.cyc-B3.cyc)/256.0);
    printf("RATIO B/A plain (vec chi*):        ticks %.2f x | core-cyc %.2f x   [%.0f / %.0f cyc]\n",Bp.med/Avp.med,Bp.cyc/Avp.cyc,Bp.cyc,Avp.cyc);
    printf("RATIO B/A gray  (vec chi*):        ticks %.2f x | core-cyc %.2f x   [%.0f / %.0f cyc]\n",Bg.med/Avg.med,Bg.cyc/Avg.cyc,Bg.cyc,Avg.cyc);
    printf("RATIO B/A plain (scalar chi*):     ticks %.2f x | core-cyc %.2f x   [%.0f / %.0f cyc]\n",Bp.med/Asp.med,Bp.cyc/Asp.cyc,Bp.cyc,Asp.cyc);
    printf("RATIO B/A allvec (both GFNI):      ticks %.2f x | core-cyc %.2f x   [%.0f / %.0f cyc]\n",Bv.med/Avv.med,Bv.cyc/Avv.cyc,Bv.cyc,Avv.cyc);
    printf("RATIO B/A compute-only (no probes): ticks %.2f x | core-cyc %.2f x   [B3 %.0f / A3 %.0f cyc];  vec-vec %.2f x\n",
           B3.med/A3.med,B3.cyc/A3.cyc,B3.cyc,A3.cyc,Bv3.cyc/AV3.cyc);
    printf("MODEL-FORM: 256 * (per-u5-iteration cost) / (our fingerprint cost) = 256 * %.1f / %.1f = %.1f  (core-cyc; fingerprint cost = A minus its peel)\n",
           (Bp.cyc-B1.cyc)/256.0,Avp.cyc-A1.cyc,(Bp.cyc-B1.cyc)/(Avp.cyc-A1.cyc));

    /* ---------------- PER-PHASE BREAKDOWN, all three rows ----------------
     * cycles per (k_-1,u_6) candidate and per sequence element; A processes
     * 255 elements per candidate; B/C process 256 x 255 elements (s-loop). */
    {
        double Aelem=255.0, Belem=256.0*255.0;
        double A_peel=A1.cyc, A_seq=A2.cyc-A1.cyc, A_fp=A3.cyc-A2.cyc, A_pr=Avp.cyc-A3.cyc;
        double B_peel=B1.cyc, B_iS=B2.cyc-B1.cyc, B_ha=B3.cyc-B2.cyc, B_pr=Bp.cyc-B3.cyc;
        double C_peel=P1.cyc, C_iS=P2.cyc-P1.cyc, C_ha=P3.cyc-P2.cyc, C_pr=Bpub.cyc-P3.cyc;
        printf("\n=== PER-PHASE BREAKDOWN per (k_-1,u_6) candidate, core cycles (median stage deltas) ===\n");
        printf("%-36s %12s %12s %16s %10s %12s\n","variant","peel","seq-build","fp/hash+s-loop","probes","total");
        printf("%-36s %12.1f %12.1f %16.1f %10.1f %12.1f\n","(A) OURS plain-peel vec-chi*",A_peel,A_seq,A_fp,A_pr,Avp.cyc);
        printf("%-36s %12s %12s %16s %10s %12.1f   [GFNI peel %.1f, GFNI peel+seq+kernel %.1f, +2 probes %.1f]\n",
               "(A) OURS allvec (GFNI)","-","-","-","-",Avv.cyc,AV1.cyc,AV3.cyc,Avv.cyc-AV3.cyc);
        printf("%-36s %12.1f %12.1f %16.1f %10.1f %12.1f   [peel shared once; 256-step s-loop: iSBOX pass + R64 hash]\n",
               "(B) DFJ-optimized plain",B_peel,B_iS,B_ha,B_pr,Bp.cyc);
        printf("%-36s %12s %12s %16s %10s %12.1f   [GFNI peel + per-s GFNI iSBOX; gather hash]\n",
               "(B) DFJ-optimized allvec","-","-","-","-",Bv.cyc);
        printf("%-36s %12.1f %12.1f %16.1f %10.1f %12.1f   [peel RECOMPUTED inside every s-iteration]\n",
               "(C) DFJ-as-published plain",C_peel,C_iS,C_ha,C_pr,Bpub.cyc);
        printf("per sequence element (A: /255; B,C: /(256x255)): A peel %.3f seq %.3f fp %.3f probes %.3f | B peel %.4f iSBOX %.3f hash %.3f probes %.3f | C peel %.3f iSBOX %.3f hash %.3f probes %.3f\n",
               A_peel/Aelem,A_seq/Aelem,A_fp/Aelem,A_pr/Aelem,
               B_peel/Belem,B_iS/Belem,B_ha/Belem,B_pr/Belem,
               C_peel/Belem,C_iS/Belem,C_ha/Belem,C_pr/Belem);
        printf("per s-iteration (B,C; /256): B %.1f cyc [iSBOX %.1f + hash %.1f + probe %.1f] | C %.1f cyc [peel %.1f + iSBOX %.1f + hash %.1f + probe %.1f]\n",
               (Bp.cyc-B1.cyc)/256.0,B_iS/256.0,B_ha/256.0,B_pr/256.0,
               Bpub.cyc/256.0,C_peel/256.0,C_iS/256.0,C_ha/256.0,C_pr/256.0);
    }

    /* ------------------------------------------------------------------
     * HEAD-TO-HEAD, interleaved round-robin: each round times one A chunk and
     * one B chunk (same config pair) back-to-back, so the per-round ratio is
     * immune to core-clock drift (no PMU in this VM to count real cycles).
     * ------------------------------------------------------------------ */
    printf("\n=== HEAD-TO-HEAD interleaved (per-round ratio, %d rounds; A chunk 2^%d, B chunk 2^%d, Bpub chunk 2^%d cands) ===\n",21,logNA-4,logNB-4,logNP-10);
    {
        enum{RR=21}; u64 nA=NA>>4, nB=NB>>4, nP=NP>>10;
        double rp[RR],rg[RR],rv[RR],rpb[RR],rpv[RR],cal[RR],aP[RR],bP[RR],aG[RR],bG[RR],aV[RR],bV[RR],pB[RR];
        double it;
        for(int r=0;r<RR;r++){
            cal[r]=calib_ratio();
            u64 t0,t1;
            t0=rdtsc_(); SINK^=bench_A(O,&KA,nA,0,0,4); t1=rdtsc_(); aP[r]=(double)(t1-t0)/nA;
            t0=rdtsc_(); SINK^=bench_B(O,&KB,nB,0,0,4,&it); t1=rdtsc_(); bP[r]=(double)(t1-t0)/nB;
            t0=rdtsc_(); SINK^=bench_A(O,&KA,nA,1,0,4); t1=rdtsc_(); aG[r]=(double)(t1-t0)/nA;
            t0=rdtsc_(); SINK^=bench_B(O,&KB,nB,1,0,4,&it); t1=rdtsc_(); bG[r]=(double)(t1-t0)/nB;
            t0=rdtsc_(); SINK^=bench_A(O,&KA,nA,2,0,4); t1=rdtsc_(); aV[r]=(double)(t1-t0)/nA;
            t0=rdtsc_(); SINK^=bench_B(O,&KB,nB,2,0,4,&it); t1=rdtsc_(); bV[r]=(double)(t1-t0)/nB;
            t0=rdtsc_(); SINK^=bench_Bpub(O,&KP,nP,4); t1=rdtsc_(); pB[r]=(double)(t1-t0)/nP;
            rp[r]=bP[r]/aP[r]; rg[r]=bG[r]/aG[r]; rv[r]=bV[r]/aV[r]; rpb[r]=pB[r]/aP[r]; rpv[r]=pB[r]/aV[r];
        }
        qsort(rp,RR,sizeof(double),cmpd); qsort(rg,RR,sizeof(double),cmpd); qsort(rv,RR,sizeof(double),cmpd);
        qsort(rpb,RR,sizeof(double),cmpd); qsort(rpv,RR,sizeof(double),cmpd);
        qsort(aP,RR,sizeof(double),cmpd); qsort(bP,RR,sizeof(double),cmpd);
        qsort(aG,RR,sizeof(double),cmpd); qsort(bG,RR,sizeof(double),cmpd);
        qsort(aV,RR,sizeof(double),cmpd); qsort(bV,RR,sizeof(double),cmpd);
        qsort(pB,RR,sizeof(double),cmpd);
        qsort(cal,RR,sizeof(double),cmpd);
        double c=cal[RR/2];
        printf("core/TSC ratio over rounds: median %.3f (min %.3f, max %.3f)\n",c,cal[0],cal[RR-1]);
        printf("plain : A %8.1f t (%8.1f cyc)  B %9.1f t (%9.1f cyc)  ratio median %6.2fx [min %.2f max %.2f]\n",
               aP[RR/2],aP[RR/2]*c,bP[RR/2],bP[RR/2]*c,rp[RR/2],rp[0],rp[RR-1]);
        printf("gray  : A %8.1f t (%8.1f cyc)  B %9.1f t (%9.1f cyc)  ratio median %6.2fx [min %.2f max %.2f]\n",
               aG[RR/2],aG[RR/2]*c,bG[RR/2],bG[RR/2]*c,rg[RR/2],rg[0],rg[RR-1]);
        printf("allvec: A %8.1f t (%8.1f cyc)  B %9.1f t (%9.1f cyc)  ratio median %6.2fx [min %.2f max %.2f]\n",
               aV[RR/2],aV[RR/2]*c,bV[RR/2],bV[RR/2]*c,rv[RR/2],rv[0],rv[RR-1]);
        printf("published: Bpub %9.1f t (%9.1f cyc)  vs ours-plain ratio median %6.2fx [min %.2f max %.2f] | vs ours-allvec ratio median %6.2fx [min %.2f max %.2f]\n",
               pB[RR/2],pB[RR/2]*c,rpb[RR/2],rpb[0],rpb[RR-1],rpv[RR/2],rpv[0],rpv[RR-1]);
    }
    printf("gates: %s\n",gates_ok?"ALL PASS":"FAIL");
    (void)AG1;
    return gates_ok?0:1;
}
