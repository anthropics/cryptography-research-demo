// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef KERNEL_V2_VEC_H
#define KERNEL_V2_VEC_H
/* kernel_v2_vec: GFNI/AVX-512 vectorization of kernel_full_v2
 * (kernel_v2_scalar.h).  Bit-identical (fp0,fp1).
 *
 * Structure:
 *   1. e = U0^U1^U2^U3 (4 ZMM), d = inv(e) via VGF2P8AFFINEINVQB.
 *   2. ONE parity scatter of the raw multiset {d_1..d_255} -> chi_d (4
 *      qwords).  In-register class biasing (t ^ 64w) + VPSLLVQ (out-of-range
 *      counts give 0) + VPTERNLOG merging: 20 vector uops / 16 elements.
 *      The scatter does not depend on the scalar corner, so the corner and
 *      the beta-solve chains hide behind it.
 *   3. power sums S1,S3,S5,S7 (GF2P8MULB powers + byte xor-folds).
 *   4. scalar corner (ALL products via L1-resident log/exp tables with a
 *      zero zone -- no 64 KB MUL[][] traffic): P7 -> a0,b0 (fp0) and the
 *      fp1 beta-solve a1, A=S3/S1^3, u=UROOT[A], {ba,bb} -- hoisted, so the
 *      tail is two short vector shifts + two FNV chains.
 *   5. chi0 = gather(chi_d, u -> a0^-1(u ^ b0)) via VPERMB; fp0 = FNV(chi0).
 *      chi1 = gather(chi_d, u -> a1^-1 u) ^ e_0;
 *      fp1 = min(FNV(chi1 shifted by ba), FNV(chi1 shifted by bb)).
 *   alpha-chain for fp1: P7' -> P11 -> scalar deep chain (P13,P19,P23,...,
 *   P61) -> exhaustive (measured: P7'=P11=0 => P13=0 on this distribution).
 *   Rare tiers (S1=0 at rate 2^-8; deep chain at rate ~2^-16)
 *   call the scalar Designer-1 canonicalizer on the explicit set.
 */
#include <immintrin.h>
#include "kernel_v2_scalar.h"   /* reference + tables + canon_even */

/* ---------------------------------------------------------------- */
static const u64 KV_IDENT = 0x0102040810204080ULL;
static __m512i KV_ONE, KV_B64, KV_B128, KV_B192, KV_SHR3, KV_LUT, KV_IDX[4];
static __m256i KV_WPERM[4], KV_BCTL256[8];
static u64 KV_BITPERM[8];
/* L1-resident scalar GF tables */
static u16 KV_LOGZ[256];        /* log, with LOGZ[0]=512 -> zero zone */
static u8  KV_EXPZ[1040];       /* exp (period 255) for 0..509, zeros 510.. */
static u8  KV_SQ[256], KV_INVCUBE[256];
static u8  KV_A73[256], KV_A116[256]; /* P_m^(-1/m) maps */

static u64 kv_mat_from_cols(const u8 col[8]){
    u64 q=0;
    for(int i=0;i<8;i++){ u8 row=0; for(int j=0;j<8;j++) if((col[j]>>i)&1) row|=(u8)(1u<<j); q|=(u64)row<<(8*(7-i)); }
    return q;
}

static void kernel_v2_vec_init(void){
    kernel_v2_init();                  /* POW1357, POW9_11_13, canon_even tables */
    KV_ONE=_mm512_set1_epi64(1);
    KV_B64=_mm512_set1_epi64(64); KV_B128=_mm512_set1_epi64(128); KV_B192=_mm512_set1_epi64(192);
    { u8 col[8]; for(int j=0;j<8;j++) col[j]=(j>=3)?(u8)(1<<(j-3)):0; KV_SHR3=_mm512_set1_epi64((long long)kv_mat_from_cols(col)); }
    for(int l=0;l<8;l++){ u8 col[8]; for(int j=0;j<8;j++) col[j]=(u8)(1<<(j^l)); KV_BITPERM[l]=kv_mat_from_cols(col); }
    for(int h2=0;h2<4;h2++) KV_WPERM[h2]=_mm256_set_epi64x(3^h2,2^h2,1^h2,0^h2);
    for(int c=0;c<8;c++){ char b[32]; for(int j=0;j<32;j++) b[j]=(char)(((j&7)^c)|(j&8)); KV_BCTL256[c]=_mm256_loadu_si256((const __m256i*)b); }
    { u8 lut[16]={1,2,4,8,16,32,64,128,0,0,0,0,0,0,0,0}; KV_LUT=_mm512_broadcast_i32x4(_mm_loadu_si128((const __m128i*)lut)); }
    for(int c=0;c<4;c++){ u8 t[64]; for(int i=0;i<64;i++) t[i]=(u8)(64*c+i); KV_IDX[c]=_mm512_loadu_si512(t); }
    /* scalar tables */
    for(int x=0;x<256;x++) KV_LOGZ[x]= x ? GF_log[x] : 512;
    for(int r=0;r<1040;r++) KV_EXPZ[r]= (r<510) ? GF_exp[r%255] : 0;
    for(int x=0;x<256;x++){ KV_SQ[x]=gmul((u8)x,(u8)x); KV_INVCUBE[x]= x ? GF_inv[gmul((u8)x,gmul((u8)x,(u8)x))] : 0; }
    for(int p=0;p<256;p++){
        if(!p){KV_A73[0]=KV_A116[0]=0;continue;}
        int l=GF_log[p];
        KV_A73[p] =GF_exp[(255-(l*73)%255)%255];
        KV_A116[p]=GF_exp[(255-(l*116)%255)%255];
    }
}

/* ---------------- small helpers ---------------- */
static inline u8 kv_mulz(u8 a,u8 b){          /* branchless GF mult, L1 tables */
    return KV_EXPZ[KV_LOGZ[a]+KV_LOGZ[b]];
}
static inline u64 kv_fnv4(u64 w0,u64 w1,u64 w2,u64 w3){
    u64 h=0xcbf29ce484222325ULL;
    h^=w0;h*=0x100000001b3ULL; h^=w1;h*=0x100000001b3ULL;
    h^=w2;h*=0x100000001b3ULL; h^=w3;h*=0x100000001b3ULL;
    return h;
}
static inline u64 kv_hxor_q(__m512i v){
    __m256i x=_mm256_xor_si256(_mm512_castsi512_si256(v),_mm512_extracti64x4_epi64(v,1));
    __m128i y=_mm_xor_si128(_mm256_castsi256_si128(x),_mm256_extracti128_si256(x,1));
    y=_mm_xor_si128(y,_mm_unpackhi_epi64(y,y));
    return (u64)_mm_cvtsi128_si64(y);
}
static inline u8 kv_hxor_b(__m512i v){
    u64 q=kv_hxor_q(v); q^=q>>32; q^=q>>16; q^=q>>8; return (u8)q;
}
/* shift_beta(chi)[t] = chi[t ^ beta] on a YMM */
static inline __m256i kv_shift256v(__m256i chi, u8 beta){
    int h2=beta>>6, c=(beta>>3)&7, l=beta&7;
    __m256i z=_mm256_permutexvar_epi64(KV_WPERM[h2],chi);
    z=_mm256_shuffle_epi8(z,KV_BCTL256[c]);
    z=_mm256_gf2p8affine_epi64_epi8(z,_mm256_set1_epi64x((long long)KV_BITPERM[l]),0);
    return z;
}
static inline u64 kv_fnv256(__m256i z){
    __m128i lo=_mm256_castsi256_si128(z),hi=_mm256_extracti128_si256(z,1);
    return kv_fnv4((u64)_mm_cvtsi128_si64(lo),(u64)_mm_extract_epi64(lo,1),
                   (u64)_mm_cvtsi128_si64(hi),(u64)_mm_extract_epi64(hi,1));
}

/* ---------------- chi_d parity scatter (255 elements) ----------------
 * in-register class bias: for class w the count vector is t ^ 64w (valid
 * 0..63 exactly for elements of class w, >=64 otherwise -> sllv gives 0);
 * position 0 is forced to count 1024 (out of range in every class). */
static inline void kv_scatter(const u8*dbuf, u64 w[4]){
    __m512i A0=_mm512_setzero_si512(),A1=A0,A2=A0,A3=A0;
    const __m512i one=KV_ONE, b1=KV_B64, b2=KV_B128, b3=KV_B192;
    for(int g=0;g<32;g+=2){
        __m512i t=_mm512_cvtepu8_epi64(_mm_loadl_epi64((const __m128i*)(dbuf+8*g)));
        __m512i s=_mm512_cvtepu8_epi64(_mm_loadl_epi64((const __m128i*)(dbuf+8*g+8)));
        if(g==0) t=_mm512_mask_blend_epi64(1,t,_mm512_set1_epi64(1024));
        __m512i v0=_mm512_sllv_epi64(one,t),                      u0=_mm512_sllv_epi64(one,s);
        __m512i v1=_mm512_sllv_epi64(one,_mm512_xor_si512(t,b1)),u1=_mm512_sllv_epi64(one,_mm512_xor_si512(s,b1));
        __m512i v2=_mm512_sllv_epi64(one,_mm512_xor_si512(t,b2)),u2=_mm512_sllv_epi64(one,_mm512_xor_si512(s,b2));
        __m512i v3=_mm512_sllv_epi64(one,_mm512_xor_si512(t,b3)),u3=_mm512_sllv_epi64(one,_mm512_xor_si512(s,b3));
        A0=_mm512_ternarylogic_epi64(A0,v0,u0,0x96);
        A1=_mm512_ternarylogic_epi64(A1,v1,u1,0x96);
        A2=_mm512_ternarylogic_epi64(A2,v2,u2,0x96);
        A3=_mm512_ternarylogic_epi64(A3,v3,u3,0x96);
    }
    w[0]=kv_hxor_q(A0); w[1]=kv_hxor_q(A1); w[2]=kv_hxor_q(A2); w[3]=kv_hxor_q(A3);
}

/* chi'[u] = chi[ainv * (u ^ beta)]: VPERMB byte-gather from the 32-byte
 * image of chi in the low half of chiv + bit LUT + VPTESTMB. */
static inline void kv_gather(__m512i chiv, u8 ainv, u8 beta, u64 out[4]){
    const __m512i ba=_mm512_set1_epi8((char)ainv), bb=_mm512_set1_epi8((char)beta);
    const __m512i m7=_mm512_set1_epi8(7);
    for(int c=0;c<4;c++){
        __m512i f =_mm512_gf2p8mul_epi8(_mm512_xor_si512(KV_IDX[c],bb),ba);
        __m512i fh=_mm512_gf2p8affine_epi64_epi8(f,KV_SHR3,0);
        __m512i g =_mm512_permutexvar_epi8(fh,chiv);
        __m512i m =_mm512_shuffle_epi8(KV_LUT,_mm512_and_si512(f,m7));
        out[c]=_cvtmask64_u64(_mm512_test_epi8_mask(g,m));
    }
}

/* ---------------------------------------------------------------- */
static u64 kvv_stat_s1zero=0, kvv_stat_p7z=0, kvv_stat_deep=0;

static inline void kernel_v2_vec(const u8*U0,const u8*U1,const u8*U2,const u8*U3,
                                 u64*fp0,u64*fp1){
    /* 1. e, d=inv(e) */
    const __m512i ID=_mm512_set1_epi64((long long)KV_IDENT);
    __m512i D[4];
    for(int i=0;i<4;i++){
        __m512i e=_mm512_ternarylogic_epi64(_mm512_loadu_si512(U0+64*i),_mm512_loadu_si512(U1+64*i),
                                            _mm512_loadu_si512(U2+64*i),0x96);
        e=_mm512_xor_si512(e,_mm512_loadu_si512(U3+64*i));
        D[i]=_mm512_gf2p8affineinv_epi64_epi8(e,ID,0);       /* inv(0)=0 */
    }
    D[0]=_mm512_mask_blend_epi8(1,D[0],_mm512_setzero_si512()); /* omega=0 excluded */
    /* 2. chi_d scatter (independent of the scalar corner) */
    u8 dbuf[256] __attribute__((aligned(64)));
    for(int i=0;i<4;i++) _mm512_store_si512(dbuf+64*i,D[i]);
    u64 cd[4]; kv_scatter(dbuf,cd);
    /* 3. power sums */
    __m512i T1=_mm512_setzero_si512(),T3=T1,T5=T1,T7=T1;
    for(int i=0;i<4;i++){
        __m512i d=D[i];
        __m512i d2=_mm512_gf2p8mul_epi8(d,d);
        __m512i d3=_mm512_gf2p8mul_epi8(d2,d);
        __m512i d4=_mm512_gf2p8mul_epi8(d2,d2);
        __m512i d5=_mm512_gf2p8mul_epi8(d4,d);
        __m512i d7=_mm512_gf2p8mul_epi8(d4,d3);
        T1=_mm512_xor_si512(T1,d); T3=_mm512_xor_si512(T3,d3);
        T5=_mm512_xor_si512(T5,d5); T7=_mm512_xor_si512(T7,d7);
    }
    u8 S1=kv_hxor_b(T1),S3=kv_hxor_b(T3),S5=kv_hxor_b(T5),S7=kv_hxor_b(T7);
    /* 4. scalar corner (fp0: identical semantics to kernel_full) */
    u8 S2=KV_SQ[S1],S4=KV_SQ[S2],S6=KV_SQ[S3];
    u8 Q=(u8)(kv_mulz(S1,S6)^kv_mulz(S2,S5)^kv_mulz(S3,S4));  /* P7' (even set) */
    u8 P7=(u8)(S7^Q);                                        /* odd set P7 */
    u8 a0;
    if(__builtin_expect(P7!=0,1)) a0=KV_A73[P7];
    else{
        /* P_11 over the multiset d_1..d_255 (n=255 odd): nS11+S1S10+S2S9+S3S8 */
        __m512i T9=_mm512_setzero_si512(),T11=T9;
        for(int i=0;i<4;i++){
            __m512i d=D[i];
            __m512i d2=_mm512_gf2p8mul_epi8(d,d), d4=_mm512_gf2p8mul_epi8(d2,d2), d8=_mm512_gf2p8mul_epi8(d4,d4);
            __m512i d9=_mm512_gf2p8mul_epi8(d8,d), d3=_mm512_gf2p8mul_epi8(d2,d), d11=_mm512_gf2p8mul_epi8(d8,d3);
            T9=_mm512_xor_si512(T9,d9); T11=_mm512_xor_si512(T11,d11);
        }
        u8 S9=kv_hxor_b(T9),S11=kv_hxor_b(T11),S8=KV_SQ[S4],S10=KV_SQ[S5];
        u8 P11=(u8)(S11^kv_mulz(S1,S10)^kv_mulz(S2,S9)^kv_mulz(S3,S8));
        a0=P11?KV_A116[P11]:(u8)1;
    }
    u8 b0=kv_mulz(a0,S1);
    /* fp1 beta-solve, hoisted (depends only on S1,S3,S5) */
    u8 a1=0, ba=0, bb=0; int s1nz=(S1!=0);
    if(__builtin_expect(s1nz,1)){
        if(__builtin_expect(Q!=0,1)) a1=KV_A73[Q];
        else{
            kvv_stat_p7z++;
            /* set power sum S9 (= multiset sum) -> P11; deeper: scalar chain */
            __m512i T9=_mm512_setzero_si512();
            for(int i=0;i<4;i++){
                __m512i d=D[i];
                __m512i d2=_mm512_gf2p8mul_epi8(d,d), d4=_mm512_gf2p8mul_epi8(d2,d2), d8=_mm512_gf2p8mul_epi8(d4,d4);
                __m512i d9=_mm512_gf2p8mul_epi8(d8,d);
                T9=_mm512_xor_si512(T9,d9);
            }
            u8 S9=kv_hxor_b(T9),S8=KV_SQ[S4],S10=KV_SQ[S5];
            u8 P11=(u8)(kv_mulz(S1,S10)^kv_mulz(S2,S9)^kv_mulz(S3,S8));
            if(P11) a1=KV_A116[P11];
            else{
                /* deep alpha chain (rate ~2^-16): scalar, shared with the ref */
                kvv_stat_deep++;
                chi_t W; W.w[0]=cd[0]^1ULL; W.w[1]=cd[1]; W.w[2]=cd[2]; W.w[3]=cd[3];
                a1=kv2_alpha_deep(&W);
                if(!a1){
                    /* fp0 still needs computing first; then exhaustive tier */
                    __m512i chivD=_mm512_castsi256_si512(_mm256_loadu_si256((const __m256i*)cd));
                    u64 c0d[4]; kv_gather(chivD,GF_inv[a0],b0,c0d);
                    *fp0=kv_fnv4(c0d[0],c0d[1],c0d[2],c0d[3]);
                    *fp1=canon_bruteAGL(&W,NULL);
                    return;
                }
            }
        }
        u8 A=kv_mulz(S3,KV_INVCUBE[S1]);          /* S3/S1^3 (alpha-free) */
        u8 u=CE_UROOT[A];
        u8 aS1=kv_mulz(a1,S1);
        ba=kv_mulz(u,aS1); bb=(u8)(ba^aS1);
    }
    /* 5. gathers + hashes */
    __m512i chiv=_mm512_castsi256_si512(_mm256_loadu_si256((const __m256i*)cd));
    u64 c0[4]; kv_gather(chiv,GF_inv[a0],b0,c0);
    *fp0=kv_fnv4(c0[0],c0[1],c0[2],c0[3]);
    if(__builtin_expect(!s1nz,0)){
        /* rate 2^-8: Designer-1 s1=0 tier (scalar) on the explicit set */
        kvv_stat_s1zero++;
        chi_t W; W.w[0]=cd[0]^1ULL; W.w[1]=cd[1]; W.w[2]=cd[2]; W.w[3]=cd[3];
        *fp1=canon_even(&W,NULL);
        return;
    }
    u64 c1[4]; kv_gather(chiv,GF_inv[a1],0,c1);
    c1[0]^=1ULL;                                  /* + e_0 (a1*0=0 -> bit 0) */
    __m256i chi1=_mm256_loadu_si256((const __m256i*)c1);
    u64 ha=kv_fnv256(kv_shift256v(chi1,ba));
    u64 hb=kv_fnv256(kv_shift256v(chi1,bb));
    *fp1=ha<hb?ha:hb;
}


/* ======================================================================
 * 4-entry batched kernel: same per-entry computation as kernel_v2_vec,
 * phases grouped across 4 independent entries so that their serial
 * latency chains (sum reduction -> scalar corner -> gather -> FNV) overlap
 * in the out-of-order window instead of serializing each entry's tail.
 * Outputs are identical to 4 kernel_v2_vec calls.
 * rows[k] = {U0,U1,U2,U3} for entry k.
 * ==================================================================== */
typedef struct {
    u64 cd[4];                 /* chi_d (raw-set parity histogram) */
    u8 S1,S3,S5,S7;
    u8 a0,b0,a1,ba,bb;
    int rare;                  /* 0 fast path, 1 s1=0, 2 exhaustive */
} kv4_st_t;

static inline void kv4_phase1(const u8*U0,const u8*U1,const u8*U2,const u8*U3,u8*dbuf,kv4_st_t*s){
    const __m512i ID=_mm512_set1_epi64((long long)KV_IDENT);
    __m512i T1=_mm512_setzero_si512(),T3=T1,T5=T1,T7=T1;
    for(int i=0;i<4;i++){
        __m512i e=_mm512_ternarylogic_epi64(_mm512_loadu_si512(U0+64*i),_mm512_loadu_si512(U1+64*i),
                                            _mm512_loadu_si512(U2+64*i),0x96);
        e=_mm512_xor_si512(e,_mm512_loadu_si512(U3+64*i));
        __m512i d=_mm512_gf2p8affineinv_epi64_epi8(e,ID,0);
        if(i==0) d=_mm512_mask_blend_epi8(1,d,_mm512_setzero_si512());
        _mm512_store_si512(dbuf+64*i,d);
        __m512i d2=_mm512_gf2p8mul_epi8(d,d);
        __m512i d3=_mm512_gf2p8mul_epi8(d2,d);
        __m512i d4=_mm512_gf2p8mul_epi8(d2,d2);
        __m512i d5=_mm512_gf2p8mul_epi8(d4,d);
        __m512i d7=_mm512_gf2p8mul_epi8(d4,d3);
        T1=_mm512_xor_si512(T1,d); T3=_mm512_xor_si512(T3,d3);
        T5=_mm512_xor_si512(T5,d5); T7=_mm512_xor_si512(T7,d7);
    }
    s->S1=kv_hxor_b(T1);s->S3=kv_hxor_b(T3);s->S5=kv_hxor_b(T5);s->S7=kv_hxor_b(T7);
}
static inline void kv4_phase_scalar(const u8*dbuf,kv4_st_t*s){
    u8 S1=s->S1,S3=s->S3,S5=s->S5,S7=s->S7;
    u8 S2=KV_SQ[S1],S4=KV_SQ[S2],S6=KV_SQ[S3];
    u8 Q=(u8)(kv_mulz(S1,S6)^kv_mulz(S2,S5)^kv_mulz(S3,S4));
    u8 P7=(u8)(S7^Q);
    u8 a0;
    if(__builtin_expect(P7!=0,1)) a0=KV_A73[P7];
    else{
        __m512i T9=_mm512_setzero_si512(),T11=T9;
        for(int i=0;i<4;i++){
            __m512i d=_mm512_load_si512(dbuf+64*i);
            __m512i d2=_mm512_gf2p8mul_epi8(d,d), d4=_mm512_gf2p8mul_epi8(d2,d2), d8=_mm512_gf2p8mul_epi8(d4,d4);
            __m512i d9=_mm512_gf2p8mul_epi8(d8,d), d3=_mm512_gf2p8mul_epi8(d2,d), d11=_mm512_gf2p8mul_epi8(d8,d3);
            T9=_mm512_xor_si512(T9,d9); T11=_mm512_xor_si512(T11,d11);
        }
        u8 S9=kv_hxor_b(T9),S11=kv_hxor_b(T11),S8=KV_SQ[S4],S10=KV_SQ[S5];
        u8 P11=(u8)(S11^kv_mulz(S1,S10)^kv_mulz(S2,S9)^kv_mulz(S3,S8));
        a0=P11?KV_A116[P11]:(u8)1;
    }
    s->a0=a0; s->b0=kv_mulz(a0,S1);
    s->rare=0;
    if(__builtin_expect(S1==0,0)){ s->rare=1; return; }
    u8 a1;
    if(__builtin_expect(Q!=0,1)) a1=KV_A73[Q];
    else{
        kvv_stat_p7z++;
        __m512i T9=_mm512_setzero_si512();
        for(int i=0;i<4;i++){
            __m512i d=_mm512_load_si512(dbuf+64*i);
            __m512i d2=_mm512_gf2p8mul_epi8(d,d), d4=_mm512_gf2p8mul_epi8(d2,d2), d8=_mm512_gf2p8mul_epi8(d4,d4);
            __m512i d9=_mm512_gf2p8mul_epi8(d8,d);
            T9=_mm512_xor_si512(T9,d9);
        }
        u8 S9=kv_hxor_b(T9),S8=KV_SQ[S4],S10=KV_SQ[S5];
        u8 P11=(u8)(kv_mulz(S1,S10)^kv_mulz(S2,S9)^kv_mulz(S3,S8));
        if(P11) a1=KV_A116[P11];
        else{
            kvv_stat_deep++;
            chi_t W; W.w[0]=s->cd[0]^1ULL; W.w[1]=s->cd[1]; W.w[2]=s->cd[2]; W.w[3]=s->cd[3];
            a1=kv2_alpha_deep(&W);
            if(!a1){ s->rare=2; return; }
        }
    }
    s->a1=a1;
    u8 A=kv_mulz(S3,KV_INVCUBE[S1]);
    u8 u=CE_UROOT[A];
    u8 aS1=kv_mulz(a1,S1);
    s->ba=kv_mulz(u,aS1); s->bb=(u8)(s->ba^aS1);
}
static inline void kernel_v2_vec4(const u8* const rows[4][4], u64 fp0[4], u64 fp1[4]){
    u8 dbuf[4][256] __attribute__((aligned(64)));
    kv4_st_t st[4];
    for(int k=0;k<4;k++) kv4_phase1(rows[k][0],rows[k][1],rows[k][2],rows[k][3],dbuf[k],&st[k]);
    for(int k=0;k<4;k++) kv_scatter(dbuf[k],st[k].cd);
    for(int k=0;k<4;k++) kv4_phase_scalar(dbuf[k],&st[k]);
    for(int k=0;k<4;k++){
        kv4_st_t*s=&st[k];
        __m512i chiv=_mm512_castsi256_si512(_mm256_loadu_si256((const __m256i*)s->cd));
        u64 c0[4]; kv_gather(chiv,GF_inv[s->a0],s->b0,c0);
        fp0[k]=kv_fnv4(c0[0],c0[1],c0[2],c0[3]);
        if(__builtin_expect(s->rare,0)){
            chi_t W; W.w[0]=s->cd[0]^1ULL; W.w[1]=s->cd[1]; W.w[2]=s->cd[2]; W.w[3]=s->cd[3];
            if(s->rare==1){ kvv_stat_s1zero++; fp1[k]=canon_even(&W,NULL); }
            else fp1[k]=canon_bruteAGL(&W,NULL);
            continue;
        }
        u64 c1[4]; kv_gather(chiv,GF_inv[s->a1],0,c1);
        c1[0]^=1ULL;
        __m256i chi1=_mm256_loadu_si256((const __m256i*)c1);
        u64 ha=kv_fnv256(kv_shift256v(chi1,s->ba));
        u64 hb=kv_fnv256(kv_shift256v(chi1,s->bb));
        fp1[k]=ha<hb?ha:hb;
    }
}

#endif
