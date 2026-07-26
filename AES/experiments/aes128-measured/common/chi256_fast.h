// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef CHI256_FAST_H
#define CHI256_FAST_H
/* Fast χ-CANON: α* via P_7(S_k), β* via S_1 (odd popcount) or brute (even).
 * χ-scale(α) via bit-matrix permute of 256-bit word (GF(2)-linear index map).
 */
#include "chi256.h"

/* accumulator: S_1,S_3,S_5,S_7 + χ + popcount-parity */
typedef struct {
    chi_t chi;
    u8 S1,S3,S5,S7;
    int npar; /* popcount parity */
} chiacc_t;
static inline void chiacc_clr(chiacc_t*A){chi_clr(&A->chi);A->S1=A->S3=A->S5=A->S7=0;A->npar=0;}
/* POW357[d] = (d^3<<0)|(d^5<<8)|(d^7<<16) precomputed */
static u32 POW357[256];
static void chifast_init(void){
    for(int d=0;d<256;d++){
        if(d==0){POW357[d]=0;continue;}
        int l=GF_log[d];
        POW357[d]=(u32)GF_exp[(3*l)%255]|((u32)GF_exp[(5*l)%255]<<8)|((u32)GF_exp[(7*l)%255]<<16);
    }
}
static inline void chiacc_add(chiacc_t*A,u8 d){
    u64 m=1ULL<<(d&63); int w=d>>6;
    A->npar ^= (int)((A->chi.w[w]>>(d&63))&1)^1; /* parity flips iff bit was 0 (now 1) or was 1 (now 0): popcount changes by ±1, parity flips always */
    A->chi.w[w]^=m;
    A->S1^=d;
    u32 p=POW357[d]; A->S3^=p&0xff; A->S5^=(p>>8)&0xff; A->S7^=(p>>16)&0xff;
}
/* npar is computed at the end as popcount(chi)&1, not tracked per add:
   popcount parity equals (#adds) mod 2 only if no bit is flipped twice, and
   when d repeats the bit toggles back with popcount unchanged. */

/* 8×8 GF(2) matrix for mul-by-α: col j = bits of α·2^j */
static inline void alpha_matrix(u8 alpha, u8 M[8]){
    u8 aj=alpha;
    for(int j=0;j<8;j++){ M[j]=aj; aj=xt(aj); }
}
/* apply χ'[t]=χ[M^{-1}t] i.e., set-bits d→Md. Compute via: for each set bit d,
   set bit Md. Fast: M is linear, so compute Md via 8 conditional XORs.
   We do it the "slow" way (256 iters) but with cheap inner (no gmul). */
static inline void chi_apply_linmap(const chi_t*in, const u8 M[8], u8 beta, chi_t*out){
    chi_clr(out);
    /* precompute Md for all d via linear combo: Md = XOR_{j: d_j=1} M[j]. Gray-walk. */
    u8 Md=0;
    for(int d=0;d<256;d++){
        if(d){int j=__builtin_ctz(d); Md^=M[j]^((d&(d-1))?0:0);} /* gray? no, direct */
    }
    /* simpler: just loop */
    for(int d=0;d<256;d++){
        if(!chi_get(in,d)) continue;
        u8 t=0; for(int j=0;j<8;j++) if(d&(1<<j)) t^=M[j];
        chi_flip(out, t^beta);
    }
}

static const u64 XMASK[6]={
    0x5555555555555555ULL,0x3333333333333333ULL,0x0f0f0f0f0f0f0f0fULL,
    0x00ff00ff00ff00ffULL,0x0000ffff0000ffffULL,0x00000000ffffffffULL};
/* delta-swap on 4×u64 for β-shift by single bit k: χ'[t]=χ[t⊕2^k]. */
static inline void chi_xorshift_bit(chi_t*c,int k){
    if(k<6){
        u64 m=XMASK[k]; int s=1<<k;
        for(int w=0;w<4;w++){
            u64 lo=c->w[w]&m, hi=c->w[w]&~m;
            c->w[w]=(lo<<s)|(hi>>s);
        }
    } else if(k==6){
        u64 t;
        t=c->w[0];c->w[0]=c->w[1];c->w[1]=t;
        t=c->w[2];c->w[2]=c->w[3];c->w[3]=t;
    } else { /* k==7 */
        u64 t;
        t=c->w[0];c->w[0]=c->w[2];c->w[2]=t;
        t=c->w[1];c->w[1]=c->w[3];c->w[3]=t;
    }
}

/* FAST canon: returns fp. */
static u64 chi_canon_fast(const chiacc_t*A){
    /* popcount parity */
    int npar=(__builtin_popcountll(A->chi.w[0])^__builtin_popcountll(A->chi.w[1])
             ^__builtin_popcountll(A->chi.w[2])^__builtin_popcountll(A->chi.w[3]))&1;
    /* P_7 from S_k */
    u8 S1=A->S1,S3=A->S3,S5=A->S5,S7=A->S7;
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    u8 P7=(npar?S7:0)^gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4);
    u8 alpha;
    if(P7){
        int l=GF_log[P7]; alpha=GF_exp[(255-(l*73)%255)%255];
    } else {
        /* fallback: P_11 brute (rare ~1/256) */
        u8 P11=chi_Pm(&A->chi,11);
        if(P11){int l=GF_log[P11];alpha=GF_exp[(255-(l*116)%255)%255];}
        else alpha=1;
    }
    u8 M[8]; alpha_matrix(alpha,M);
    if(npar){
        /* odd popcount: β*=α·S_1 unique. Permute once. */
        u8 beta=gmul(alpha,S1);
        chi_t cs; chi_apply_linmap(&A->chi,M,beta,&cs);
        return chi_hash(&cs);
    } else {
        /* even popcount: S_1(αD⊕β)=αS_1 (β-free). Use S_3:
           S_3(αD⊕β)=α³S_3⊕α²S_1²β⊕αS_1β². Set=0:
           β²⊕αS_1β = α²S_3/S_1  ⟹  γ²⊕γ=S_3/S_1³, β=γ·αS_1.
           2 roots γ,γ⊕1 → 2 β's. If S_1=0: brute. */
        if(S1){
            u8 aS1=gmul(alpha,S1);
            u8 rhs=gmul(S3,GF_inv[gmul(S1,gmul(S1,S1))]);
            /* Tr(rhs)==0 required for solution; if not, fallback brute */
            u8 tr=rhs; for(int i=0;i<7;i++)tr=gmul(tr,tr)^rhs; /* wrong — compute Tr */
            /* Tr(c)=c+c²+c⁴+...+c^{128} */
            u8 t=0,cc=rhs; for(int i=0;i<8;i++){t^=cc;cc=gmul(cc,cc);}
            if(t==0){
                u8 g0=GF_htr[rhs], g1=g0^1;
                u8 b0=gmul(g0,aS1), b1=gmul(g1,aS1);
                chi_t c0,c1;
                chi_apply_linmap(&A->chi,M,b0,&c0);
                chi_apply_linmap(&A->chi,M,b1,&c1);
                u64 h0=chi_hash(&c0),h1=chi_hash(&c1);
                return h0<h1?h0:h1;
            }
        }
        /* fallback: brute β via Gray walk */
        chi_t cs; chi_apply_linmap(&A->chi,M,0,&cs);
        u64 best=chi_hash(&cs);
        for(int i=1;i<256;i++){
            int k=__builtin_ctz(i);
            chi_xorshift_bit(&cs,k);
            u64 h=chi_hash(&cs);
            if(h<best)best=h;
        }
        return best;
    }
}

/* convenience: build acc from e[] */
static inline void chiacc_from_e(const u8 e[256], chiacc_t*A){
    chiacc_clr(A);
    for(int om=1;om<256;om++){
        u8 d=GF_inv[e[om]];
        chi_flip(&A->chi,d);
        A->S1^=d; u32 p=POW357[d]; A->S3^=p&0xff;A->S5^=(p>>8)&0xff;A->S7^=(p>>16)&0xff;
    }
}

#endif
