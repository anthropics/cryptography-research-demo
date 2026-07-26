// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef CHI256_H
#define CHI256_H
/* χ-CANON at q=256: AGL(1,256)-invariant fingerprint of χ∈GF(2)^{256}.
 *   χ_D[d] = parity of #{ω: d_ω=d}
 *   AGL action: π_{α,β}(χ)[t]=χ[α^{-1}(t⊕β)], equiv: set-bits map d↦αd⊕β.
 *   canon(χ): fix α* via P_7 (scale-cov, trans-inv), then min over β.
 *     fp(χ) = min_{β} FNV64(π_{α*,β}(χ))
 *   This is AGL-invariant (proof: β-orbit unchanged by pre-composition).
 * EITHER: fp0=canon(χ), fp1=canon(χ⊕e_0). Bridge matches fp0==fp0 OR fp1==fp1.
 */
#include "aes_core.h"

typedef struct { u64 w[4]; } chi_t;
static inline void chi_clr(chi_t*c){c->w[0]=c->w[1]=c->w[2]=c->w[3]=0;}
static inline void chi_flip(chi_t*c,u8 d){c->w[d>>6]^=1ULL<<(d&63);}
static inline int  chi_get(const chi_t*c,u8 d){return (int)((c->w[d>>6]>>(d&63))&1);}
static inline u64  chi_hash(const chi_t*c){
    u64 h=0xcbf29ce484222325ULL;
    for(int i=0;i<4;i++){h^=c->w[i];h*=0x100000001b3ULL;}
    return h;
}

/* P_m over set-bits of χ (trans-inv, scale-cov deg m). Brute O(n²), n=popcount. */
static u8 chi_Pm(const chi_t*chi, int m){
    u8 idx[256]; int n=0;
    for(int d=0;d<256;d++) if(chi_get(chi,d)) idx[n++]=(u8)d;
    u8 P=0;
    for(int i=0;i<n;i++){
        int li=GF_log[idx[i]]; /* unused if idx[i]==0, but idx[i]^idx[j]≠0 below */
        (void)li;
        for(int j=i+1;j<n;j++){
            u8 x=idx[i]^idx[j];
            P^=GF_exp[(GF_log[x]*m)%255];
        }
    }
    return P;
}
/* Fast P_7 via power-sums S_k of set-bits (O(n)). */
static u8 chi_P7_fast(const chi_t*chi, u8*S1_out){
    u8 S1=0,S3=0,S5=0,S7=0; int n=0;
    for(int d=0;d<256;d++) if(chi_get(chi,d)){
        n++;
        if(d){int l=GF_log[d];S1^=d;S3^=GF_exp[(3*l)%255];S5^=GF_exp[(5*l)%255];S7^=GF_exp[(7*l)%255];}
    }
    if(S1_out)*S1_out=S1;
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    /* P_7 = Σ_{i<j}(d_i⊕d_j)^7. Expand (Lucas, all C(7,k) odd):
       = (n-1 mod2)S_7 ⊕ (S_1S_6⊕S_7)⊕(S_2S_5⊕S_7)⊕(S_3S_4⊕S_7)
       pairs (0,7): Σ_{i<j}(d_i^7+d_j^7)=(n-1)S_7.
       pairs (k,7-k) k=1,2,3: Σ_{i<j}(d_i^k d_j^{7-k}+d_i^{7-k}d_j^k)=S_k S_{7-k}⊕S_7.
       total: (n-1)S_7 ⊕ S_1S_6⊕S_2S_5⊕S_3S_4 ⊕ 3·S_7
            = (n-1)S_7 ⊕ S_7 ⊕ S_1S_6⊕S_2S_5⊕S_3S_4   (3 mod2=1)
            = n·S_7 ⊕ S_1S_6⊕S_2S_5⊕S_3S_4. */
    u8 P7 = ((n&1)?S7:0) ^ gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4);
    return P7;
}

/* α* from P_m chain (coprime to 255=3·5·17). */
static u8 chi_alpha(const chi_t*chi){
    static const int ms[]={7,11,13,19,23,29,0};
    static const int invm[]={73,116,157,94,111,44,0}; /* m^{-1} mod 255 */
    for(int k=0;ms[k];k++){
        u8 Pm = (ms[k]==7)? chi_P7_fast(chi,NULL) : chi_Pm(chi,ms[k]);
        if(Pm){
            int l=GF_log[Pm];
            int la=(255 - (int)((long)l*invm[k]%255))%255;
            return GF_exp[la];
        }
    }
    return 1; /* degenerate */
}

/* Apply π_{α,0}: set-bits d→αd. */
static void chi_scale(const chi_t*in, u8 alpha, chi_t*out){
    chi_clr(out);
    for(int d=0;d<256;d++) if(chi_get(in,d)) chi_flip(out, gmul(alpha,(u8)d));
}
/* Apply π_{1,β}: set-bits d→d⊕β. Equiv: XOR-index-permute. */
static void chi_xlat(const chi_t*in, u8 beta, chi_t*out){
    chi_clr(out);
    for(int d=0;d<256;d++) if(chi_get(in,d)) chi_flip(out, (u8)(d^beta));
}

/* canonical fp (brute-β). */
static u64 chi_canon(const chi_t*chi){
    u8 alpha=chi_alpha(chi);
    chi_t cs; chi_scale(chi,alpha,&cs);
    u64 best=~0ULL;
    for(int beta=0;beta<256;beta++){
        chi_t ct; chi_xlat(&cs,(u8)beta,&ct);
        u64 h=chi_hash(&ct);
        if(h<best) best=h;
    }
    return best;
}

/* from e[1..255] → d=1/e → χ → (fp0,fp1). */
static void chi_fp_from_e(const u8 e[256], u64*fp0, u64*fp1){
    chi_t chi; chi_clr(&chi);
    for(int om=1;om<256;om++) chi_flip(&chi, GF_inv[e[om]]);
    *fp0=chi_canon(&chi);
    chi_flip(&chi,0);
    *fp1=chi_canon(&chi);
}
/* from d[1..255] directly */
static void chi_fp_from_d(const u8 d[256], u64*fp0, u64*fp1){
    chi_t chi; chi_clr(&chi);
    for(int om=1;om<256;om++) chi_flip(&chi, d[om]);
    *fp0=chi_canon(&chi);
    chi_flip(&chi,0);
    *fp1=chi_canon(&chi);
}

#endif
