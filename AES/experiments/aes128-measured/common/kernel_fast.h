// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef KERNEL_FAST_H
#define KERNEL_FAST_H
/* Optimized A1 offline inner kernel: UU-assemble + POWINV + χ-CANON (2-pass).
 * Pass 1: e_ω → d_ω (store) + accumulate S_1357 (u32 XOR).
 * Pass 2: compute (α*,β*), build χ* via amul[d_ω]⊕β, hash.
 * fp1 (EITHER): P_7'=P_7⊕S_7, redo α/β/pass2.
 */
#include "chi256_fast.h"

static u32 POW1357[256]; /* d | d³<<8 | d⁵<<16 | d⁷<<24 */
static void kernel_init(void){
    for(int d=0;d<256;d++){
        if(!d){POW1357[d]=0;continue;}
        int l=GF_log[d];
        POW1357[d]=(u32)d|((u32)GF_exp[(3*l)%255]<<8)|((u32)GF_exp[(5*l)%255]<<16)|((u32)GF_exp[(7*l)%255]<<24);
    }
}

/* α* from (P7,npar,S_k). P7_eff = (npar?S7:0)^S1S6^S2S5^S3S4. */
static inline u8 alpha_from_S(u8 S1,u8 S3,u8 S5,u8 S7,int npar,const u8*dbuf,int n){
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    u8 P7=(npar?S7:0)^gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4);
    if(P7){int l=GF_log[P7];return GF_exp[(255-(l*73)%255)%255];}
    /* fallback P_11 brute over dbuf (rare) */
    u8 P11=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){u8 x=dbuf[i]^dbuf[j];if(x)P11^=GF_exp[(11*GF_log[x])%255];}
    if(P11){int l=GF_log[P11];return GF_exp[(255-(l*116)%255)%255];}
    return 1;
}

/* Full kernel: given UU row-ptrs, compute (fp0,fp1). Returns cyc breakdown via out params. */
static inline void kernel_full(const u8*U0,const u8*U1,const u8*U2,const u8*U3,
                               u64*fp0,u64*fp1){
    u8 d[256]; u32 Sk=0;
    /* pass 1 */
    for(int om=1;om<256;om++){
        u8 e=U0[om]^U1[om]^U2[om]^U3[om];
        u8 dv=GF_inv[e];
        d[om]=dv;
        Sk^=POW1357[dv];
    }
    u8 S1=Sk&0xff,S3=(Sk>>8)&0xff,S5=(Sk>>16)&0xff,S7=(Sk>>24)&0xff;
    /* popcount parity: always ODD for 255-element input (see paper, POW-packing appendix). */
    /* α*,β* for fp0 (npar=1) */
    u8 a0=alpha_from_S(S1,S3,S5,S7,1,d+1,255);
    u8 b0=gmul(a0,S1);
    /* pass 2: χ* */
    chi_t c0; chi_clr(&c0);
    const u8*am0=MUL[a0];
    for(int om=1;om<256;om++) chi_flip(&c0, am0[d[om]]^b0);
    *fp0=chi_hash(&c0);
    /* fp1: χ⊕e_0 → npar=0. S_k unchanged (0 contributes 0). P_7 flips S_7 term.
     * BUT: alpha_from_S's P_11 fallback uses dbuf=d+1 which is the MULTISET,
     * not set-bits of χ⊕e_0. If P_7=0 fallback is hit, P_11 would be wrong.
     * For fp1 we need P_m over set-bits(χ⊕e_0). Compute via χ if fallback needed. */
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    u8 P7_1=gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4); /* npar=0 */
    u8 a1;
    if(P7_1){int l=GF_log[P7_1];a1=GF_exp[(255-(l*73)%255)%255];}
    else{
        /* rare: compute χ⊕e_0 and brute P_11 */
        chi_t ch;chi_clr(&ch);for(int om=1;om<256;om++)chi_flip(&ch,d[om]);chi_flip(&ch,0);
        u8 P11=chi_Pm(&ch,11);
        if(P11){int l=GF_log[P11];a1=GF_exp[(255-(l*116)%255)%255];}else a1=1;
    }
    /* even popcount: β via S_3-solve. */
    u8 b1a,b1b; int nb=0;
    if(S1){
        u8 aS1=gmul(a1,S1);
        u8 rhs=gmul(S3,GF_inv[gmul(S1,gmul(S1,S1))]);
        u8 t=0,cc=rhs;for(int i=0;i<8;i++){t^=cc;cc=gmul(cc,cc);}
        if(t==0){
            u8 g0=GF_htr[rhs]; b1a=gmul(g0,aS1); b1b=gmul(g0^1,aS1); nb=2;
        }
    }
#ifdef KERNEL_FP1_BRUTE
    nb=0;
#endif
    if(nb==0){
        /* S_1=0 or Tr≠0: brute-β. */
        chi_t cs; chi_clr(&cs);
        const u8*am1=MUL[a1];
        for(int om=1;om<256;om++) chi_flip(&cs,am1[d[om]]);
        chi_flip(&cs,0); /* the ⊕e_0, scaled: α·0=0 */
        u64 best=chi_hash(&cs);
        for(int i=1;i<256;i++){int k=__builtin_ctz(i);chi_xorshift_bit(&cs,k);u64 h=chi_hash(&cs);if(h<best)best=h;}
        *fp1=best; return;
    }
    const u8*am1=MUL[a1];
    chi_t c1a,c1b; chi_clr(&c1a);chi_clr(&c1b);
    for(int om=1;om<256;om++){u8 ad=am1[d[om]]; chi_flip(&c1a,ad^b1a); chi_flip(&c1b,ad^b1b);}
    /* ⊕e_0 scaled+shifted: bit (α·0⊕β)=β */
    chi_flip(&c1a,b1a); chi_flip(&c1b,b1b);
    u64 h1a=chi_hash(&c1a),h1b=chi_hash(&c1b);
    *fp1=h1a<h1b?h1a:h1b;
}

#endif
