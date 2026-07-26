// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef KERNEL_V2_SCALAR_H
#define KERNEL_V2_SCALAR_H
/* kernel_full_v2: the A1 offline table kernel with the even-parity
 * brute-beta REMOVED.
 *
 * Identical to kernel_full (kernel_fast.h) except the fp1 beta-stage:
 *   fp1 set = W1 = (chi XOR e_0) = parity set of {d_1..d_255} with bit 0
 *   flipped; |W1| is always EVEN.  beta-rule (Designer 1, canon_even.h):
 *     S1 != 0 : A = S3/S1^3 (AGL-invariant), u = CE_UROOT[A]
 *               (u^2+u = A if Tr(A)=0 else A+tau), candidates
 *               { u*alpha*S1, (u+1)*alpha*S1 } -> min of 2 hashes.
 *               (old code: only the Tr(A)=0 half took this path, with the
 *                SAME candidate pair {g0*aS1,(g0+1)*aS1}; Tr(A)=1 half and
 *                the S1=0 sets went through a 256-translate brute min.)
 *     S1 == 0 : (rate 2^-8) Designer 1's s1=0 tier verbatim:
 *               alpha from P23->P29->P31, beta in c0+ker(L7); called via
 *               canon_even() on the explicit set.
 *   alpha-chain for fp1 (S1!=0): P7' (even-n P7, no S7 term) -> P11 -> P13
 *   via power-sum closed forms (Designer 1 tiers), deepest tier exhaustive
 *   AGL brute (canon_bruteAGL).  The OLD kernel's chain was
 *   P7' -> P11(O(n^2) on the set) -> alpha=1 (non-invariant degenerate
 *   fallback); the closed-form P11 equals the O(n^2) one, so fp1 is
 *   bit-identical to the old kernel on every entry the old kernel did NOT
 *   brute-force AND where its alpha came from P7' or P11.
 * fp0 is byte-for-byte the old code path (npar=1, odd set, unique beta*).
 */
#include "kernel_fast.h"
#include "canon_even.h"

/* d^9 | d^11<<8 | d^13<<16 | d^17<<24 */
static u32 POW9_11_13[256];
static void kernel_v2_init(void){
    kernel_init();
    canon_even_init();
    for(int d=0;d<256;d++){
        if(!d){POW9_11_13[d]=0;continue;}
        int l=GF_log[d];
        POW9_11_13[d]=(u32)GF_exp[(9*l)%255]|((u32)GF_exp[(11*l)%255]<<8)|((u32)GF_exp[(13*l)%255]<<16)
                     |((u32)GF_exp[(17*l)%255]<<24);
    }
}

/* optional path statistics */
static u64 kv2_stat_s1zero=0, kv2_stat_p7z=0, kv2_stat_p11z=0, kv2_stat_deepx=0, kv2_stat_trA1=0;

/* Deep alpha chain for the (S1!=0, even n, P7'=P11=0) tier, rate ~2^-16.
 * Measured on this entry distribution: P7'=P11=0 always forces P13=0, so the
 * real next rung is P19; P19 itself vanishes at ~2^-3 given P7'=P11=0.
 * Chain of the remaining m coprime to 255 with weight>=3 (P_m via the
 * closed power-sum form, set sums recomputed from the explicit set);
 * returns alpha* or 0 (-> exhaustive AGL minimum, Designer 1's last tier). */
static const int KV2_DEEP_M[]={13,19,23,29,31,37,41,43,47,49,53,59,61,0};
static u8 kv2_alpha_deep(const chi_t*W){
    u8 S[64]; ce_sums(W,62,S);
    int n=ce_popcount(W);
    for(int k=0;KV2_DEEP_M[k];k++){
        int m=KV2_DEEP_M[k];
        u8 P=ce_Pm(S,n,m);
        if(P){ int l=GF_log[P]; return GF_exp[(255-(l*CE_MINV[m])%255)%255]; }
    }
    return 0;
}

static inline void kernel_full_v2(const u8*U0,const u8*U1,const u8*U2,const u8*U3,
                                  u64*fp0,u64*fp1){
    u8 d[256]; u32 Sk=0;
    /* pass 1 (identical to kernel_full) */
    for(int om=1;om<256;om++){
        u8 e=U0[om]^U1[om]^U2[om]^U3[om];
        u8 dv=GF_inv[e];
        d[om]=dv;
        Sk^=POW1357[dv];
    }
    u8 S1=Sk&0xff,S3=(Sk>>8)&0xff,S5=(Sk>>16)&0xff,S7=(Sk>>24)&0xff;
    /* fp0 (identical to kernel_full: npar=1) */
    u8 a0=alpha_from_S(S1,S3,S5,S7,1,d+1,255);
    u8 b0=gmul(a0,S1);
    chi_t c0; chi_clr(&c0);
    const u8*am0=MUL[a0];
    for(int om=1;om<256;om++) chi_flip(&c0, am0[d[om]]^b0);
    *fp0=chi_hash(&c0);

    /* fp1: even set W1 = parity{d} XOR e_0 */
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    if(__builtin_expect(S1==0,0)){
        /* rate 2^-8: Designer 1's s1=0 tier (P23/P29/P31 + L7 coset),
         * trivial sizes, and exhaustive deepest tier -- all inside canon_even. */
        kv2_stat_s1zero++;
        chi_t ch; chi_clr(&ch);
        for(int om=1;om<256;om++) chi_flip(&ch,d[om]);
        chi_flip(&ch,0);
        *fp1=canon_even(&ch,NULL);
        return;
    }
    u8 P7_1=gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4); /* even n: no S7 term */
    u8 a1;
    if(__builtin_expect(P7_1!=0,1)){int l=GF_log[P7_1];a1=GF_exp[(255-(l*73)%255)%255];}
    else{
        kv2_stat_p7z++;
        /* S9 over the set (multiset sums: duplicates cancel, 0 adds 0) */
        u32 Tk=0; for(int om=1;om<256;om++) Tk^=POW9_11_13[d[om]];
        u8 S9=(u8)Tk;
        u8 S8=gmul(S4,S4),S10=gmul(S5,S5);
        /* n even: P11 = S1S10 + S2S9 + S3S8 */
        u8 P11=gmul(S1,S10)^gmul(S2,S9)^gmul(S3,S8);
        if(P11){int l=GF_log[P11];a1=GF_exp[(255-(l*116)%255)%255];}
        else{
            /* deep tier (~2^-16): continue the weight>=3 chain on the set */
            kv2_stat_p11z++;
            chi_t ch; chi_clr(&ch);
            for(int om=1;om<256;om++) chi_flip(&ch,d[om]);
            chi_flip(&ch,0);
            a1=kv2_alpha_deep(&ch);
            if(!a1){
                /* deepest tier: exhaustive AGL minimum (Designer 1) */
                kv2_stat_deepx++;
                *fp1=canon_bruteAGL(&ch,NULL);
                return;
            }
        }
    }
    /* two-candidate beta (tau-trick): always, no Tr test, no brute */
    u8 A=gmul(S3,GF_inv[gmul(S1,S2)]);     /* S3/S1^3 (alpha-free) */
    if(CE_TR[A]) kv2_stat_trA1++;          /* old code would have bruted these */
    u8 u=CE_UROOT[A];
    u8 aS1=gmul(a1,S1), ba=gmul(u,aS1), bb=(u8)(ba^aS1);
    const u8*am1=MUL[a1];
    chi_t ca; chi_clr(&ca);
    for(int om=1;om<256;om++) chi_flip(&ca,(u8)(am1[d[om]]^ba));
    chi_flip(&ca,ba);                      /* the (+) e_0 element: a1*0 + ba */
    chi_t cb=ca;
    for(int k=0;k<8;k++) if(aS1&(1<<k)) chi_xorshift_bit(&cb,k); /* cb = ca translated by aS1 */
    u64 ha=chi_hash(&ca),hb=chi_hash(&cb);
    *fp1=ha<hb?ha:hb;
}

#endif
