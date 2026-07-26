// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef CANON_EVEN_H
#define CANON_EVEN_H
/* ======================================================================
 * canon_even.h -- COMPLETE, cheap AGL(1,256)-canonicalisation of subsets of
 * GF(256), BOTH parities, with NO brute-beta step in the generic even case.
 *
 * See NOTES.md for derivation/proofs.  Summary of the beta-rule:
 *
 *   odd  n            : beta* = s1                      (1 candidate)
 *   even n, s1 != 0   : A = s3/s1^3 (AGL-invariant);
 *                       u = root of u^2+u = A      if Tr(A)=0
 *                         = root of u^2+u = A+tau  if Tr(A)=1  (Tr(tau)=1)
 *                       beta* in { u*s1, (u+1)*s1 }      (2 candidates)
 *   even n, s1 == 0   : L7(c) := s3 c^4 + s5 c^2 + s3^2 c   (GF(2)-linear);
 *                       S_7(W+c) = s7 + L7(c); choose c0 s.t. s7+L7(c0) is
 *                       the canonical rep of the coset s7 + Im(L7);
 *                       beta* in c0 + ker(L7)  (|ker| in {1,2,4}; the
 *                       ker=all case only when s3=s5=0, which also forces
 *                       the alpha-chain to the full-group brute tier).
 *   The canonical representative is the (hash, then lexicographic) minimum
 *   of  C_b = { alpha* d + b : d in V } over the candidate set.
 *
 * alpha* = P_m^{-1/m} with the standard fallback chain (P7,P11,P13 for
 * s1!=0; P23,P29,P31 for even-n & s1=0); deepest fallback = exhaustive
 * min over all 65280 group elements (always correct; rate ~2^-16/2^-24).
 * ====================================================================== */
#include <stdlib.h>
#include "chi256_fast.h"   /* aes_core.h, chi_t, chi_hash, chi_xorshift_bit */

/* ---------------- tables ---------------- */
static u8  CE_TR[256];          /* field trace GF(2^8)->GF(2) */
static u8  CE_TAU;              /* fixed element with Tr(tau)=1 */
static u8  CE_UROOT[256];       /* u with u^2+u = c (+tau if Tr(c)=1) */
static u8  CE_POW[256][256];    /* CE_POW[k][d] = d^k (k>=1); [0][d]=1 */
static int CE_MINV[256];        /* m^{-1} mod 255, or 0 if gcd(m,255)!=1 */
static u8 *CE_TC  = NULL;       /* 2^24 : (s3<<16|s5<<8|s7) -> c0 */
static u32 CE_TKER[65536];      /* (s3<<8|s5) -> dim | b1<<8 | b2<<16 */

static u64 ce_stat_s1zero = 0, ce_stat_fullbrute = 0, ce_stat_chain[8] = {0};

static void canon_even_init(void){
    aes_core_init();
    chifast_init();
    /* trace */
    for(int c=0;c<256;c++){
        u8 t=0,x=(u8)c; for(int i=0;i<8;i++){t^=x;x=gmul(x,x);} CE_TR[c]=t; /* t in {0,1} */
    }
    /* tau: smallest element with Tr=1 */
    CE_TAU=0; for(int c=1;c<256;c++) if(CE_TR[c]==1){CE_TAU=(u8)c;break;}
    /* UROOT */
    for(int c=0;c<256;c++){
        u8 target = CE_TR[c] ? (u8)(c ^ CE_TAU) : (u8)c;
        int found=0;
        for(int u=0;u<256;u++){ if((u8)(gmul((u8)u,(u8)u)^u)==target){CE_UROOT[c]=(u8)u;found=1;break;} }
        if(!found){fprintf(stderr,"UROOT build fail c=%d\n",c);abort();}
    }
    /* powers */
    for(int d=0;d<256;d++){ CE_POW[0][d]=1; for(int k=1;k<256;k++) CE_POW[k][d]=(d==0)?0:GF_exp[(GF_log[d]*k)%255]; }
    /* m^{-1} mod 255 */
    for(int m=0;m<256;m++){ CE_MINV[m]=0; for(int x=1;x<255;x++) if((m*x)%255==1){CE_MINV[m]=x;break;} }
    /* L7 tables */
    CE_TC = (u8*)malloc(1u<<24);
    if(!CE_TC){fprintf(stderr,"malloc CE_TC fail\n");abort();}
    for(int s3=0;s3<256;s3++){
        u8 s3sq=gmul((u8)s3,(u8)s3);
        for(int s5=0;s5<256;s5++){
            /* L7(2^i) */
            u8 img[8];
            for(int i=0;i<8;i++){
                u8 c=(u8)(1<<i), c2=gmul(c,c), c4=gmul(c2,c2);
                img[i]=(u8)(gmul((u8)s3,c4)^gmul((u8)s5,c2)^gmul(s3sq,c));
            }
            /* Gaussian elimination over GF(2): pivots by bit position */
            u8 pivvec[8], pivpre[8]; int haspiv[8]; for(int i=0;i<8;i++){pivvec[i]=0;pivpre[i]=0;haspiv[i]=0;}
            u8 kerb[8]; int kdim=0;
            for(int i=0;i<8;i++){
                u8 v=img[i], pre=(u8)(1<<i);
                for(int b=7;b>=0;b--){
                    if(!((v>>b)&1)) continue;
                    if(haspiv[b]){ v^=pivvec[b]; pre^=pivpre[b]; }
                    else { pivvec[b]=v; pivpre[b]=pre; haspiv[b]=1; v=0; pre=0; break; }
                }
                if(v==0 && pre!=0){ kerb[kdim++]=pre; } /* kernel vector */
            }
            int dimmark = (kdim>2)?8:kdim;
            if(kdim>2 && !(s3==0 && s5==0)){ fprintf(stderr,"kernel dim %d at s3=%d s5=%d\n",kdim,s3,s5); abort(); }
            CE_TKER[(s3<<8)|s5] = (u32)dimmark | ((u32)(kdim>0?kerb[0]:0)<<8) | ((u32)(kdim>1?kerb[1]:0)<<16);
            /* for each s7: reduce, c0 = combo of preimages used */
            for(int s7=0;s7<256;s7++){
                u8 y=(u8)s7, c0=0;
                for(int b=7;b>=0;b--){
                    if(((y>>b)&1) && haspiv[b]){ y^=pivvec[b]; c0^=pivpre[b]; }
                }
                CE_TC[((u32)s3<<16)|((u32)s5<<8)|(u32)s7]=c0;
            }
        }
    }
}

/* ---------------- power sums from chi ---------------- */
/* S[k] = sum over set bits d of d^k, k=1..K (K<=254). */
static void ce_sums(const chi_t*chi,int K,u8*S){
    for(int k=0;k<=K;k++) S[k]=0;
    for(int w=0;w<4;w++){
        u64 x=chi->w[w];
        while(x){
            int j=__builtin_ctzll(x); x&=x-1;
            u8 d=(u8)(64*w+j);
            if(!d) continue;
            int l=GF_log[d];
            for(int k=1;k<=K;k++) S[k]^=GF_exp[(l*k)%255];
        }
    }
}
static inline int ce_popcount(const chi_t*c){
    return __builtin_popcountll(c->w[0])+__builtin_popcountll(c->w[1])
          +__builtin_popcountll(c->w[2])+__builtin_popcountll(c->w[3]);
}
/* generic pairwise power-sum P_m = sum_{i<j}(x_i+x_j)^m from single-index
 * sums (S[0..] with S[0] unused; n = |V|):
 *   P_m = [n odd] S_m + sum_{k: k proper nonzero submask of m, k < m^k} S_k S_{m^k}
 */
static u8 ce_Pm(const u8*S,int n,int m){
    u8 P=(n&1)?S[m]:0;
    for(int k=(m-1)&m;k>0;k=(k-1)&m){
        int k2=m^k;
        if(k<k2) P^=gmul(S[k],S[k2]);
    }
    return P;
}

/* ---------------- set transforms ---------------- */
/* T = alpha * V (elementwise) */
static inline void ce_scale(const chi_t*in,u8 alpha,chi_t*out){
    const u8*am=MUL[alpha];
    out->w[0]=out->w[1]=out->w[2]=out->w[3]=0;
    for(int w=0;w<4;w++){
        u64 x=in->w[w];
        while(x){int j=__builtin_ctzll(x);x&=x-1;u8 t=am[64*w+j];out->w[t>>6]^=1ULL<<(t&63);} /* bijection: XOR==OR */
    }
}
/* in-place translate by b: bits d -> d^b (composition of delta-swaps) */
static inline void ce_xlat_inplace(chi_t*c,u8 b){
    for(int k=0;k<8;k++) if(b&(1<<k)) chi_xorshift_bit(c,k);
}
static inline int ce_lexlt(const chi_t*a,const chi_t*b){
    for(int i=3;i>=0;i--){ if(a->w[i]!=b->w[i]) return a->w[i]<b->w[i]; }
    return 0;
}
static inline int ce_equal(const chi_t*a,const chi_t*b){
    return a->w[0]==b->w[0]&&a->w[1]==b->w[1]&&a->w[2]==b->w[2]&&a->w[3]==b->w[3];
}

/* candidate bookkeeping: keep the (hash,lex)-minimal set */
typedef struct { u64 h; chi_t set; int have; } ce_best_t;
static inline void ce_consider(ce_best_t*B,const chi_t*c){
    u64 h=chi_hash(c);
    if(!B->have || h<B->h || (h==B->h && ce_lexlt(c,&B->set))){ B->h=h; B->set=*c; B->have=1; }
}

/* ---------------- deepest fallback: exhaustive over AGL ---------------- */
static u64 canon_bruteAGL(const chi_t*chi,chi_t*out){
    ce_best_t B; B.have=0;
    for(int a=1;a<256;a++){
        chi_t T; ce_scale(chi,(u8)a,&T);
        ce_consider(&B,&T);
        for(int i=1;i<256;i++){ int k=__builtin_ctz(i); chi_xorshift_bit(&T,k); ce_consider(&B,&T); }
    }
    if(out)*out=B.set;
    return B.h;
}

/* ---------------- the new canonicalisation ---------------- */
static u64 ce_stat_trivial=0;
/* targeted extra odd power sums: S[k] for k in ks[] (plus squares filled by caller) */
static void ce_sums_list(const chi_t*chi,const int*ks,int nk,u8*S){
    for(int i=0;i<nk;i++) S[ks[i]]=0;
    for(int w=0;w<4;w++){
        u64 x=chi->w[w];
        while(x){
            int j=__builtin_ctzll(x); x&=x-1;
            u8 d=(u8)(64*w+j); if(!d) continue;
            int l=GF_log[d];
            for(int i=0;i<nk;i++) S[ks[i]]^=GF_exp[(l*ks[i])%255];
        }
    }
}
static u64 canon_even(const chi_t*chi, chi_t*out){
    int n=ce_popcount(chi);
    /* trivial orbits: n=0 (empty), n=1 (->{0}), n=255 (->GF(256)\{0}),
     * n=256 (full set) are single orbits each; return fixed representatives.
     * (These sizes are also exactly where ALL weight-3 P_m vanish identically.) */
    if(n<=1 || n>=255){
        chi_t C; chi_clr(&C);
        if(n==1) chi_flip(&C,0);
        else if(n==256){ C.w[0]=C.w[1]=C.w[2]=C.w[3]=~0ULL; }
        else if(n==255){ C.w[0]=C.w[1]=C.w[2]=C.w[3]=~0ULL; chi_flip(&C,0); }
        ce_stat_trivial++;
        if(out)*out=C; return chi_hash(&C);
    }
    /* fast sums S1,S3,S5,S7 */
    u8 S1=0,S3=0,S5=0,S7=0;
    for(int w=0;w<4;w++){
        u64 x=chi->w[w];
        while(x){int j=__builtin_ctzll(x);x&=x-1;u8 d=(u8)(64*w+j);S1^=d;u32 p=POW357[d];S3^=(u8)p;S5^=(u8)(p>>8);S7^=(u8)(p>>16);}
    }
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    u8 alpha=0; /* 0 = not yet determined */
    int odd=n&1;

    /* ---- alpha stage ---- */
    if(odd || S1){
        u8 P7=(u8)((odd?S7:0)^gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4));
        if(P7){ alpha=GF_exp[(255-(GF_log[P7]*73)%255)%255]; ce_stat_chain[0]++; }
        else {
            /* need S9,S11,S13 for P11,P13 (popcount*3 lookups, rate ~2^-8) */
            u8 S[16]; S[1]=S1;S[2]=S2;S[3]=S3;S[4]=S4;S[5]=S5;S[6]=S6;S[7]=S7;
            static const int ks[3]={9,11,13}; ce_sums_list(chi,ks,3,S);
            S[8]=gmul(S4,S4); S[10]=gmul(S5,S5); S[12]=gmul(S6,S6);
            static const int chain[2]={11,13};
            for(int t=0;t<2&&!alpha;t++){
                u8 P=ce_Pm(S,n,chain[t]);
                if(P){ int l=GF_log[P]; alpha=GF_exp[(255-(l*CE_MINV[chain[t]])%255)%255]; ce_stat_chain[1+t]++; }
            }
        }
    } else {
        /* even n, S1==0: weight-3 P_m all vanish; use weight-4/5 pair sums.
         * At S1=0, n even (explicit, terms with S1 factors dropped):
         *   P23 = S3 S5^4 + S5 S9^2 + S3^2 S17
         *   P29 = S5 S3^8 + S9 S5^4 + S3^4 S17
         *   P31 = S3 S7^4 + S5 S13^2 + S3^2 S25 + S7 S3^8 + S9 S11^2
         *         + S5^2 S21 + S11 S5^4 + S3^4 S19 + S13 S9^2 + S7^2 S17
         *   (pairs (k,m-k) with S_k, S_{m-k} having no S_1=S_2=S_4=S_8=S_16 factor;
         *    even-index sums replaced by Frobenius squares S_{2j}=S_j^2)
         */
        ce_stat_s1zero++;
        u8 S[32]; static const int ks[2]={9,17}; ce_sums_list(chi,ks,2,S);
        u8 S9=S[9],S17=S[17];
        u8 S5_2=gmul(S5,S5),S5_4=gmul(S5_2,S5_2),S3_4=gmul(S6,S6),S3_8=gmul(S3_4,S3_4),S9_2=gmul(S9,S9);
        u8 P23=(u8)(gmul(S3,S5_4)^gmul(S5,S9_2)^gmul(S6,S17));
        if(P23){ alpha=GF_exp[(255-(GF_log[P23]*CE_MINV[23])%255)%255]; ce_stat_chain[3]++; }
        else {
            u8 P29=(u8)(gmul(S5,S3_8)^gmul(S9,S5_4)^gmul(S3_4,S17));
            if(P29){ alpha=GF_exp[(255-(GF_log[P29]*CE_MINV[29])%255)%255]; ce_stat_chain[4]++; }
            else {
                static const int ks2[5]={11,13,19,21,25}; ce_sums_list(chi,ks2,5,S);
                u8 S11=S[11],S13=S[13],S19=S[19],S21=S[21],S25=S[25];
                u8 S7_2=gmul(S7,S7),S7_4=gmul(S7_2,S7_2),S11_2=gmul(S11,S11),S13_2=gmul(S13,S13);
                u8 P31=(u8)(gmul(S3,S7_4)^gmul(S5,S13_2)^gmul(S6,S25)^gmul(S7,S3_8)^gmul(S9,S11_2)
                            ^gmul(gmul(S5,S5),S21)^gmul(S11,S5_4)^gmul(S3_4,S19)^gmul(S13,S9_2)^gmul(S7_2,S17));
                if(P31){ alpha=GF_exp[(255-(GF_log[P31]*CE_MINV[31])%255)%255]; ce_stat_chain[5]++; }
            }
        }
    }
    if(!alpha){ ce_stat_fullbrute++; return canon_bruteAGL(chi,out); }

    /* ---- scaled frame ---- */
    chi_t T; ce_scale(chi,alpha,&T);
    ce_best_t B; B.have=0;
    if(odd){
        u8 b=gmul(alpha,S1);
        chi_t C=T; ce_xlat_inplace(&C,b); ce_consider(&B,&C);
    } else if(S1){
        u8 s1=gmul(alpha,S1);
        u8 A=gmul(S3,GF_inv[gmul(S1,S2)]);          /* S3/S1^3, alpha-free */
        u8 u=CE_UROOT[A];
        u8 b0=gmul(u,s1);
        chi_t C=T; ce_xlat_inplace(&C,b0); ce_consider(&B,&C);
        ce_xlat_inplace(&C,s1); ce_consider(&B,&C);   /* b0 ^ s1 */
    } else {
        /* s1=0 path: L7 coset */
        u8 a3=gmul(gmul(alpha,alpha),alpha), a5=gmul(gmul(a3,alpha),alpha), a7=gmul(gmul(a5,alpha),alpha);
        u8 s3=gmul(a3,S3), s5=gmul(a5,S5), s7=gmul(a7,S7);
        u8 c0=CE_TC[((u32)s3<<16)|((u32)s5<<8)|(u32)s7];
        u32 kr=CE_TKER[((u32)s3<<8)|(u32)s5];
        int kdim=(int)(kr&0xff); u8 k1=(u8)(kr>>8), k2=(u8)(kr>>16);
        if(kdim>2){ /* cannot happen when alpha came from P23/P29/P31 (s3,s5 not both 0), but stay total */
            ce_stat_fullbrute++; return canon_bruteAGL(chi,out);
        }
        chi_t C=T; ce_xlat_inplace(&C,c0); ce_consider(&B,&C);
        if(kdim>=1){ chi_t C1=C; ce_xlat_inplace(&C1,k1); ce_consider(&B,&C1);
            if(kdim>=2){ chi_t C2=C; ce_xlat_inplace(&C2,k2); ce_consider(&B,&C2);
                         ce_xlat_inplace(&C2,k1); ce_consider(&B,&C2); } }
    }
    if(out)*out=B.set;
    return B.h;
}

/* ---------------- the EXISTING even rule (for timing/cross-check) ----------
 * replicates chi_canon_fast's even branch: Tr=0 -> 2 candidates, else 256-brute.
 * (alpha stage: P7 else P11 chi_Pm brute, as in chi256_fast.h) */
static u64 canon_even_old(const chiacc_t*A, int*was_brute){
    u8 S1=A->S1,S3=A->S3,S5=A->S5;
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    u8 P7=(u8)(gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4));
    u8 alpha;
    if(P7){alpha=GF_exp[(255-(GF_log[P7]*73)%255)%255];}
    else {u8 P11=chi_Pm(&A->chi,11); alpha = P11 ? GF_exp[(255-(GF_log[P11]*116)%255)%255] : 1;}
    chi_t T; ce_scale(&A->chi,alpha,&T);
    if(S1){
        u8 rhs=gmul(S3,GF_inv[gmul(S1,gmul(S1,S1))]);
        if(CE_TR[rhs]==0){
            *was_brute=0;
            u8 aS1=gmul(alpha,S1), g0=GF_htr[rhs];
            u8 b0=gmul(g0,aS1);
            chi_t C=T; ce_xlat_inplace(&C,b0); u64 h0=chi_hash(&C);
            ce_xlat_inplace(&C,aS1); u64 h1=chi_hash(&C);
            return h0<h1?h0:h1;
        }
    }
    *was_brute=1;
    u64 best=chi_hash(&T);
    for(int i=1;i<256;i++){int k=__builtin_ctz(i);chi_xorshift_bit(&T,k);u64 h=chi_hash(&T);if(h<best)best=h;}
    return best;
}

/* fast-path-only version of the NEW even rule on an accumulator (timing) */
static u64 canon_even_new_acc(const chiacc_t*A, chi_t*outset){
    u8 S1=A->S1,S3=A->S3,S5=A->S5;
    if(!S1){ return canon_even(&A->chi,outset); } /* rare tier */
    u8 S2=gmul(S1,S1),S4=gmul(S2,S2),S6=gmul(S3,S3);
    u8 P7=(u8)(gmul(S1,S6)^gmul(S2,S5)^gmul(S3,S4));
    if(!P7){ return canon_even(&A->chi,outset); }    /* rare tier (P11 via sums) */
    u8 alpha=GF_exp[(255-(GF_log[P7]*73)%255)%255];
    chi_t T; ce_scale(&A->chi,alpha,&T);
    u8 s1=gmul(alpha,S1);
    u8 Aq=gmul(S3,GF_inv[gmul(S1,S2)]);
    u8 b0=gmul(CE_UROOT[Aq],s1);
    ce_xlat_inplace(&T,b0); u64 h0=chi_hash(&T);
    chi_t T2=T; ce_xlat_inplace(&T2,s1); u64 h1=chi_hash(&T2);
    if(h0<h1||(h0==h1&&ce_lexlt(&T,&T2))){ if(outset)*outset=T; return h0; }
    if(outset)*outset=T2; return h1;
}

#endif
