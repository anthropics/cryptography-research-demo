// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* saes.h — Small-scale AES SR*(4,4,4): 4×4 state of GF(2^4) nibbles.
 * GF(16) = GF(2)[x]/(x^4+x+1).  S-box = A·inv(x) ⊕ c (affine∘inverse).
 * MC = circ(2,3,1,1), iMC = circ(E,B,D,9).  Column-major nibble layout.
 * 7 rounds: whitening + 6 full rounds + 1 final (no MC).
 * Key schedule: AES-128 analogue on 4 nibble-columns.
 */
#ifndef SAES_H
#define SAES_H
#include <stdint.h>
#include <string.h>

/* ---------- GF(16) ---------- */
static uint8_t sgf_mul_tab[16][16];
static uint8_t sgf_inv_tab[16];
static uint8_t sgf_log[16], sgf_alog[16];

static inline uint8_t sgf_mul(uint8_t a,uint8_t b){ return sgf_mul_tab[a&15][b&15]; }
static inline uint8_t sgf_inv(uint8_t a){ return sgf_inv_tab[a&15]; }
static inline uint8_t sgf_pow(uint8_t a,int e){
    if(!a) return e==0?1:0;
    return sgf_alog[((sgf_log[a]*(e%15))%15+15)%15];
}

/* ---------- S-box: A·inv(x)⊕c over GF(2)^4 ---------- */
/* Cid-Murphy-Robshaw small-scale AES S-box (one common choice):
 * Use the 4-bit affine map that makes S-box have no fixed points and is
 * differentially 4-uniform. We pick A,c so SB(0)=c≠0 and SB is bijective.
 * A (rows, LSB first): [1 0 1 1; 1 1 0 1; 1 1 1 0; 0 1 1 1], c=0x6.
 * (Any invertible A works for the bridge; this one matches CMR'05.)
 */
static uint8_t sSB[16], siSB[16];
static uint8_t sA_mat[4] = {0x0D,0x0B,0x07,0x0E}; /* row i as 4-bit mask (LSB=bit0) */
static uint8_t sA_c = 0x6;
static uint8_t sMinv[16]; /* M^{-1} lookup (GF(2)-linear part of S-box) */

static inline uint8_t saff(uint8_t x){
    uint8_t y=0;
    for(int i=0;i<4;i++) y |= (__builtin_parity(x & sA_mat[i])&1)<<i;
    return y;
}

/* ---------- round ops (16-nibble state in uint8_t[16], low nibble used) ---------- */
static inline void ssb(uint8_t*s){for(int i=0;i<16;i++)s[i]=sSB[s[i]];}
static inline void sisb(uint8_t*s){for(int i=0;i<16;i++)s[i]=siSB[s[i]];}
static inline void ssr(uint8_t*s){
    uint8_t t[16]; for(int c=0;c<4;c++)for(int r=0;r<4;r++)t[4*c+r]=s[4*((c+r)&3)+r];
    memcpy(s,t,16);
}
static inline void sisr(uint8_t*s){
    uint8_t t[16]; for(int c=0;c<4;c++)for(int r=0;r<4;r++)t[4*c+r]=s[4*((c+4-r)&3)+r];
    memcpy(s,t,16);
}
static inline void smc_col(uint8_t*c){
    uint8_t a0=c[0],a1=c[1],a2=c[2],a3=c[3];
    c[0]=sgf_mul(2,a0)^sgf_mul(3,a1)^a2^a3;
    c[1]=a0^sgf_mul(2,a1)^sgf_mul(3,a2)^a3;
    c[2]=a0^a1^sgf_mul(2,a2)^sgf_mul(3,a3);
    c[3]=sgf_mul(3,a0)^a1^a2^sgf_mul(2,a3);
}
static inline void simc_col(uint8_t*c){
    uint8_t a0=c[0],a1=c[1],a2=c[2],a3=c[3];
    c[0]=sgf_mul(0xE,a0)^sgf_mul(0xB,a1)^sgf_mul(0xD,a2)^sgf_mul(0x9,a3);
    c[1]=sgf_mul(0x9,a0)^sgf_mul(0xE,a1)^sgf_mul(0xB,a2)^sgf_mul(0xD,a3);
    c[2]=sgf_mul(0xD,a0)^sgf_mul(0x9,a1)^sgf_mul(0xE,a2)^sgf_mul(0xB,a3);
    c[3]=sgf_mul(0xB,a0)^sgf_mul(0xD,a1)^sgf_mul(0x9,a2)^sgf_mul(0xE,a3);
}
static inline void smc(uint8_t*s){for(int c=0;c<4;c++)smc_col(s+4*c);}
static inline void simc(uint8_t*s){for(int c=0;c<4;c++)simc_col(s+4*c);}
static inline void sark(uint8_t*s,const uint8_t*k){for(int i=0;i<16;i++)s[i]^=k[i];}

/* ---------- key schedule (AES-128 analogue) ---------- */
static const uint8_t sRcon[11]={0,1,2,4,8,3,6,12,11,5,10}; /* powers of 2 in GF(16) */
static void skey_expand(const uint8_t key[16], uint8_t rk[11][16]){
    memcpy(rk[0],key,16);
    for(int r=1;r<=10;r++){
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        uint8_t tt=t[0]; t[0]=sSB[t[1]];t[1]=sSB[t[2]];t[2]=sSB[t[3]];t[3]=sSB[tt];
        t[0]^=sRcon[r];
        for(int i=0;i<4;i++) rk[r][i]=rk[r-1][i]^t[i];
        for(int c=1;c<4;c++)for(int i=0;i<4;i++) rk[r][4*c+i]=rk[r-1][4*c+i]^rk[r][4*(c-1)+i];
    }
}
/* derive full schedule from rk[r0] */
static void skey_from(int r0, const uint8_t kr0[16], uint8_t rk[11][16]){
    memcpy(rk[r0],kr0,16);
    for(int r=r0;r>=1;r--){
        for(int c=3;c>=1;c--)for(int i=0;i<4;i++) rk[r-1][4*c+i]=rk[r][4*c+i]^rk[r][4*(c-1)+i];
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        uint8_t tt=t[0]; t[0]=sSB[t[1]];t[1]=sSB[t[2]];t[2]=sSB[t[3]];t[3]=sSB[tt];
        t[0]^=sRcon[r];
        for(int i=0;i<4;i++) rk[r-1][i]=rk[r][i]^t[i];
    }
    for(int r=r0+1;r<=10;r++){
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        uint8_t tt=t[0]; t[0]=sSB[t[1]];t[1]=sSB[t[2]];t[2]=sSB[t[3]];t[3]=sSB[tt];
        t[0]^=sRcon[r];
        for(int i=0;i<4;i++) rk[r][i]=rk[r-1][i]^t[i];
        for(int c=1;c<4;c++)for(int i=0;i<4;i++) rk[r][4*c+i]=rk[r-1][4*c+i]^rk[r][4*(c-1)+i];
    }
}
/* invert key schedule: rk[7] → rk[0] */
static void skey_invert_from7(const uint8_t rk7[16], uint8_t rk[11][16]){
    memcpy(rk[7],rk7,16);
    for(int r=7;r>=1;r--){
        for(int c=3;c>=1;c--)for(int i=0;i<4;i++) rk[r-1][4*c+i]=rk[r][4*c+i]^rk[r][4*(c-1)+i];
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        uint8_t tt=t[0]; t[0]=sSB[t[1]];t[1]=sSB[t[2]];t[2]=sSB[t[3]];t[3]=sSB[tt];
        t[0]^=sRcon[r];
        for(int i=0;i<4;i++) rk[r-1][i]=rk[r][i]^t[i];
    }
}

/* ---------- 7-round encrypt with trace ---------- */
/* trace[i] = x_i (i=0..6), trace[7]=C */
static void senc7_trace(const uint8_t rk[11][16], const uint8_t P[16], uint8_t tr[8][16]){
    uint8_t s[16]; memcpy(s,P,16); sark(s,rk[0]);
    memcpy(tr[0],s,16);
    for(int r=0;r<6;r++){ ssb(s); ssr(s); smc(s); sark(s,rk[r+1]); memcpy(tr[r+1],s,16); }
    ssb(s); ssr(s); sark(s,rk[7]); memcpy(tr[7],s,16);
}
static void senc7(const uint8_t rk[11][16], const uint8_t P[16], uint8_t C[16]){
    uint8_t tr[8][16]; senc7_trace(rk,P,tr); memcpy(C,tr[7],16);
}

/* ---------- init ---------- */
static void saes_init(void){
    /* GF(16) mul */
    for(int a=0;a<16;a++)for(int b=0;b<16;b++){
        int aa=a,bb=b,p=0;
        for(int i=0;i<4;i++){ if(bb&1)p^=aa; bb>>=1; int hi=aa&8; aa=(aa<<1)&15; if(hi)aa^=3; }
        sgf_mul_tab[a][b]=p;
    }
    /* log/alog (generator 2? check: 2^1=2,2^2=4,2^3=8,2^4=3,2^5=6,2^6=12,2^7=11,
       2^8=5,2^9=10,2^10=7,2^11=14,2^12=15,2^13=13,2^14=9,2^15=1 — yes, 2 is primitive) */
    {uint8_t x=1; for(int i=0;i<15;i++){sgf_alog[i]=x;sgf_log[x]=i;x=sgf_mul(x,2);} }
    /* inv */
    sgf_inv_tab[0]=0; for(int a=1;a<16;a++)for(int b=1;b<16;b++) if(sgf_mul(a,b)==1){sgf_inv_tab[a]=b;break;}
    /* S-box */
    for(int x=0;x<16;x++){ sSB[x]=saff(sgf_inv_tab[x])^sA_c; }
    for(int x=0;x<16;x++) siSB[sSB[x]]=x;
    /* Minv: SB(x)=A(inv(x))^c ⟹ A^{-1}(w)=inv(iSB[w^c]) */
    for(int w=0;w<16;w++) sMinv[w]=sgf_inv_tab[siSB[w^sA_c]];
}

#endif
