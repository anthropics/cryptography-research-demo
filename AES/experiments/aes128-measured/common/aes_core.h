// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef AES_CORE_H
#define AES_CORE_H
/* AES-128 primitives for A1 e2e attack.
 * State indexing: col-major, p = 4*col + row (FIPS-197 column-major).
 *   col0={0,1,2,3}  diag0={0,5,10,15}  adiag={0,13,10,7}
 * Round r: x_r --SB--> y_r --SR--> z_r --MC--> w_r --ARK(k_r)--> x_{r+1}
 * Last round (r=NROUNDS): no MC.
 */
#include <stdint.h>
#include <string.h>
#include <stdio.h>

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

/* ---------- GF(256) ---------- */
static u8 GF_log[256], GF_exp[510], GF_inv[256], GF_sqrt[256], GF_htr[256];
static u8 SBOX[256], iSBOX[256];
static u8 MUL[256][256];
/* Affine layer of SBOX: SBOX(x)=AFF(inv(x))=L(inv(x))⊕0x63.
 * Linv = inverse of LINEAR part L (for bridge: Δy=L(Δinv)). */
static u8 AFF[256], iAFF[256], Lfwd[256], Linv[256];

static const u8 MCc[4][4]  = {{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
static const u8 iMCc[4][4] = {{14,11,13,9},{9,14,11,13},{13,9,14,11},{11,13,9,14}};

static inline u8 xt(u8 a){ return (u8)((a<<1)^((a&0x80)?0x1b:0)); }
static inline u8 gf_mul_slow(u8 a,u8 b){u8 r=0;while(b){if(b&1)r^=a;a=xt(a);b>>=1;}return r;}
static inline u8 gmul(u8 a,u8 b){ return MUL[a][b]; }
static inline u8 gpow(u8 a,int e){
    if(a==0) return e==0?1:0;
    int l=GF_log[a]; return GF_exp[((long)l*e)%255 + ((((long)l*e)%255)<0?255:0)];
}

static void aes_core_init(void){
    /* log/exp with generator 3 */
    u8 x=1; for(int i=0;i<255;i++){GF_exp[i]=x;GF_log[x]=i;x=gf_mul_slow(x,3);}
    for(int i=255;i<510;i++) GF_exp[i]=GF_exp[i-255];
    GF_log[0]=0; /* sentinel */
    GF_inv[0]=0; for(int i=1;i<256;i++) GF_inv[i]=GF_exp[255-GF_log[i]];
    for(int i=0;i<256;i++) GF_sqrt[gf_mul_slow(i,i)]=i;
    /* solve x²+x=c: brute (GF(2^8) even-n, no simple half-trace). GF_htr[c]=one root (or 0 if none). */
    for(int c=0;c<256;c++){GF_htr[c]=0;for(int x=0;x<256;x++)if((gf_mul_slow(x,x)^x)==c){GF_htr[c]=x;break;}}
    for(int a=0;a<256;a++) for(int b=0;b<256;b++) MUL[a][b]=gf_mul_slow(a,b);
    /* SBOX = affine(inv(x)) */
    for(int x=0;x<256;x++){
        u8 s=GF_inv[x], r=s;
        for(int i=0;i<4;i++){ s=(u8)((s<<1)|(s>>7)); r^=s; }
        r^=0x63;
        SBOX[x]=r; iSBOX[r]=x;
    }
    /* AFF(t)=L(t)⊕0x63;  L(t)=t⊕rol(t,1)⊕rol(t,2)⊕rol(t,3)⊕rol(t,4) */
    for(int t=0;t<256;t++){
        u8 s=t,r=t; for(int i=0;i<4;i++){s=(u8)((s<<1)|(s>>7));r^=s;}
        Lfwd[t]=r; Linv[r]=t;
        AFF[t]=r^0x63; iAFF[r^0x63]=t;
    }
    /* sanity */
    if(SBOX[0]!=0x63||SBOX[1]!=0x7c||iSBOX[0x63]!=0){fprintf(stderr,"SBOX FAIL\n");}
}

/* ---------- state ops ---------- */
static inline void st_cpy(u8*d,const u8*s){memcpy(d,s,16);}
static inline void st_xor(u8*d,const u8*a,const u8*b){for(int i=0;i<16;i++)d[i]=a[i]^b[i];}
static inline void sub_bytes(u8*s){for(int i=0;i<16;i++)s[i]=SBOX[s[i]];}
static inline void isub_bytes(u8*s){for(int i=0;i<16;i++)s[i]=iSBOX[s[i]];}
static inline void shift_rows(u8*s){
    u8 t[16]; for(int c=0;c<4;c++)for(int r=0;r<4;r++) t[4*c+r]=s[4*((c+r)&3)+r];
    memcpy(s,t,16);
}
static inline void ishift_rows(u8*s){
    u8 t[16]; for(int c=0;c<4;c++)for(int r=0;r<4;r++) t[4*c+r]=s[4*((c-r+4)&3)+r];
    memcpy(s,t,16);
}
static inline void mix_cols(u8*s){
    for(int c=0;c<4;c++){u8 a[4]={s[4*c],s[4*c+1],s[4*c+2],s[4*c+3]};
        for(int r=0;r<4;r++) s[4*c+r]=gmul(MCc[r][0],a[0])^gmul(MCc[r][1],a[1])^gmul(MCc[r][2],a[2])^gmul(MCc[r][3],a[3]);}
}
static inline void imix_cols(u8*s){
    for(int c=0;c<4;c++){u8 a[4]={s[4*c],s[4*c+1],s[4*c+2],s[4*c+3]};
        for(int r=0;r<4;r++) s[4*c+r]=gmul(iMCc[r][0],a[0])^gmul(iMCc[r][1],a[1])^gmul(iMCc[r][2],a[2])^gmul(iMCc[r][3],a[3]);}
}

/* ---------- key schedule (bytes) ---------- */
typedef struct { u8 rk[11][16]; int nr; } aesctx;
static void aes_set_key(aesctx*ctx,const u8 key[16],int nr){
    static const u8 Rcon[11]={0,1,2,4,8,16,32,64,128,0x1b,0x36};
    ctx->nr=nr; memcpy(ctx->rk[0],key,16);
    for(int r=1;r<=10;r++){
        u8*p=ctx->rk[r-1];u8*q=ctx->rk[r];
        /* words are columns (4 bytes each) */
        u8 t[4]={SBOX[p[13]],SBOX[p[14]],SBOX[p[15]],SBOX[p[12]]};
        t[0]^=Rcon[r];
        for(int i=0;i<4;i++) q[i]=p[i]^t[i];
        for(int i=4;i<16;i++) q[i]=p[i]^q[i-4];
    }
}

/* ---------- encrypt with full trace ---------- */
typedef struct { u8 x[12][16],y[12][16],z[12][16],w[12][16]; } trace_t;
static void aes_enc_trace(const aesctx*ctx,const u8 P[16],u8 C[16],trace_t*T){
    u8 s[16]; st_cpy(s,P);
    for(int i=0;i<16;i++) s[i]^=ctx->rk[0][i];
    for(int r=1;r<=ctx->nr;r++){
        if(T) st_cpy(T->x[r],s);
        sub_bytes(s); if(T) st_cpy(T->y[r],s);
        shift_rows(s); if(T) st_cpy(T->z[r],s);
        if(r<ctx->nr){ mix_cols(s); if(T) st_cpy(T->w[r],s); }
        else { if(T) st_cpy(T->w[r],s); } /* w_nr = z_nr */
        for(int i=0;i<16;i++) s[i]^=ctx->rk[r][i];
    }
    st_cpy(C,s);
}
static inline void aes_enc(const aesctx*ctx,const u8 P[16],u8 C[16]){aes_enc_trace(ctx,P,C,NULL);}

/* ---------- position sets ---------- */
static const int DIAG0[4]={0,5,10,15};
static const int ADIAG[4]={0,13,10,7};   /* (r,(4-r)%4) */
static const int P_r[4][4]={            /* P_r[r][i]=pos of (i,(r+i)%4) */
    {0,5,10,15},{4,9,14,3},{8,13,2,7},{12,1,6,11}};

/* ---------- DDT ---------- */
static u8 DDTsol[256][256][4]; static u8 DDTn[256][256];
static void ddt_init(void){
    for(int dx=0;dx<256;dx++)for(int dy=0;dy<256;dy++)DDTn[dx][dy]=0;
    for(int x=0;x<256;x++)for(int dx=0;dx<256;dx++){
        u8 dy=SBOX[x]^SBOX[x^dx];
        if(DDTn[dx][dy]<4) DDTsol[dx][dy][DDTn[dx][dy]]=x;
        DDTn[dx][dy]++;
    }
}

#endif
