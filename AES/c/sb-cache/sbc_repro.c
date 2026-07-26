// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* sbc_repro.c — Verification of the W/UU S-box cache (XOR-separability of e_ω)
 *
 * Verifies the claim e_ω = ⊕_r U_r^{b_r}[ω] (z_4-XOR-separable),
 * hence 8 U-tables per x_3-config → 16 z_4-combos with 0 SB/ω.
 *
 * Derived from Derbez--Fouque--Jean'13 Prop. 2.
 *
 * Tests per random base ("key"):
 *   (a) ddtgray_walk  : BRGC 2^20, 1.0625 SB/(entry,ω) expected
 *   (b) sbcache_walk  : x_3-gray × 16 z_4, 0.1875 SB/(entry,ω) expected
 *   (c) fp[a][i] == fp[b][i] for all i∈[0,2^20)
 *   (d) 1024 random entries: fp matches cold_e() reference (24 SB/ω)
 *
 * Output: per-key c_SB, match counts; summary over N_KEYS.
 *
 * Build:  gcc -O3 -march=native -o sbc_repro sbc_repro.c
 * Run:    ./sbc_repro [seed] [nkeys]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

/* ==================================================================== */
/* AES S-box + GF(256)                                                  */
/* ==================================================================== */

static const uint8_t SB[256] = {
  0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
  0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
  0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
  0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
  0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
  0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
  0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
  0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
  0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
  0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
  0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
  0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
  0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
  0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
  0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
  0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16
};

static uint8_t iSB[256];
static uint8_t GLOG[256], GEXP[510];     /* log/exp base generator 3 */
static uint8_t GINV[256];
static uint8_t GMUL2[256], GMUL3[256];   /* for e_ω assembly */

static const uint8_t MC_[4][4]  = {{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
static const uint8_t iMC_[4][4] = {{0x0e,0x0b,0x0d,0x09},{0x09,0x0e,0x0b,0x0d},
                                   {0x0d,0x09,0x0e,0x0b},{0x0b,0x0d,0x09,0x0e}};
static uint8_t MCmul[4][4][256];         /* MCmul[r][c][x] = MC[r][c]·x */

/* POWINV: 12 Frob-rep power-sums of d=1/e */
#define NREP 12
static const int REPS[NREP] = {1,3,5,7,9,11,13,17,19,21,25,37};
static uint8_t POWINV[256][NREP];

static inline uint8_t gmul(uint8_t a, uint8_t b){
    if(!a||!b) return 0;
    return GEXP[(int)GLOG[a]+(int)GLOG[b]];
}
static inline uint8_t gpow(uint8_t a, int k){
    if(!a) return 0;
    return GEXP[((long)GLOG[a]*(k%255))%255];
}

static void init_tables(void){
    for(int i=0;i<256;i++) iSB[SB[i]]=(uint8_t)i;
    /* log/exp over GF(256), poly 0x11b, gen 3 */
    uint8_t x=1;
    for(int i=0;i<255;i++){
        GEXP[i]=x; GLOG[x]=(uint8_t)i;
        uint8_t x2=(uint8_t)((x<<1)^((x&0x80)?0x1b:0));
        x ^= x2; /* ×3 */
    }
    for(int i=255;i<510;i++) GEXP[i]=GEXP[i-255];
    GLOG[0]=0; /* unused */
    GINV[0]=0; for(int i=1;i<256;i++) GINV[i]=GEXP[255-GLOG[i]];
    for(int i=0;i<256;i++){ GMUL2[i]=gmul(2,(uint8_t)i); GMUL3[i]=gmul(3,(uint8_t)i); }
    for(int r=0;r<4;r++)for(int c=0;c<4;c++)for(int i=0;i<256;i++)
        MCmul[r][c][i]=gmul(MC_[r][c],(uint8_t)i);
    for(int e=0;e<256;e++){
        uint8_t d=GINV[e];
        for(int j=0;j<NREP;j++) POWINV[e][j]=gpow(d,REPS[j]);
    }
}

/* ==================================================================== */
/* Instrumented SB                                                      */
/* ==================================================================== */

static uint64_t sb_ctr;
static inline uint8_t SBc(uint8_t x){ sb_ctr++; return SB[x]; }

/* ==================================================================== */
/* Position bookkeeping                                                 */
/* ==================================================================== */

#define POS(row,col) ((row)+4*(col))
static int r_of_p[16], row_of_p[16], coef_p[16];
static int adiag[4][4];     /* adiag[r][i] = POS(i,(r+i)%4) */
static int ic_of_col[4];    /* i_c=(4-c)%4 */

static void init_pos(void){
    for(int p=0;p<16;p++){
        int row=p&3, col=p>>2;
        row_of_p[p]=row;
        r_of_p[p]=(col-row)&3;
        coef_p[p]=MC_[r_of_p[p]][row];
    }
    for(int r=0;r<4;r++)for(int i=0;i<4;i++) adiag[r][i]=POS(i,(r+i)&3);
    for(int c=0;c<4;c++) ic_of_col[c]=(4-c)&3;
}

/* ==================================================================== */
/* PRNG (xorshift64*)                                                   */
/* ==================================================================== */

static uint64_t rng_s;
static inline uint64_t rng(void){
    rng_s^=rng_s>>12; rng_s^=rng_s<<25; rng_s^=rng_s>>27;
    return rng_s*0x2545F4914F6CDD1DULL;
}
static inline uint8_t rng8(void){ return (uint8_t)(rng()>>56); }
static inline uint8_t rng8nz(void){ uint8_t x; do x=rng8(); while(!x); return x; }

/* ==================================================================== */
/* Base (one DFJ-10-tuple-like configuration)                           */
/* ==================================================================== */

typedef struct {
    uint8_t x2[4], d1, d5, z4[4];
    /* derived */
    uint8_t dz4[4];          /* iMC[r,0]·d5 */
    uint8_t v[4][2];         /* x_4[r,r] pair: iSB(z4[r]), iSB(z4[r]⊕dz4[r]) */
    uint8_t sbv[4][2];       /* SB(v[r][b]) — ω-indep */
    uint8_t x3[16][2];       /* x_3[p] pair */
    uint8_t sbx3[16][2];     /* SB(x_3[p]) — ω-indep */
    uint8_t dx3[256][16];    /* Δx_3[p]^{(ω)} — base-fixed */
} Base;

/* DDT first-solution table: DDTSOL[din][dout] = first x with SB(x)^SB(x^din)=dout, or 0xFF sentinel; DDTN = count */
static uint8_t DDTSOL[256][256], DDTN[256][256];
static void init_ddt(void){
    memset(DDTN,0,sizeof DDTN);
    memset(DDTSOL,0xFF,sizeof DDTSOL);
    for(int din=1;din<256;din++) for(int x=0;x<256;x++){
        uint8_t dout=(uint8_t)(SB[x]^SB[x^din]);
        if(DDTN[din][dout]==0) DDTSOL[din][dout]=(uint8_t)x;
        DDTN[din][dout]++;
    }
}

/* mode 0: random x_3 pairs (structural test)
 * mode 1: rebound-valid x_3 pairs (DDT solutions) */
static int gen_base(Base *B, int valid_mode){
    B->d1=rng8nz(); B->d5=rng8nz();
    for(int i=0;i<4;i++) B->x2[i]=rng8();
    for(int r=0;r<4;r++) B->z4[r]=rng8();
    for(int r=0;r<4;r++){
        B->dz4[r]=gmul(iMC_[r][0],B->d5);
        B->v[r][0]=iSB[B->z4[r]];
        B->v[r][1]=iSB[(uint8_t)(B->z4[r]^B->dz4[r])];
        B->sbv[r][0]=SB[B->v[r][0]];
        B->sbv[r][1]=SB[B->v[r][1]];
    }
    /* precompute dx3[ω][p]; needs Δy_2[i]^{(ω)} */
    uint8_t sbx2[4]; for(int i=0;i<4;i++) sbx2[i]=SB[B->x2[i]];
    for(int w=0;w<256;w++){
        uint8_t dy2[4];
        for(int i=0;i<4;i++)
            dy2[i]=(uint8_t)(SB[(uint8_t)(B->x2[i]^gmul(MC_[i][0],(uint8_t)w))]^sbx2[i]);
        for(int p=0;p<16;p++){
            int row=p&3, col=p>>2, ic=ic_of_col[col];
            B->dx3[w][p]=MCmul[row][ic][dy2[ic]];
        }
    }
    if(!valid_mode){
        /* x_3 pairs: random base, partner = ⊕ Δx_3[p]^{(δ1)} */
        for(int p=0;p<16;p++){
            B->x3[p][0]=rng8();
            B->x3[p][1]=(uint8_t)(B->x3[p][0]^B->dx3[B->d1][p]);
            B->sbx3[p][0]=SB[B->x3[p][0]];
            B->sbx3[p][1]=SB[B->x3[p][1]];
        }
        return 1;
    }
    /* rebound-valid: for each r, search z_4[r] s.t. all 4 adiag_r DDT entries nonzero */
    for(int r=0;r<4;r++){
        int ok=0;
        for(int tries=0;tries<256;tries++){
            uint8_t z=rng8();
            uint8_t v0=iSB[z], v1=iSB[(uint8_t)(z^B->dz4[r])];
            uint8_t dx4=(uint8_t)(v0^v1);
            if(!dx4) continue;
            int allok=1;
            uint8_t x3v[4];
            for(int i=0;i<4;i++){
                int p=adiag[r][i];
                uint8_t din=B->dx3[B->d1][p];
                uint8_t dout=gmul(iMC_[i][r],dx4); /* Δy_3[i,(r+i)%4]=iMC[i,r]·Δx_4[r,r] */
                if(DDTN[din][dout]==0){allok=0;break;}
                x3v[i]=DDTSOL[din][dout];
            }
            if(!allok) continue;
            /* accept */
            B->z4[r]=z; B->v[r][0]=v0; B->v[r][1]=v1;
            B->sbv[r][0]=SB[v0]; B->sbv[r][1]=SB[v1];
            for(int i=0;i<4;i++){
                int p=adiag[r][i];
                B->x3[p][0]=x3v[i];
                B->x3[p][1]=(uint8_t)(x3v[i]^B->dx3[B->d1][p]);
                B->sbx3[p][0]=SB[B->x3[p][0]];
                B->sbx3[p][1]=SB[B->x3[p][1]];
            }
            ok=1; break;
        }
        if(!ok) return 0; /* caller retries whole base */
    }
    return 1;
}

/* ==================================================================== */
/* Fingerprint                                                          */
/* ==================================================================== */

static inline uint64_t fp_of_e(const uint8_t e[256]){
    uint64_t h=0xcbf29ce484222325ULL;
    uint8_t Sk[NREP]={0};
    for(int w=1;w<256;w++){
        h=(h^e[w])*0x100000001b3ULL;
        const uint8_t *pv=POWINV[e[w]];
        for(int j=0;j<NREP;j++) Sk[j]^=pv[j];
    }
    for(int j=0;j<NREP;j++) h=(h^Sk[j])*0x100000001b3ULL;
    return h;
}

/* ==================================================================== */
/* COLD reference: e_ω from scratch (24 SB/ω, no shortcuts)             */
/* ==================================================================== */

static void cold_e(const Base *B, uint32_t x3b, uint32_t z4b, uint8_t e[256]){
    uint8_t x3v[16], sbx3v[16], x4v[4], sbx4v[4];
    for(int p=0;p<16;p++){int b=(x3b>>p)&1; x3v[p]=B->x3[p][b]; sbx3v[p]=B->sbx3[p][b];}
    for(int r=0;r<4;r++){int b=(z4b>>r)&1;  x4v[r]=B->v[r][b];  sbx4v[r]=B->sbv[r][b];}
    e[0]=0;
    for(int w=1;w<256;w++){
        uint8_t dy3[16];
        for(int p=0;p<16;p++)
            dy3[p]=(uint8_t)(SB[(uint8_t)(x3v[p]^B->dx3[w][p])]^sbx3v[p]);
        uint8_t dy4[4];
        for(int r=0;r<4;r++){
            uint8_t dx4=0;
            for(int i=0;i<4;i++) dx4^=MCmul[r][i][dy3[adiag[r][i]]];
            dy4[r]=(uint8_t)(SB[(uint8_t)(x4v[r]^dx4)]^sbx4v[r]);
        }
        e[w]=(uint8_t)(GMUL2[dy4[0]]^GMUL3[dy4[1]]^dy4[2]^dy4[3]);
    }
}

/* ==================================================================== */
/* DDT-GRAY walk: BRGC over 20 bits, z_4 innermost                      */
/* ==================================================================== */

static void ddtgray_walk(const Base *B, uint64_t *fp_out, uint64_t *sb_out){
    /* state arrays */
    static uint8_t dy3[16][256], dx4[4][256], dy4[4][256], e[256];
    uint8_t curx3[16], cursbx3[16], curx4[4], cursbx4[4];
    uint32_t x3b=0, z4b=0;

    /* cold init at (0,0) — NOT counted in sb_ctr (amortized init) */
    for(int p=0;p<16;p++){curx3[p]=B->x3[p][0]; cursbx3[p]=B->sbx3[p][0];}
    for(int r=0;r<4;r++){curx4[r]=B->v[r][0];   cursbx4[r]=B->sbv[r][0];}
    for(int w=0;w<256;w++){
        for(int p=0;p<16;p++) dy3[p][w]=(uint8_t)(SB[(uint8_t)(curx3[p]^B->dx3[w][p])]^cursbx3[p]);
        for(int r=0;r<4;r++){
            uint8_t a=0; for(int i=0;i<4;i++) a^=MCmul[r][i][dy3[adiag[r][i]][w]];
            dx4[r][w]=a;
            dy4[r][w]=(uint8_t)(SB[(uint8_t)(curx4[r]^a)]^cursbx4[r]);
        }
        e[w]=(uint8_t)(GMUL2[dy4[0][w]]^GMUL3[dy4[1][w]]^dy4[2][w]^dy4[3][w]);
    }
    sb_ctr=0;
    fp_out[0]=fp_of_e(e);

    for(uint32_t t=1;t<(1u<<20);t++){
        uint32_t g=t^(t>>1), gp=(t-1)^((t-1)>>1), bit=__builtin_ctz(g^gp);
        if(bit<4){
            /* z_4[bit] flip */
            int r=bit; z4b^=(1u<<r); int b=(z4b>>r)&1;
            curx4[r]=B->v[r][b]; cursbx4[r]=B->sbv[r][b];
            uint8_t cv=curx4[r], sv=cursbx4[r];
            const uint8_t *MCR=(r==0?GMUL2:(r==1?GMUL3:NULL));
            for(int w=1;w<256;w++){
                uint8_t nd=(uint8_t)(SBc((uint8_t)(cv^dx4[r][w]))^sv);
                uint8_t dd=(uint8_t)(nd^dy4[r][w]);
                dy4[r][w]=nd;
                e[w]^=(MCR?MCR[dd]:dd);
            }
        } else {
            /* x_3[p] flip, p=bit-4 */
            int p=bit-4; x3b^=(1u<<p); int b=(x3b>>p)&1;
            curx3[p]=B->x3[p][b]; cursbx3[p]=B->sbx3[p][b];
            int r=r_of_p[p]; uint8_t coef=coef_p[p];
            uint8_t cx3=curx3[p], sx3=cursbx3[p];
            uint8_t cx4=curx4[r], sx4=cursbx4[r];
            const uint8_t *MCR=(r==0?GMUL2:(r==1?GMUL3:NULL));
            for(int w=1;w<256;w++){
                uint8_t nd3=(uint8_t)(SBc((uint8_t)(cx3^B->dx3[w][p]))^sx3);
                uint8_t dd3=(uint8_t)(nd3^dy3[p][w]); dy3[p][w]=nd3;
                dx4[r][w]^=gmul(coef,dd3);
                uint8_t nd4=(uint8_t)(SBc((uint8_t)(cx4^dx4[r][w]))^sx4);
                uint8_t dd4=(uint8_t)(nd4^dy4[r][w]); dy4[r][w]=nd4;
                e[w]^=(MCR?MCR[dd4]:dd4);
            }
        }
        uint32_t idx=(x3b<<4)|z4b;
        fp_out[idx]=fp_of_e(e);
    }
    *sb_out=sb_ctr;
}

/* ==================================================================== */
/* SB-CACHE walk: x_3-gray outer, 16 z_4-combos inner via U-tables      */
/* ==================================================================== */

static void sbcache_walk(const Base *B, uint64_t *fp_out, uint64_t *sb_out){
    static uint8_t dy3[16][256], dx4[4][256];
    static uint8_t U[4][2][256];   /* U_r^b[ω] = MC[0,r]·(SB(v_r^b⊕dx4[r])⊕SB(v_r^b)) */
    uint8_t curx3[16], cursbx3[16];
    uint32_t x3b=0;

    /* cold init at x3b=0 — NOT counted */
    for(int p=0;p<16;p++){curx3[p]=B->x3[p][0]; cursbx3[p]=B->sbx3[p][0];}
    for(int w=0;w<256;w++){
        for(int p=0;p<16;p++) dy3[p][w]=(uint8_t)(SB[(uint8_t)(curx3[p]^B->dx3[w][p])]^cursbx3[p]);
        for(int r=0;r<4;r++){
            uint8_t a=0; for(int i=0;i<4;i++) a^=MCmul[r][i][dy3[adiag[r][i]][w]];
            dx4[r][w]=a;
            for(int b=0;b<2;b++){
                uint8_t dy4=(uint8_t)(SB[(uint8_t)(B->v[r][b]^a)]^B->sbv[r][b]);
                U[r][b][w]=(r==0?GMUL2[dy4]:(r==1?GMUL3[dy4]:dy4));
            }
        }
    }
    sb_ctr=0;

    for(uint32_t tx=0;;tx++){
        /* assemble 16 z_4-configs for current x_3-config */
        uint8_t e[256];
        for(uint32_t z4b=0;z4b<16;z4b++){
            int b0=z4b&1,b1=(z4b>>1)&1,b2=(z4b>>2)&1,b3=(z4b>>3)&1;
            e[0]=0;
            for(int w=1;w<256;w++)
                e[w]=(uint8_t)(U[0][b0][w]^U[1][b1][w]^U[2][b2][w]^U[3][b3][w]);
            uint32_t idx=(x3b<<4)|z4b;
            fp_out[idx]=fp_of_e(e);
        }
        if(tx==(1u<<16)-1) break;

        /* x_3-gray step: flip bit p */
        uint32_t t=tx+1, g=t^(t>>1), gp=tx^(tx>>1);
        int p=__builtin_ctz(g^gp);
        x3b^=(1u<<p); int b=(x3b>>p)&1;
        curx3[p]=B->x3[p][b]; cursbx3[p]=B->sbx3[p][b];
        int r=r_of_p[p]; uint8_t coef=coef_p[p];
        uint8_t cx3=curx3[p], sx3=cursbx3[p];
        for(int w=1;w<256;w++){
            uint8_t nd3=(uint8_t)(SBc((uint8_t)(cx3^B->dx3[w][p]))^sx3);
            uint8_t dd3=(uint8_t)(nd3^dy3[p][w]); dy3[p][w]=nd3;
            uint8_t ndx=(uint8_t)(dx4[r][w]^gmul(coef,dd3)); dx4[r][w]=ndx;
            /* rebuild BOTH U_r^0, U_r^1 for this ω */
            for(int bb=0;bb<2;bb++){
                uint8_t dy4=(uint8_t)(SBc((uint8_t)(B->v[r][bb]^ndx))^B->sbv[r][bb]);
                U[r][bb][w]=(r==0?GMUL2[dy4]:(r==1?GMUL3[dy4]:dy4));
            }
        }
    }
    *sb_out=sb_ctr;
}

/* ==================================================================== */
/* main                                                                 */
/* ==================================================================== */

#define NENT (1u<<20)

/* ==================================================================== */
/* UU-CACHE walk: FULL 20-bit separability                              */
/*                                                                      */
/* Observation: the 20 bits partition into 4 DISJOINT groups of 5:     */
/*   G_r = {x_3-bits on adiag_r} ∪ {z_4-bit r},  r=0..3                */
/* and e_ω = ⊕_r UU_r[g_r][ω], where UU_r depends ONLY on group-r bits.*/
/* Proof: Δ_ω^{(r)}=⊕_{p∈adiag_r} coef_p·dy3_p^{(ω)} depends only on   */
/* x_3-bits in adiag_r; U_r^{b_r}[ω]=f(Δ_ω^{(r)},v_r^{b_r}) depends    */
/* only on those + z_4[r]. QED.                                        */
/*                                                                      */
/* Precompute UU[4][32][256] once per base (≤4·32·255 ≈ 32640 SB).     */
/* Then ALL 2^20 entries: 0 SB, 4 load + 3 XOR per (entry,ω).          */
/* c_SB ≈ 32640/(2^20·255) ≈ 1.2e-4.                                   */
/* ==================================================================== */

static int adiag_bitpos[4][4]; /* which x_3-bit index (0..15) is adiag[r][i] */
static void init_uu_map(void){
    for(int r=0;r<4;r++) for(int i=0;i<4;i++) adiag_bitpos[r][i]=adiag[r][i];
}
static inline uint32_t group_idx(uint32_t x3b,uint32_t z4b,int r){
    uint32_t g=0;
    for(int i=0;i<4;i++) g|=((x3b>>adiag_bitpos[r][i])&1u)<<i;
    g|=((z4b>>r)&1u)<<4;
    return g;
}

static void uucache_walk(const Base *B, uint64_t *fp_out, uint64_t *sb_out){
    static uint8_t UU[4][32][256];
    /* --- precompute (counted) --- */
    sb_ctr=0;
    /* dy3 both branches: 32 tables */
    static uint8_t DY3[16][2][256];
    for(int p=0;p<16;p++) for(int b=0;b<2;b++){
        uint8_t xv=B->x3[p][b], sv=B->sbx3[p][b];
        for(int w=1;w<256;w++)
            DY3[p][b][w]=(uint8_t)(SBc((uint8_t)(xv^B->dx3[w][p]))^sv);
    }
    /* UU[r][g][ω] */
    for(int r=0;r<4;r++){
        for(int g=0;g<32;g++){
            int bz=(g>>4)&1; uint8_t vv=B->v[r][bz], sv=B->sbv[r][bz];
            for(int w=1;w<256;w++){
                uint8_t dx4=0;
                for(int i=0;i<4;i++){
                    int p=adiag[r][i], bi=(g>>i)&1;
                    dx4^=gmul((uint8_t)coef_p[p],DY3[p][bi][w]);
                }
                uint8_t dy4=(uint8_t)(SBc((uint8_t)(vv^dx4))^sv);
                UU[r][g][w]=(r==0?GMUL2[dy4]:(r==1?GMUL3[dy4]:dy4));
            }
        }
    }
    uint64_t sb_pre=sb_ctr;

    /* --- enumerate all 2^20 entries: structured as SB-CACHE (x_3 outer, z_4 inner) --- */
    /* For each x_3-config, select 8 UU-rows (2 per r: bz=0,1), then 16 z_4-combos.    */
    for(uint32_t x3b=0;x3b<(1u<<16);x3b++){
        const uint8_t *R[4][2];
        for(int r=0;r<4;r++){
            uint32_t gx=0; for(int i=0;i<4;i++) gx|=((x3b>>adiag_bitpos[r][i])&1u)<<i;
            R[r][0]=UU[r][gx]; R[r][1]=UU[r][gx|16];
        }
        uint8_t e[256];
        for(uint32_t z4b=0;z4b<16;z4b++){
            int b0=z4b&1,b1=(z4b>>1)&1,b2=(z4b>>2)&1,b3=(z4b>>3)&1;
            const uint8_t *p0=R[0][b0],*p1=R[1][b1],*p2=R[2][b2],*p3=R[3][b3];
            e[0]=0;
            for(int w=1;w<256;w++) e[w]=(uint8_t)(p0[w]^p1[w]^p2[w]^p3[w]);
            fp_out[(x3b<<4)|z4b]=fp_of_e(e);
        }
    }
    *sb_out=sb_pre; /* zero SB in enumeration phase */
}

static uint64_t bench_uucache(const Base *B){
    static uint8_t UU[4][32][256];
    static uint8_t DY3[16][2][256];
    for(int p=0;p<16;p++) for(int b=0;b<2;b++){
        uint8_t xv=B->x3[p][b], sv=B->sbx3[p][b];
        for(int w=1;w<256;w++) DY3[p][b][w]=(uint8_t)(SB[(uint8_t)(xv^B->dx3[w][p])]^sv);
    }
    for(int r=0;r<4;r++) for(int g=0;g<32;g++){
        int bz=(g>>4)&1; uint8_t vv=B->v[r][bz], sv=B->sbv[r][bz];
        for(int w=1;w<256;w++){
            uint8_t dx4=0;
            for(int i=0;i<4;i++){int p=adiag[r][i]; dx4^=gmul((uint8_t)coef_p[p],DY3[p][(g>>i)&1][w]);}
            uint8_t dy4=(uint8_t)(SB[(uint8_t)(vv^dx4)]^sv);
            UU[r][g][w]=(r==0?GMUL2[dy4]:(r==1?GMUL3[dy4]:dy4));
        }
    }
    uint64_t chk=0;
    for(uint32_t x3b=0;x3b<(1u<<16);x3b++){
        const uint8_t *R[4][2];
        for(int r=0;r<4;r++){
            uint32_t gx=0; for(int i=0;i<4;i++) gx|=((x3b>>adiag_bitpos[r][i])&1u)<<i;
            R[r][0]=UU[r][gx]; R[r][1]=UU[r][gx|16];
        }
        for(int w=1;w<256;w++){
            uint8_t u0a=R[0][0][w],u0b=R[0][1][w],u1a=R[1][0][w],u1b=R[1][1][w];
            uint8_t u2a=R[2][0][w],u2b=R[2][1][w],u3a=R[3][0][w],u3b=R[3][1][w];
            chk+=(u0a^u1a^u2a^u3a); chk+=(u0b^u1a^u2a^u3a);
            chk+=(u0a^u1b^u2a^u3a); chk+=(u0b^u1b^u2a^u3a);
            chk+=(u0a^u1a^u2b^u3a); chk+=(u0b^u1a^u2b^u3a);
            chk+=(u0a^u1b^u2b^u3a); chk+=(u0b^u1b^u2b^u3a);
            chk+=(u0a^u1a^u2a^u3b); chk+=(u0b^u1a^u2a^u3b);
            chk+=(u0a^u1b^u2a^u3b); chk+=(u0b^u1b^u2a^u3b);
            chk+=(u0a^u1a^u2b^u3b); chk+=(u0b^u1a^u2b^u3b);
            chk+=(u0a^u1b^u2b^u3b); chk+=(u0b^u1b^u2b^u3b);
        }
    }
    return chk;
}

/* ==================================================================== */
/* e-only benchmark (no fp): measure raw build cost                     */
/* ==================================================================== */

static uint64_t bench_ddtgray(const Base *B){
    static uint8_t dy3[16][256], dx4[4][256], dy4[4][256], e[256];
    uint8_t curx3[16], cursbx3[16], curx4[4], cursbx4[4];
    uint32_t x3b=0, z4b=0; uint64_t chk=0;
    for(int p=0;p<16;p++){curx3[p]=B->x3[p][0]; cursbx3[p]=B->sbx3[p][0];}
    for(int r=0;r<4;r++){curx4[r]=B->v[r][0];   cursbx4[r]=B->sbv[r][0];}
    for(int w=0;w<256;w++){
        for(int p=0;p<16;p++) dy3[p][w]=(uint8_t)(SB[(uint8_t)(curx3[p]^B->dx3[w][p])]^cursbx3[p]);
        for(int r=0;r<4;r++){
            uint8_t a=0; for(int i=0;i<4;i++) a^=MCmul[r][i][dy3[adiag[r][i]][w]];
            dx4[r][w]=a; dy4[r][w]=(uint8_t)(SB[(uint8_t)(curx4[r]^a)]^cursbx4[r]);
        }
        e[w]=(uint8_t)(GMUL2[dy4[0][w]]^GMUL3[dy4[1][w]]^dy4[2][w]^dy4[3][w]);
    }
    for(int w=1;w<256;w++) chk+=e[w];
    for(uint32_t t=1;t<(1u<<20);t++){
        int bit=__builtin_ctz(t);
        if(bit<4){
            int r=bit; z4b^=(1u<<r); int b=(z4b>>r)&1;
            uint8_t cv=B->v[r][b], sv=B->sbv[r][b];
            const uint8_t *MCR=(r==0?GMUL2:(r==1?GMUL3:NULL));
            for(int w=1;w<256;w++){
                uint8_t nd=(uint8_t)(SB[(uint8_t)(cv^dx4[r][w])]^sv);
                uint8_t dd=(uint8_t)(nd^dy4[r][w]); dy4[r][w]=nd;
                e[w]^=(MCR?MCR[dd]:dd);
            }
        } else {
            int p=bit-4; x3b^=(1u<<p); int b=(x3b>>p)&1;
            int r=r_of_p[p]; uint8_t coef=coef_p[p];
            uint8_t cx3=B->x3[p][b], sx3=B->sbx3[p][b];
            uint8_t cx4=B->v[r][(z4b>>r)&1], sx4=B->sbv[r][(z4b>>r)&1];
            const uint8_t *MCR=(r==0?GMUL2:(r==1?GMUL3:NULL));
            for(int w=1;w<256;w++){
                uint8_t nd3=(uint8_t)(SB[(uint8_t)(cx3^B->dx3[w][p])]^sx3);
                uint8_t dd3=(uint8_t)(nd3^dy3[p][w]); dy3[p][w]=nd3;
                dx4[r][w]^=gmul(coef,dd3);
                uint8_t nd4=(uint8_t)(SB[(uint8_t)(cx4^dx4[r][w])]^sx4);
                uint8_t dd4=(uint8_t)(nd4^dy4[r][w]); dy4[r][w]=nd4;
                e[w]^=(MCR?MCR[dd4]:dd4);
            }
        }
        for(int w=1;w<256;w++) chk+=e[w]; /* prevent DCE; mimics 1 op/ω read */
    }
    return chk;
}

static uint64_t bench_sbcache(const Base *B){
    static uint8_t dy3[16][256], dx4[4][256], U[4][2][256];
    uint32_t x3b=0; uint64_t chk=0;
    for(int w=0;w<256;w++){
        for(int p=0;p<16;p++) dy3[p][w]=(uint8_t)(SB[(uint8_t)(B->x3[p][0]^B->dx3[w][p])]^B->sbx3[p][0]);
        for(int r=0;r<4;r++){
            uint8_t a=0; for(int i=0;i<4;i++) a^=MCmul[r][i][dy3[adiag[r][i]][w]];
            dx4[r][w]=a;
            for(int b=0;b<2;b++){
                uint8_t dy4=(uint8_t)(SB[(uint8_t)(B->v[r][b]^a)]^B->sbv[r][b]);
                U[r][b][w]=(r==0?GMUL2[dy4]:(r==1?GMUL3[dy4]:dy4));
            }
        }
    }
    for(uint32_t tx=0;;tx++){
        for(int w=1;w<256;w++){
            uint8_t u0a=U[0][0][w],u0b=U[0][1][w],u1a=U[1][0][w],u1b=U[1][1][w];
            uint8_t u2a=U[2][0][w],u2b=U[2][1][w],u3a=U[3][0][w],u3b=U[3][1][w];
            /* 16 combos */
            chk+=(u0a^u1a^u2a^u3a); chk+=(u0b^u1a^u2a^u3a);
            chk+=(u0a^u1b^u2a^u3a); chk+=(u0b^u1b^u2a^u3a);
            chk+=(u0a^u1a^u2b^u3a); chk+=(u0b^u1a^u2b^u3a);
            chk+=(u0a^u1b^u2b^u3a); chk+=(u0b^u1b^u2b^u3a);
            chk+=(u0a^u1a^u2a^u3b); chk+=(u0b^u1a^u2a^u3b);
            chk+=(u0a^u1b^u2a^u3b); chk+=(u0b^u1b^u2a^u3b);
            chk+=(u0a^u1a^u2b^u3b); chk+=(u0b^u1a^u2b^u3b);
            chk+=(u0a^u1b^u2b^u3b); chk+=(u0b^u1b^u2b^u3b);
        }
        if(tx==(1u<<16)-1) break;
        int p=__builtin_ctz(tx+1);
        x3b^=(1u<<p); int b=(x3b>>p)&1;
        int r=r_of_p[p]; uint8_t coef=coef_p[p];
        uint8_t cx3=B->x3[p][b], sx3=B->sbx3[p][b];
        for(int w=1;w<256;w++){
            uint8_t nd3=(uint8_t)(SB[(uint8_t)(cx3^B->dx3[w][p])]^sx3);
            uint8_t dd3=(uint8_t)(nd3^dy3[p][w]); dy3[p][w]=nd3;
            uint8_t ndx=(uint8_t)(dx4[r][w]^gmul(coef,dd3)); dx4[r][w]=ndx;
            for(int bb=0;bb<2;bb++){
                uint8_t dy4=(uint8_t)(SB[(uint8_t)(B->v[r][bb]^ndx)]^B->sbv[r][bb]);
                U[r][bb][w]=(r==0?GMUL2[dy4]:(r==1?GMUL3[dy4]:dy4));
            }
        }
    }
    return chk;
}

int main(int argc, char **argv){
    uint64_t seed = (argc>1)?strtoull(argv[1],0,0):0xC0FFEEULL;
    int nkeys    = (argc>2)?atoi(argv[2]):20;
    int nspot    = (argc>3)?atoi(argv[3]):1024;
    int vmode    = (argc>4)?atoi(argv[4]):0;  /* 0=random x_3, 1=rebound-valid */
    int bench    = (argc>5)?atoi(argv[5]):0;  /* 1=also run e-only bench */

    init_tables(); init_pos(); init_ddt(); init_uu_map();
    rng_s=seed?seed:1;

    uint64_t *fpA=malloc(NENT*sizeof(uint64_t));
    uint64_t *fpB=malloc(NENT*sizeof(uint64_t));
    uint64_t *fpC=malloc(NENT*sizeof(uint64_t));
    if(!fpA||!fpB||!fpC){fprintf(stderr,"oom\n");return 1;}

    printf("# sbc_repro seed=0x%llx nkeys=%d nspot=%d vmode=%s bench=%d\n",
           (unsigned long long)seed,nkeys,nspot,vmode?"REBOUND-VALID":"random-x3",bench);
    printf("# per-key: c_SB(gray) c_SB(cache) fp_match/%u cold_match/%d\n",NENT,nspot);

    uint64_t tot_match=0, tot_matchC=0, tot_spot=0;
    double sum_cg=0, sum_cc=0, sum_cu=0;
    uint64_t tot_sbA=0, tot_sbB=0, tot_sbC=0;
    double sum_bg=0, sum_bc=0, sum_bu=0;

    for(int k=0;k<nkeys;k++){
        static Base B;
        int tries=0; while(!gen_base(&B,vmode)) tries++;
        if(vmode && k==0) printf("# (rebound-valid gen: %d retries)\n",tries);

        clock_t t0=clock();
        uint64_t sbA; ddtgray_walk(&B,fpA,&sbA);
        clock_t t1=clock();
        uint64_t sbB; sbcache_walk(&B,fpB,&sbB);
        clock_t t2=clock();
        uint64_t sbC; uucache_walk(&B,fpC,&sbC);
        clock_t t3=clock();

        /* compare all 2^20 */
        uint64_t nmatch=0, nmatchC=0;
        for(uint32_t i=0;i<NENT;i++){nmatch+=(fpA[i]==fpB[i]); nmatchC+=(fpA[i]==fpC[i]);}

        /* spot-check vs cold */
        uint64_t nspot_ok=0;
        for(int s=0;s<nspot;s++){
            uint32_t idx=(uint32_t)(rng()%NENT);
            uint8_t e[256]; cold_e(&B,idx>>4,idx&15,e);
            uint64_t h=fp_of_e(e);
            if(h==fpA[idx] && h==fpB[idx]) nspot_ok++;
        }

        /* c_SB = SB / (N_entries * 255) */
        double cg=(double)sbA/((double)NENT*255.0);
        double cc=(double)sbB/((double)NENT*255.0);
        double cu=(double)sbC/((double)NENT*255.0);
        sum_cg+=cg; sum_cc+=cc; sum_cu+=cu; tot_sbA+=sbA; tot_sbB+=sbB; tot_sbC+=sbC;
        tot_match+=nmatch; tot_matchC+=nmatchC; tot_spot+=nspot_ok;

        printf("K%02d  c_SB: gray=%.5f cache=%.5f UU=%.6f  match(cache)=%llu match(UU)=%llu  cold=%llu/%d\n",
               k,cg,cc,cu,(unsigned long long)nmatch,(unsigned long long)nmatchC,
               (unsigned long long)nspot_ok,nspot);
        printf("     times: gray=%.2fs cache=%.2fs UU=%.2fs\n",
               (double)(t1-t0)/CLOCKS_PER_SEC,(double)(t2-t1)/CLOCKS_PER_SEC,(double)(t3-t2)/CLOCKS_PER_SEC);

        if(bench){
            clock_t b0=clock(); uint64_t ca=bench_ddtgray(&B); clock_t b1=clock();
            uint64_t cb=bench_sbcache(&B); clock_t b2=clock();
            uint64_t cc2=bench_uucache(&B); clock_t b3=clock();
            double tg=(double)(b1-b0)/CLOCKS_PER_SEC, tc=(double)(b2-b1)/CLOCKS_PER_SEC, tu=(double)(b3-b2)/CLOCKS_PER_SEC;
            double ops=(double)NENT*255.0;
            printf("     bench(e-only): gray=%.3fs(%.2fns/ω) cache=%.3fs(%.2fns/ω) UU=%.3fs(%.2fns/ω)  chk=%s/%s\n",
                   tg,tg*1e9/ops,tc,tc*1e9/ops,tu,tu*1e9/ops,
                   (ca==cb)?"ok":"FAIL",(ca==cc2)?"ok":"FAIL");
            sum_bg+=tg; sum_bc+=tc; sum_bu+=tu;
        }
        fflush(stdout);
    }

    printf("\n# ==== SUMMARY (%d keys, mode=%s) ====\n",nkeys,vmode?"REBOUND-VALID":"random-x3");
    printf("# c_SB(DDT-GRAY)  mean=%.6f   (expect 17/16 = 1.0625)\n",sum_cg/nkeys);
    printf("# c_SB(SB-CACHE)  mean=%.6f   (expect  3/16 = 0.1875)\n",sum_cc/nkeys);
    printf("# c_SB(UU-CACHE)  mean=%.7f  (expect ~160·255/(2^20·255) = 1.5e-4)\n",sum_cu/nkeys);
    printf("# fp match SB-CACHE  %llu / %llu  (%s)\n",
           (unsigned long long)tot_match,(unsigned long long)nkeys*NENT,
           tot_match==(uint64_t)nkeys*NENT?"EXACT":"MISMATCH");
    printf("# fp match UU-CACHE  %llu / %llu  (%s)\n",
           (unsigned long long)tot_matchC,(unsigned long long)nkeys*NENT,
           tot_matchC==(uint64_t)nkeys*NENT?"EXACT":"MISMATCH");
    printf("# cold spot          %llu / %llu  (%s)\n",
           (unsigned long long)tot_spot,(unsigned long long)nkeys*nspot,
           tot_spot==(uint64_t)nkeys*nspot?"EXACT":"MISMATCH");
    printf("# SB-ratio cache/gray = %.4f  UU/gray = %.6f\n",
           (double)tot_sbB/(double)tot_sbA,(double)tot_sbC/(double)tot_sbA);
    if(bench){
        double ops=(double)nkeys*NENT*255.0;
        printf("# e-only wallclock ns/ω:  gray=%.2f  cache=%.2f  UU=%.2f   (speedup %.2fx / %.2fx)\n",
               sum_bg*1e9/ops,sum_bc*1e9/ops,sum_bu*1e9/ops,sum_bg/sum_bc,sum_bg/sum_bu);
    }
    int pass=(tot_match==(uint64_t)nkeys*NENT && tot_matchC==(uint64_t)nkeys*NENT
              && tot_spot==(uint64_t)nkeys*nspot);
    printf("# VERDICT: %s\n",pass?
           "PASS — SB-CACHE & UU-CACHE reproduce DDT-GRAY fp exactly":"FAIL");

    free(fpA); free(fpB); free(fpC);
    return pass?0:1;
}
