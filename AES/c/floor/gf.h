// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// GF(256) + AES S-box toolkit for floor-proof verification
#ifndef GF_H
#define GF_H
#include <stdint.h>
#include <string.h>

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

static u8 LOG[256], EXP[510], INV[256], SB[256], iSB[256];
static u8 MULT[256][256];

static inline u8 xt(u8 a){ return (a<<1) ^ ((a>>7)*0x1b); }

static void gf_init(void){
    u8 a=1;
    for(int i=0;i<255;i++){ EXP[i]=a; LOG[a]=i; a^=xt(a); }
    for(int i=255;i<510;i++) EXP[i]=EXP[i-255];
    EXP[255]=1; // wrap
    LOG[0]=0; INV[0]=0;
    for(int i=1;i<256;i++) INV[i]=EXP[255-LOG[i]];
    for(int i=0;i<256;i++) for(int j=0;j<256;j++)
        MULT[i][j] = (i&&j) ? EXP[LOG[i]+LOG[j]] : 0;
    // AES S-box
    for(int i=0;i<256;i++){
        u8 x = INV[i], r=0;
        for(int k=0;k<5;k++){ r ^= x; x = (x<<1)|(x>>7); }
        SB[i] = r ^ 0x63;
        iSB[SB[i]] = i;
    }
}

static inline u8 gmul(u8 a, u8 b){ return MULT[a][b]; }
static inline u8 ginv(u8 a){ return INV[a]; }
static inline u8 gpow(u8 a, int k){
    if(a==0) return (k==0)?1:0;
    return EXP[((long)LOG[a]*k) % 255];
}

// xorshift64 PRNG
static u64 rng_s;
static inline u64 rnd(void){
    rng_s ^= rng_s<<13; rng_s ^= rng_s>>7; rng_s ^= rng_s<<17;
    return rng_s;
}
static inline u8 rnd8(void){ return (u8)(rnd()>>40); }
static inline u8 rnd8nz(void){ u8 x; do{x=rnd8();}while(!x); return x; }

#endif
