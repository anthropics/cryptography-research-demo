// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0

#ifndef LEA_H
#define LEA_H
#include <stdint.h>
#include <string.h>

static const uint32_t LEA_DELTA[4] = {0xC3EFE9DBu, 0x44626B02u, 0x79E27C8Au, 0x78DF30ECu};

static inline uint32_t rol32(uint32_t x, int r){ r&=31; return r? ((x<<r)|(x>>(32-r))) : x; }
static inline uint32_t ror32(uint32_t x, int r){ r&=31; return r? ((x>>r)|(x<<(32-r))) : x; }

/* round keys: rk[i][0..3] = t0,t1,t2,t3  (LEA rk6 = (t0,t1,t2,t1,t3,t1)) */
typedef struct { uint32_t rk[32][4]; int rounds; } lea_ctx;

static inline void lea_setkey(lea_ctx *c, const uint32_t k[4], int rounds){
    uint32_t t[4] = {k[0],k[1],k[2],k[3]};
    c->rounds = rounds;
    for(int i=0;i<rounds;i++){
        uint32_t d = LEA_DELTA[i&3];
        t[0] = rol32(t[0] + rol32(d, i  ), 1);
        t[1] = rol32(t[1] + rol32(d, i+1), 3);
        t[2] = rol32(t[2] + rol32(d, i+2), 6);
        t[3] = rol32(t[3] + rol32(d, i+3), 11);
        c->rk[i][0]=t[0]; c->rk[i][1]=t[1]; c->rk[i][2]=t[2]; c->rk[i][3]=t[3];
    }
}

static inline void lea_enc(const lea_ctx *c, uint32_t x[4]){
    for(int i=0;i<c->rounds;i++){
        const uint32_t *rk = c->rk[i];
        uint32_t y0 = rol32((x[0]^rk[0]) + (x[1]^rk[1]), 9);
        uint32_t y1 = ror32((x[1]^rk[2]) + (x[2]^rk[1]), 5);
        uint32_t y2 = ror32((x[2]^rk[3]) + (x[3]^rk[1]), 3);
        x[3]=x[0]; x[0]=y0; x[1]=y1; x[2]=y2;
    }
}

static inline void lea_enc_nr(const lea_ctx *c, uint32_t x[4], int nr){
    for(int i=0;i<nr;i++){
        const uint32_t *rk = c->rk[i];
        uint32_t y0 = rol32((x[0]^rk[0]) + (x[1]^rk[1]), 9);
        uint32_t y1 = ror32((x[1]^rk[2]) + (x[2]^rk[1]), 5);
        uint32_t y2 = ror32((x[2]^rk[3]) + (x[3]^rk[1]), 3);
        x[3]=x[0]; x[0]=y0; x[1]=y1; x[2]=y2;
    }
}

static inline void lea_dec_nr(const lea_ctx *c, uint32_t x[4], int nr){
    /* decrypt rounds nr-1 down to 0 */
    for(int i=nr-1;i>=0;i--){
        const uint32_t *rk = c->rk[i];
        uint32_t p0 = x[3];
        uint32_t p1 = (ror32(x[0],9) - (p0^rk[0])) ^ rk[1];
        uint32_t p2 = (rol32(x[1],5) - (p1^rk[2])) ^ rk[1];
        uint32_t p3 = (rol32(x[2],3) - (p2^rk[3])) ^ rk[1];
        x[0]=p0;x[1]=p1;x[2]=p2;x[3]=p3;
    }
}
#endif
