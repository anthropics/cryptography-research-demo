// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// bc_impl_b.c — Bridge-Consistency false-accept test, IMPLEMENTATION B (independent)
// GF(256) via direct 64KB multiplication table (carryless reduce), RNG = PCG64-DXSM.
// Structurally different from impl A: vote-per-ω (solve quadratic) then histogram max.
// This is a *different algorithm* computing the *same statistic* — stronger independence.
//
// build: gcc -O3 -march=native -o bc_b bc_impl_b.c
// run:   ./bc_b <Ntrials> <seed_hi> <seed_lo>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- GF(256) via direct table (no logs) ---- */
static uint8_t GMUL[256][256];
static uint8_t GINV[256];
static uint8_t HTR[256]; /* half-trace: HTR[c] is a root of x^2+x=c when Tr(c)=0 */
static uint8_t TR[256];

static uint8_t gf_mul_slow(uint8_t a, uint8_t b){
    uint16_t r=0;
    for(int i=0;i<8;i++){ if(b&1) r^=a; b>>=1;
        a = (a<<1) ^ ((a&0x80)?0x11B:0); }
    return (uint8_t)r;
}
static void gf_init(void){
    for(int a=0;a<256;a++) for(int b=0;b<256;b++) GMUL[a][b]=gf_mul_slow(a,b);
    GINV[0]=0;
    for(int a=1;a<256;a++) for(int b=1;b<256;b++) if(GMUL[a][b]==1){GINV[a]=b;break;}
    /* trace: Tr(x)=x+x^2+x^4+...+x^128 */
    for(int a=0;a<256;a++){
        uint8_t t=0,x=a; for(int i=0;i<8;i++){ t^=x; x=GMUL[x][x]; } TR[a]=t&1;
    }
    /* root table: HTR[c] = any x with x^2+x=c (when Tr(c)=0); brute force */
    for(int c=0;c<256;c++){
        HTR[c]=0;
        for(int x=0;x<256;x++) if((GMUL[x][x]^x)==(uint8_t)c){HTR[c]=(uint8_t)x;break;}
    }
    /* self-check */
    for(int c=0;c<256;c++) if(TR[c]==0){
        uint8_t x=HTR[c]; if((GMUL[x][x]^x)!=(uint8_t)c){fprintf(stderr,"HTR fail c=%d\n",c);exit(1);}
    }
}

/* ---- PCG64-DXSM ---- */
typedef struct { __uint128_t state; } pcg_t;
static pcg_t P;
static inline uint64_t pcg64(void){
    const __uint128_t MUL = ((__uint128_t)2549297995355413924ULL<<64)|4865540595714422341ULL;
    __uint128_t s = P.state;
    P.state = s*MUL + 1442695040888963407ULL;
    uint64_t hi = (uint64_t)(s>>64), lo = (uint64_t)s | 1;
    hi ^= hi>>32; hi *= 0xDA942042E4DD58B5ULL; hi ^= hi>>48; hi *= lo;
    return hi;
}
static void seed_rng(uint64_t sh, uint64_t sl){
    P.state = ((__uint128_t)sh<<64)|sl;
    pcg64(); pcg64();
}

int main(int argc, char **argv){
    if (argc<4){fprintf(stderr,"usage: %s N sh sl\n",argv[0]); return 1;}
    uint64_t N = strtoull(argv[1],0,0);
    uint64_t sh = strtoull(argv[2],0,0);
    uint64_t sl = strtoull(argv[3],0,0);
    gf_init();
    seed_rng(sh, sl);

    uint64_t hist[17]={0};

    for(uint64_t t=0;t<N;t++){
        uint64_t r0=pcg64(),r1=pcg64(),r2=pcg64(),r3=pcg64();
        uint8_t g[16],d[16];
        for(int i=0;i<8;i++){g[i]=(r0>>(8*i));g[8+i]=(r1>>(8*i));}
        for(int i=0;i<8;i++){d[i]=(r2>>(8*i));d[8+i]=(r3>>(8*i));}

        /* Vote-per-ω: solve s^2 * d^{-1} + s = g  ⟺  s^2 + s*d = g*d  (if d≠0)
           Let u=s/d: u^2+u = g/d. Solvable iff Tr(g/d)=0; roots u∈{HTR(g/d), HTR(g/d)+1}
           ⟹ s∈{u*d, (u+1)*d}.  If d=0: RHS=s, so match iff s=g (one s). */
        uint8_t votes[256]; memset(votes,0,256);
        for(int w=0;w<16;w++){
            uint8_t dw=d[w], gw=g[w];
            if(dw==0){
                votes[gw]++;  /* s=gw matches (incl. s=0 ignored later) */
            } else {
                uint8_t q = GMUL[gw][GINV[dw]]; /* g/d */
                if(TR[q]==0){
                    uint8_t u = HTR[q];
                    uint8_t s1 = GMUL[u][dw];
                    uint8_t s2 = GMUL[u^1][dw];
                    votes[s1]++; votes[s2]++;
                }
            }
        }
        int best=0;
        for(int s=1;s<256;s++) if(votes[s]>best) best=votes[s];
        hist[best]++;
    }

    printf("# bc_impl_b N=%llu seed=(0x%llX,0x%llX)\n",
           (unsigned long long)N,(unsigned long long)sh,(unsigned long long)sl);
    uint64_t ge13=0;
    for(int k=0;k<=16;k++){
        printf("max=%2d  %12llu\n",k,(unsigned long long)hist[k]);
        if(k>=13) ge13+=hist[k];
    }
    printf("# false-accepts @ tau=13: %llu  (rate = %.3e)\n",
           (unsigned long long)ge13,(double)ge13/(double)N);
    return 0;
}
