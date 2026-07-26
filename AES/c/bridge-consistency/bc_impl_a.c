// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// bc_impl_a.c — Bridge-Consistency false-accept test, IMPLEMENTATION A (reference)
// GF(256) via log/alog; RNG = xoshiro256**; scalar loop.
// Protocol: draw (g[1..16], d[1..16]) uniform iid; compute max_{s=1..255} score(s);
//           where score(s) = #{ω : g[ω] == s^2 * d[ω]^{-1} XOR s}  (AES conv 0^{-1}=0)
// Output: histogram of max-score over N trials; #accepts at τ=13.
//
// build: gcc -O3 -march=native -o bc_a bc_impl_a.c
// run:   ./bc_a <Ntrials> <seed_hi> <seed_lo>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ---- GF(256) arithmetic, AES polynomial x^8+x^4+x^3+x+1, generator 0x03 ---- */
static uint8_t LOG[256], ALOG[510], INV[256], SQR[256];

static void gf_init(void) {
    uint16_t x = 1;
    for (int i = 0; i < 255; i++) {
        ALOG[i] = (uint8_t)x;
        ALOG[i+255] = (uint8_t)x;
        LOG[x] = (uint8_t)i;
        x ^= (x << 1) ^ ((x & 0x80) ? 0x11B : 0);
        x &= 0xFF;
    }
    LOG[0] = 0; /* sentinel; never used for nonzero-only paths */
    INV[0] = 0;
    for (int a = 1; a < 256; a++) INV[a] = ALOG[255 - LOG[a]];
    for (int a = 0; a < 256; a++) {
        SQR[a] = (a==0)?0:ALOG[(2*LOG[a]) % 255];
    }
}
static inline uint8_t gf_mul(uint8_t a, uint8_t b) {
    if (!a || !b) return 0;
    return ALOG[LOG[a] + LOG[b]];
}

/* ---- xoshiro256** (explicit-state) ---- */
static inline uint64_t rotl(uint64_t x, int k){ return (x<<k)|(x>>(64-k)); }
static inline uint64_t xoshiro(uint64_t *xs){
    uint64_t r = rotl(xs[1]*5,7)*9;
    uint64_t t = xs[1]<<17;
    xs[2]^=xs[0]; xs[3]^=xs[1]; xs[1]^=xs[2]; xs[0]^=xs[3]; xs[2]^=t; xs[3]=rotl(xs[3],45);
    return r;
}
static void seed_rng(uint64_t *xs, uint64_t hi, uint64_t lo){
    uint64_t z = hi;
    for(int i=0;i<4;i++){
        z += 0x9E3779B97F4A7C15ULL;
        uint64_t r = z;
        r = (r ^ (r>>30)) * 0xBF58476D1CE4E5B9ULL;
        r = (r ^ (r>>27)) * 0x94D049BB133111EBULL;
        r ^= r>>31;
        xs[i] = r ^ (i==0?lo:0);
    }
    for(int i=0;i<16;i++) xoshiro(xs);
}

int main(int argc, char **argv){
    if (argc < 4){ fprintf(stderr,"usage: %s N seed_hi seed_lo\n",argv[0]); return 1; }
    uint64_t N = strtoull(argv[1],0,0);
    uint64_t sh = strtoull(argv[2],0,0);
    uint64_t sl = strtoull(argv[3],0,0);
    gf_init();

    uint64_t hist[17] = {0};
    uint64_t sumscore_hist[33] = {0};

    #pragma omp parallel
    {
        uint64_t lh[17]={0}, ls[33]={0};
        uint8_t g[16], d[16], e[16];
        #ifdef _OPENMP
        int tid=omp_get_thread_num();
        #else
        int tid=0;
        #endif
        uint64_t xs[4];
        seed_rng(xs, sh ^ (0x9E3779B97F4A7C15ULL*(uint64_t)tid), sl);
        #pragma omp for schedule(static)
        for (uint64_t t = 0; t < N; t++){
        /* draw 32 random bytes */
        uint64_t r0 = xoshiro(xs), r1 = xoshiro(xs), r2 = xoshiro(xs), r3 = xoshiro(xs);
        for(int i=0;i<8;i++){ g[i]=(r0>>(8*i))&0xFF; g[8+i]=(r1>>(8*i))&0xFF; }
        for(int i=0;i<8;i++){ d[i]=(r2>>(8*i))&0xFF; d[8+i]=(r3>>(8*i))&0xFF; }
        for(int i=0;i<16;i++) e[i] = INV[d[i]];

        int best = 0;
        int sumscore = 0;
        for (int s = 1; s < 256; s++){
            uint8_t s2 = SQR[s];
            int sc = 0;
            for (int w = 0; w < 16; w++){
                uint8_t rhs = gf_mul(s2, e[w]) ^ (uint8_t)s;
                sc += (rhs == g[w]);
            }
            sumscore += sc;
            if (sc > best) best = sc;
        }
        lh[best]++;
        if (sumscore <= 32) ls[sumscore]++;
        }
        #pragma omp critical
        { for(int k=0;k<17;k++)hist[k]+=lh[k]; for(int k=0;k<33;k++)sumscore_hist[k]+=ls[k]; }
    }

    printf("# bc_impl_a N=%llu seed=(0x%llX,0x%llX)\n",
           (unsigned long long)N,(unsigned long long)sh,(unsigned long long)sl);
    printf("# max-score histogram:\n");
    uint64_t ge13 = 0;
    for(int k=0;k<=16;k++){
        printf("max=%2d  %12llu\n", k, (unsigned long long)hist[k]);
        if (k>=13) ge13 += hist[k];
    }
    printf("# false-accepts @ tau=13: %llu  (rate = %.3e)\n",
           (unsigned long long)ge13, (double)ge13/(double)N);
    printf("# sumscore (sanity; expect concentrated near 16):\n");
    for(int k=0;k<=32;k++) if(sumscore_hist[k])
        printf("sum=%2d  %12llu\n", k, (unsigned long long)sumscore_hist[k]);
    return 0;
}
