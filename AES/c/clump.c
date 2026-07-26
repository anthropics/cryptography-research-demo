// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// clump.c — high-throughput clumping test for the I_{m,n} fingerprint.
//
// Build:  gcc -O3 -march=native -funroll-loops -o clump clump.c
// Run:    ./clump [N] [seed]        (default N=3e8)
//
// Per sample: draw a random 255-multiset over GF(256), compute the 12-byte
// I_{m,n} fingerprint (with the 15-exponent fallback), pack the first 8 bytes
// big-endian into a uint64 key.  After N samples, sort and report collision
// counts / H2 estimates for k=5,6,7-byte prefixes.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// -------------------------------------------------------------- GF(2^8) ---

static uint8_t EXP[512], LOG[256], SQ[256];
static uint8_t MUL[256][256];

static void gf_init(void) {
    unsigned x = 1;
    for (int i = 0; i < 255; i++) {
        EXP[i] = (uint8_t)x;
        LOG[x] = (uint8_t)i;
        x ^= x << 1;                 // *3 (primitive)
        if (x & 0x100) x ^= 0x11B;
    }
    for (int i = 255; i < 512; i++) EXP[i] = EXP[i - 255];
    LOG[0] = 0;  // unused sentinel
    for (int a = 0; a < 256; a++) {
        SQ[a] = a ? EXP[(2 * LOG[a]) % 255] : 0;
        for (int b = 0; b < 256; b++)
            MUL[a][b] = (a && b) ? EXP[(LOG[a] + LOG[b]) % 255] : 0;
    }
}

// --------------------------------------------------- exponents & tables ---

static const int ALL_EXP[15] =
    {7,11,13,19,23,29,31,37,43,47,53,59,61,91,127};

// POW_ODD[v] = (v^1, v^3, ..., v^127) packed into 8 uint64 (64 bytes).
static uint64_t POW_ODD[256][8] __attribute__((aligned(64)));

// submask pairs {k, m-k} with 0<k<m-k, k⊂m, flattened per exponent
static int SP_off[16];               // SP_off[j]..SP_off[j+1] is exponent j's slice
static uint8_t SP_k[4096], SP_mk[4096];

static void tables_init(void) {
    for (int v = 0; v < 256; v++) {
        uint8_t *row = (uint8_t *)POW_ODD[v];
        for (int i = 0; i < 64; i++) {
            int r = 2*i + 1;
            row[i] = v ? EXP[((unsigned)LOG[v] * r) % 255] : 0;
        }
    }
    int pos = 0;
    for (int j = 0; j < 15; j++) {
        SP_off[j] = pos;
        int m = ALL_EXP[j];
        for (int k = 1; k < m; k++)
            if ((k & m) == k && k < m - k) {
                SP_k[pos] = (uint8_t)k;
                SP_mk[pos] = (uint8_t)(m - k);
                pos++;
            }
    }
    SP_off[15] = pos;
}

// ------------------------------------------------------------ xoshiro256**

static uint64_t s[4];
static inline uint64_t rotl(uint64_t x, int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t rng_u64(void){
    uint64_t r = rotl(s[1]*5,7)*9;
    uint64_t t = s[1]<<17;
    s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=rotl(s[3],45);
    return r;
}
static void rng_seed(uint64_t seed){
    // splitmix64 to fill state
    for(int i=0;i<4;i++){
        seed += 0x9E3779B97F4A7C15ULL;
        uint64_t z=seed;
        z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;
        z=(z^(z>>27))*0x94D049BB133111EBULL;
        s[i]=z^(z>>31);
    }
}

// -------------------------------------------------------- fingerprint ---

// returns 0 on sentinel (<13 nonzero P_m), else the packed 8-byte key
static inline uint64_t fp_one(void) {
    // 1. draw 256 random bytes (use 255)
    uint64_t rb[32];
    for (int i = 0; i < 32; i++) rb[i] = rng_u64();
    const uint8_t *g = (const uint8_t *)rb;

    // 2. p_odd[0..63] = XOR_w POW_ODD[g[w]]
    uint64_t p_odd[8] __attribute__((aligned(64))) = {0};
    for (int w = 0; w < 255; w++) {
        const uint64_t *row = POW_ODD[g[w]];
        p_odd[0]^=row[0]; p_odd[1]^=row[1]; p_odd[2]^=row[2]; p_odd[3]^=row[3];
        p_odd[4]^=row[4]; p_odd[5]^=row[5]; p_odd[6]^=row[6]; p_odd[7]^=row[7];
    }
    const uint8_t *po = (const uint8_t *)p_odd;      // po[i] = p_{2i+1}

    // 3. expand to p[1..127]
    uint8_t p[128];
    for (int i = 0; i < 64; i++) p[2*i+1] = po[i];
    for (int r = 1; r < 64; r++) p[2*r] = SQ[p[r]];

    // 4. P_m for the 15 exponents
    uint8_t P[15];
    for (int j = 0; j < 15; j++) {
        int m = ALL_EXP[j];
        uint8_t acc = p[m];
        for (int t = SP_off[j]; t < SP_off[j+1]; t++)
            acc ^= MUL[p[SP_k[t]]][p[SP_mk[t]]];
        P[j] = acc;
    }

    // 5. first 13 nonzero (fallback)
    int idx[13], c = 0;
    for (int j = 0; j < 15 && c < 13; j++)
        if (P[j]) idx[c++] = j;
    if (c < 13) return 0;                             // sentinel

    // 6. I_j = EXP[(log Pm0 * n_j - log Pn_j * m0) mod 255], j=0..7
    int m0 = ALL_EXP[idx[0]];
    int lPm0 = LOG[P[idx[0]]];
    uint64_t key = 0;
    for (int jj = 0; jj < 8; jj++) {
        int n  = ALL_EXP[idx[jj+1]];
        int lPn = LOG[P[idx[jj+1]]];
        // (a*n - b*m0) mod 255, keep nonnegative
        int e = ((lPm0 * n - lPn * m0) % 255 + 255) % 255;
        key = (key << 8) | EXP[e];
    }
    return key;
}

// --------------------------------------------------- sort & collisions ---

static int cmp_u64(const void *a, const void *b){
    uint64_t x=*(const uint64_t*)a, y=*(const uint64_t*)b;
    return (x>y)-(x<y);
}

static uint64_t count_prefix_collisions(const uint64_t *k, size_t M, int bytes){
    int sh = 8*(8-bytes);
    uint64_t C = 0, run = 1;
    uint64_t prev = k[0] >> sh;
    for (size_t i = 1; i < M; i++) {
        uint64_t v = k[i] >> sh;
        if (v == prev) { run++; }
        else { C += run*(run-1)/2; run = 1; prev = v; }
    }
    C += run*(run-1)/2;
    return C;
}

// ------------------------------------------------------------------- main

int main(int argc, char **argv) {
    size_t N = (argc>1)? strtoull(argv[1],0,10) : 300000000ULL;
    uint64_t seed = (argc>2)? strtoull(argv[2],0,10) : 0xDEADBEEFULL;

    gf_init(); tables_init(); rng_seed(seed);

    // quick self-test: GF sanity
    if (MUL[0x57][0x83] != 0xC1) { fprintf(stderr,"GF self-test FAIL\n"); return 1; }

    // deterministic cross-check vector: g[w] = w % 100
    // expected (from validated Python): 31 28 6f aa c1 f9 99 79 a9 82 c6 0b ba fb b6
    {
        static const uint8_t expect[15]=
            {0x31,0x28,0x6f,0xaa,0xc1,0xf9,0x99,0x79,0xa9,0x82,0xc6,0x0b,0xba,0xfb,0xb6};
        uint64_t p_odd[8]={0};
        for(int w=0;w<255;w++){
            const uint64_t *row=POW_ODD[w%100];
            for(int q=0;q<8;q++) p_odd[q]^=row[q];
        }
        const uint8_t *po=(const uint8_t*)p_odd;
        uint8_t p[128];
        for(int i=0;i<64;i++)p[2*i+1]=po[i];
        for(int r=1;r<64;r++)p[2*r]=SQ[p[r]];
        fprintf(stderr,"[*] selftest P_m(w%%100): ");
        int bad=0;
        for(int j=0;j<15;j++){
            int m=ALL_EXP[j]; uint8_t acc=p[m];
            for(int t=SP_off[j];t<SP_off[j+1];t++) acc^=MUL[p[SP_k[t]]][p[SP_mk[t]]];
            fprintf(stderr,"%02x ",acc);
            if(acc!=expect[j]) bad++;
        }
        fprintf(stderr," %s\n", bad?"FAIL":"OK");
        if(bad){fprintf(stderr,"    (%d mismatches vs Python reference)\n",bad);return 1;}
    }
    if (argc>1 && !strcmp(argv[1],"selftest")) return 0;

    fprintf(stderr,"[*] allocating %.2f GiB for keys ...\n", N*8.0/(1<<30));
    uint64_t *keys = (uint64_t*)malloc(N*sizeof(uint64_t));
    if(!keys){fprintf(stderr,"alloc fail\n");return 1;}

    fprintf(stderr,"[*] generating N=%zu fingerprints ...\n", N);
    size_t M=0, sentinels=0;
    clock_t t0=clock();
    for(size_t i=0;i<N;i++){
        uint64_t k = fp_one();
        if(k) keys[M++]=k; else sentinels++;
        if((i&((1<<24)-1))==0 && i){
            double dt=(double)(clock()-t0)/CLOCKS_PER_SEC;
            fprintf(stderr,"    %zu/%zu  %.2f M/s  ETA %.1f min\r",
                    i,N,i/dt/1e6,(N-i)/(i/dt)/60);
        }
    }
    double dt=(double)(clock()-t0)/CLOCKS_PER_SEC;
    fprintf(stderr,"\n[*] %zu fingerprints in %.1fs (%.2f M/s); sentinels=%zu\n",
            M,dt,M/dt/1e6,sentinels);

    fprintf(stderr,"[*] sorting ...\n");
    t0=clock();
    qsort(keys,M,sizeof(uint64_t),cmp_u64);
    fprintf(stderr,"    sort: %.1fs\n",(double)(clock()-t0)/CLOCKS_PER_SEC);

    double t_sort=(double)(clock()-t0)/CLOCKS_PER_SEC;

    FILE *f = fopen("clump_results.tsv","a");
    if (f && ftell(f)==0)
        fprintf(f,"N\tM\tseed\tk\tcollisions\tH2_hat\tideal\texpected_C\tgen_s\tsort_s\n");

    printf("\nPREFIX collision-entropy:\n");
    printf(" k        collisions C      H2_hat    ideal    expected C\n");
    for(int k=1;k<=8;k++){
        uint64_t C = count_prefix_collisions(keys,M,k);
        double ideal = k*log2(255.0);
        double expC = (double)M*(M-1)/2.0/pow(255.0,k);
        double H2 = C? log2((double)M)+log2((double)(M-1))-log2(2.0*C) : INFINITY;
        printf(" %d  %18llu   %8.3f  %7.3f   %.3e\n",
               k,(unsigned long long)C,H2,ideal,expC);
        if (f) fprintf(f,"%zu\t%zu\t%llu\t%d\t%llu\t%.6f\t%.6f\t%.6e\t%.2f\t%.2f\n",
                       N,M,(unsigned long long)seed,k,(unsigned long long)C,
                       H2,ideal,expC,dt,t_sort);
    }
    if (f) fclose(f);
    fprintf(stderr,"[*] results appended to clump_results.tsv\n");
    free(keys);
    return 0;
}
