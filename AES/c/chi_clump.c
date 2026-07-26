// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// chi_clump.c — collision-entropy test for the χ*-canonicalization fingerprint
// (report.tex §3.3 / sec:chicanon), analogous to clump.c for I_{m,n}.
//
// Build:  gcc -O3 -march=native -funroll-loops -fopenmp -o chi_clump chi_clump.c -lm
// Run:    ./chi_clump [N] [seed] [threads]   (default N=3e8)
//
// Per sample: draw a random 255-multiset over GF(256), compute the 32-byte
// χ* fingerprint via the moving-frame canonical (α*=P_7^{-1/7}, β*=α*S_1,
// fallback P_11/P_13; n=255 odd), pack the first 8 bytes big-endian into a
// uint64 key.  After N samples, sort and report prefix collision counts /
// H2 estimates for k=1..8.  Ideal line is 8k (bytes of χ* ∈ {0..255}).
//
// AGL-invariance is self-tested on 10^4 random instances at startup.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// -------------------------------------------------------------- GF(2^8) ---
static uint8_t EXP[512], LOG[256], SQ[256];
static uint8_t MUL[256][256];

static void gf_init(void) {
    unsigned x = 1;
    for (int i = 0; i < 255; i++) { EXP[i]=(uint8_t)x; LOG[x]=(uint8_t)i;
        x ^= x<<1; if (x&0x100) x^=0x11B; }
    for (int i=255;i<512;i++) EXP[i]=EXP[i-255];
    LOG[0]=0;
    for (int a=0;a<256;a++){ SQ[a]=a?EXP[(2*LOG[a])%255]:0;
        for(int b=0;b<256;b++) MUL[a][b]=(a&&b)?EXP[(LOG[a]+LOG[b])%255]:0; }
}
static inline uint8_t gf_pow(uint8_t a,int e){
    return a?EXP[((unsigned)LOG[a]*(e%255))%255]:(e?0:1);
}

// packed POW_ODD[v] = (v^1,v^3,v^5,v^7) in 4 bytes; we only need r≤7 here
static uint32_t POW4[256];
static void tables_init(void){
    for(int v=0;v<256;v++){
        uint8_t r1=(uint8_t)v, r3=gf_pow(v,3), r5=gf_pow(v,5), r7=gf_pow(v,7);
        POW4[v] = r1 | (r3<<8) | (r5<<16) | ((uint32_t)r7<<24);
    }
}

// --------------------------------------------------------- χ* (n odd) ---
// Compute χ* of vals[0..254] (255 elements).  Output: 4×uint64 (256 bits).
// Returns 0 on sentinel (P_7=P_11=P_13=0).
static inline int chi_star(const uint8_t *vals, uint64_t out[4]){
    uint64_t chi[4]={0,0,0,0};
    uint32_t p4=0;                         // p1|p3|p5|p7
    for(int w=0;w<255;w++){
        uint8_t v=vals[w];
        chi[v>>6] ^= 1ULL<<(v&63);
        p4 ^= POW4[v];
    }
    uint8_t p1=p4, p3=p4>>8, p5=p4>>16, p7=p4>>24;
    uint8_t p2=SQ[p1], p4_=SQ[p2], p6=SQ[p3];
    // P_7 = p1 p6 ^ p2 p5 ^ p3 p4 ^ p7   (|V| odd)
    uint8_t P7 = MUL[p1][p6]^MUL[p2][p5]^MUL[p3][p4_]^p7;
    // α* from first nonzero of P_7, P_11, P_13; exponent = -inv(m,255) mod 255
    uint8_t alpha;
    if(P7){
        alpha = gf_pow(P7,182);            // 7·73 ≡ 1 (255) ⇒ -73 ≡ 182
    }else{
        // need p_r up to 13 for P_11,P_13; recompute via chi support
        // (rare: P_7=0 ~1/256).  Derive p_r from χ (odd-mult values).
        uint8_t p[14]={0};
        for(int d=0;d<256;d++) if(chi[d>>6]&(1ULL<<(d&63)))
            for(int r=1;r<14;r++) p[r]^=gf_pow((uint8_t)d,r);
        // Lucas for m=11 (1011b): pairs {1,10},{2,9},{3,8}
        uint8_t P11 = MUL[p[1]][p[10]]^MUL[p[2]][p[9]]^MUL[p[3]][p[8]]^p[11];
        if(P11){ alpha = gf_pow(P11,139); }     // inv(11,255)=116 ⇒ -116≡139
        else{
            // m=13 (1101b): pairs {1,12},{4,9},{5,8}
            uint8_t P13 = MUL[p[1]][p[12]]^MUL[p[4]][p[9]]^MUL[p[5]][p[8]]^p[13];
            if(!P13) return 0;
            alpha = gf_pow(P13,98);             // inv(13,255)=? 13·157=2041=8·255+1 ⇒ -157≡98
        }
    }
    uint8_t beta = MUL[alpha][p1];
    // permute: out_bit[α·d ⊕ β] = chi_bit[d]
    out[0]=out[1]=out[2]=out[3]=0;
    const uint8_t *row = MUL[alpha];
    for(int q=0;q<4;q++){
        uint64_t w=chi[q]; int base=q<<6;
        while(w){
            int d = base + __builtin_ctzll(w);
            w &= w-1;
            int t = row[d]^beta;
            out[t>>6] |= 1ULL<<(t&63);
        }
    }
    return 1;
}

// ------------------------------------------------------------ xoshiro256**
typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t rng_u64(rng_t *R){
    uint64_t r=rotl(R->s[1]*5,7)*9, t=R->s[1]<<17;
    R->s[2]^=R->s[0];R->s[3]^=R->s[1];R->s[1]^=R->s[2];R->s[0]^=R->s[3];
    R->s[2]^=t;R->s[3]=rotl(R->s[3],45); return r;
}
static void rng_seed(rng_t *R,uint64_t seed){
    for(int i=0;i<4;i++){ seed+=0x9E3779B97F4A7C15ULL; uint64_t z=seed;
        z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL;
        R->s[i]=z^(z>>31); }
}

// --------------------------------------------------- sort & collisions ---
static int cmp_u64(const void *a,const void *b){
    uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b; return (x>y)-(x<y);
}
static uint64_t prefix_coll(const uint64_t *k,size_t M,int bytes){
    int sh=8*(8-bytes); uint64_t C=0,run=1,prev=k[0]>>sh;
    for(size_t i=1;i<M;i++){ uint64_t v=k[i]>>sh;
        if(v==prev)run++; else{C+=run*(run-1)/2;run=1;prev=v;} }
    return C+run*(run-1)/2;
}

// ---------------------------------------------------------------- self-test
static void selftest(void){
    rng_t R; rng_seed(&R,1);
    int ok=0,sent=0;
    for(int t=0;t<10000;t++){
        uint8_t V[255],W[255]; uint64_t rb[32];
        for(int i=0;i<32;i++) rb[i]=rng_u64(&R);
        memcpy(V,rb,255);
        uint8_t a=1+(rng_u64(&R)%255), b=(uint8_t)rng_u64(&R);
        for(int i=0;i<255;i++) W[i]=MUL[a][V[i]]^b;
        uint64_t cV[4],cW[4];
        int rv=chi_star(V,cV), rw=chi_star(W,cW);
        if(!rv||!rw){sent++;continue;}
        if(!memcmp(cV,cW,32)) ok++;
        else { fprintf(stderr,"[!] AGL-invariance FAIL at t=%d\n",t); exit(1); }
    }
    fprintf(stderr,"[*] χ* AGL-invariance selftest: %d/%d OK (%d sentinels)\n",
            ok,10000,sent);
}

// ------------------------------------------------------------------- main
int main(int argc,char**argv){
    size_t N=(argc>1)?strtoull(argv[1],0,10):300000000ULL;
    uint64_t seed=(argc>2)?strtoull(argv[2],0,10):0xDEADBEEFULL;
    int nthr=(argc>3)?atoi(argv[3]):1;
    gf_init(); tables_init();

    selftest();
    if(argc>1 && !strcmp(argv[1],"selftest")) return 0;

    fprintf(stderr,"[*] allocating %.2f GiB for keys ...\n",N*8.0/(1<<30));
    uint64_t *keys=(uint64_t*)malloc(N*sizeof(uint64_t));
    if(!keys){fprintf(stderr,"alloc fail\n");return 1;}

#ifdef _OPENMP
    omp_set_num_threads(nthr);
#else
    (void)nthr;
#endif
    fprintf(stderr,"[*] generating N=%zu χ* fingerprints on %d thread(s)...\n",N,nthr);
    size_t sentinels=0;
    double t0=(double)clock()/CLOCKS_PER_SEC;
    #pragma omp parallel reduction(+:sentinels)
    {
        rng_t R;
#ifdef _OPENMP
        rng_seed(&R, seed + 0x1234567ULL*omp_get_thread_num());
        #pragma omp for schedule(static)
#else
        rng_seed(&R, seed);
#endif
        for(size_t i=0;i<N;i++){
            uint64_t rb[32]; for(int j=0;j<32;j++) rb[j]=rng_u64(&R);
            uint64_t out[4];
            if(!chi_star((const uint8_t*)rb,out)){ keys[i]=0; sentinels++; continue; }
            // first 8 bytes big-endian (bits 0..63 of χ*, byte 0 = bits 0..7 MSB)
            uint64_t k=0;
            for(int b=0;b<8;b++) k=(k<<8)|((out[0]>>(8*b))&0xFF);
            keys[i]=k;
        }
    }
    double dt=(double)clock()/CLOCKS_PER_SEC-t0;
    // compact out sentinels (key==0; genuine 0 key has prob 2^-64, ignore)
    size_t M=0; for(size_t i=0;i<N;i++) if(keys[i]) keys[M++]=keys[i];
    fprintf(stderr,"[*] %zu fps in %.1fs CPU (%.2f M/s/thr); sentinels=%zu\n",
            M,dt,N/dt/1e6,sentinels);

    fprintf(stderr,"[*] sorting ...\n");
    double ts0=(double)clock()/CLOCKS_PER_SEC;
    qsort(keys,M,sizeof(uint64_t),cmp_u64);
    double ts=(double)clock()/CLOCKS_PER_SEC-ts0;
    fprintf(stderr,"    sort: %.1fs\n",ts);

    FILE *f=fopen("chi_clump_results.tsv","a");
    if(f&&ftell(f)==0)
        fprintf(f,"N\tM\tseed\tk\tcollisions\tH2_hat\tideal\texpected_C\tgen_s\tsort_s\n");
    printf("\nPREFIX collision-entropy (χ*):\n");
    printf(" k        collisions C      H2_hat    ideal    expected C\n");
    for(int k=1;k<=8;k++){
        uint64_t C=prefix_coll(keys,M,k);
        double ideal=8.0*k;
        double expC=(double)M*(M-1)/2.0/pow(256.0,k);
        double H2=C?log2((double)M)+log2((double)(M-1))-log2(2.0*C):INFINITY;
        printf(" %d  %18llu   %8.4f  %7.3f   %.3e\n",
               k,(unsigned long long)C,H2,ideal,expC);
        if(f)fprintf(f,"%zu\t%zu\t%llu\t%d\t%llu\t%.6f\t%.6f\t%.6e\t%.2f\t%.2f\n",
                     N,M,(unsigned long long)seed,k,(unsigned long long)C,
                     H2,ideal,expC,dt,ts);
    }
    if(f)fclose(f);
    fprintf(stderr,"[*] appended chi_clump_results.tsv\n");
    free(keys);
    return 0;
}
