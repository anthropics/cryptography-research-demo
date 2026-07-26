// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// wrongkey.c — exhaustive wrong-key randomization test for the I_{m,n}
// fingerprint (report.tex §3.1).
//
// Build:  gcc -O3 -march=native -funroll-loops -o wrongkey wrongkey.c
// Run:    ./wrongkey [N_keys] [seed]
//
// For each of N_keys random AES-128 keys:
//   - build a δ-set at x_1 (256 states, byte 0 varies), encrypt rounds 1..6;
//   - compute the correct offline fingerprint I_off (both parities);
//   - exhaustively enumerate ALL 2^32 guesses of k_6[0,13,10,7], compute the
//     online fingerprint, and record:
//        * ghost matches to I_off (want: exactly 1, the correct guess)
//        * Hamming-distance-to-I_off histogram
//        * exact prefix collision-entropy H2(k) for k=1,2,3
//
// Shares GF / fingerprint kernel with clump.c (duplicated here to keep the
// file self-contained).

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// =============================================== GF(2^8) and tables =======

static uint8_t EXP[512], LOG[256], SQ[256], MUL[256][256];
static uint8_t SBOX[256], SBOX_INV[256], A_LIN[256], A_LIN_INV[256], FIELD_INV[256];
static uint8_t GINV[256];                // GINV[x] = FIELD_INV[A_LIN_INV[x]]

static void gf_init(void) {
    unsigned x = 1;
    for (int i = 0; i < 255; i++) { EXP[i]=(uint8_t)x; LOG[x]=(uint8_t)i;
        x ^= x<<1; if (x&0x100) x^=0x11B; }
    for (int i=255;i<512;i++) EXP[i]=EXP[i-255];
    LOG[0]=0;
    for (int a=0;a<256;a++){ SQ[a]=a?EXP[(2*LOG[a])%255]:0;
        for(int b=0;b<256;b++) MUL[a][b]=(a&&b)?EXP[(LOG[a]+LOG[b])%255]:0; }
    // field inverse, S-box, affine L
    FIELD_INV[0]=0;
    for (int a=1;a<256;a++) FIELD_INV[a]=EXP[(255-LOG[a])%255];
    for (int a=0;a<256;a++){
        uint8_t b=FIELD_INV[a], r=b;
        for(int sh=1;sh<=4;sh++) r ^= (uint8_t)((b<<sh)|(b>>(8-sh)));
        SBOX[a]=r^0x63;
        // A_LIN[x]=A·x (linear part, no +0x63)
        uint8_t rb=a, rl=a;
        for(int sh=1;sh<=4;sh++) rl ^= (uint8_t)((rb<<sh)|(rb>>(8-sh)));
        A_LIN[a]=rl;
    }
    for (int a=0;a<256;a++){ SBOX_INV[SBOX[a]]=(uint8_t)a; A_LIN_INV[A_LIN[a]]=(uint8_t)a; }
    for (int a=0;a<256;a++) GINV[a]=FIELD_INV[A_LIN_INV[a]];
    // sanity
    if (SBOX[0x00]!=0x63 || SBOX[0x53]!=0xED) { fprintf(stderr,"SBOX FAIL\n"); exit(1); }
    if (MUL[0x57][0x83]!=0xC1) { fprintf(stderr,"GF FAIL\n"); exit(1); }
}

// =============================================== fingerprint kernel =======

static const int ALL_EXP[15]={7,11,13,19,23,29,31,37,43,47,53,59,61,91,127};
static uint64_t POW_ODD[256][8] __attribute__((aligned(64)));
static int SP_off[16]; static uint8_t SP_k[4096], SP_mk[4096];

static void tables_init(void){
    for(int v=0;v<256;v++){ uint8_t*row=(uint8_t*)POW_ODD[v];
        for(int i=0;i<64;i++){int r=2*i+1; row[i]=v?EXP[((unsigned)LOG[v]*r)%255]:0;} }
    int pos=0;
    for(int j=0;j<15;j++){ SP_off[j]=pos; int m=ALL_EXP[j];
        for(int k=1;k<m;k++) if((k&m)==k && k<m-k){SP_k[pos]=(uint8_t)k;SP_mk[pos]=(uint8_t)(m-k);pos++;} }
    SP_off[15]=pos;
}

static inline int fp_from_P(const uint8_t P[15], uint8_t I[12]){
    int idx[13],c=0; for(int j=0;j<15&&c<13;j++) if(P[j]) idx[c++]=j;
    if(c<13){memset(I,0,12);return 0;}
    int m0=ALL_EXP[idx[0]], lPm0=LOG[P[idx[0]]];
    for(int jj=0;jj<12;jj++){int n=ALL_EXP[idx[jj+1]], lPn=LOG[P[idx[jj+1]]];
        int e=((lPm0*n - lPn*m0)%255+255)%255; I[jj]=EXP[e]; }
    return 1;
}

// fingerprint of g[0..254] for BOTH parities at once (shared p_odd).
// Returns bitmask: bit0 = raw ok, bit1 = append-0 ok.
static inline int fp_both(const uint8_t *g, uint8_t I0[12], uint8_t I1[12]){
    uint64_t po[8] __attribute__((aligned(64)))={0};
    for(int w=0;w<255;w++){ const uint64_t*r=POW_ODD[g[w]];
        po[0]^=r[0];po[1]^=r[1];po[2]^=r[2];po[3]^=r[3];
        po[4]^=r[4];po[5]^=r[5];po[6]^=r[6];po[7]^=r[7]; }
    const uint8_t*pod=(const uint8_t*)po;
    uint8_t p[128]; for(int i=0;i<64;i++)p[2*i+1]=pod[i];
    for(int r=1;r<64;r++)p[2*r]=SQ[p[r]];
    uint8_t P0[15],P1[15];
    for(int j=0;j<15;j++){int m=ALL_EXP[j]; uint8_t acc=0;
        for(int t=SP_off[j];t<SP_off[j+1];t++) acc^=MUL[p[SP_k[t]]][p[SP_mk[t]]];
        P1[j]=acc;               // |V|=256 even: no lone p_m
        P0[j]=acc^p[m]; }        // |V|=255 odd:  plus lone p_m
    return fp_from_P(P0,I0) | (fp_from_P(P1,I1)<<1);
}

// single-parity helper (used for offline)
static inline int fp_of(const uint8_t *g, int len, uint8_t I[12]){
    uint8_t I0[12],I1[12]; int m=fp_both(g,I0,I1);
    // len==256 with g[255]==0 is equivalent to the even-parity output of g[0..254]
    if(len==255){memcpy(I,I0,12);return m&1;}
    memcpy(I,I1,12);return (m>>1)&1;
}

// =============================================== AES rounds 1..6 ==========

static inline uint8_t xt(uint8_t a){return (uint8_t)((a<<1)^((a>>7)*0x1B));}
static const uint8_t RCON[10]={0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1B,0x36};
static const int SR_IDX[16]={0,5,10,15, 4,9,14,3, 8,13,2,7, 12,1,6,11};

static void key_schedule(const uint8_t key[16], uint8_t rk[11][16]){
    memcpy(rk[0],key,16);
    for(int i=1;i<11;i++){
        uint8_t t[4]={rk[i-1][12],rk[i-1][13],rk[i-1][14],rk[i-1][15]};
        uint8_t r[4]={(uint8_t)(SBOX[t[1]]^RCON[i-1]),SBOX[t[2]],SBOX[t[3]],SBOX[t[0]]};
        for(int j=0;j<4;j++) rk[i][j]=rk[i-1][j]^r[j];
        for(int j=4;j<16;j++) rk[i][j]=rk[i-1][j]^rk[i][j-4];
    }
}

static void aes_round(uint8_t s[16], const uint8_t rk[16], int last){
    uint8_t t[16]; for(int i=0;i<16;i++) t[i]=SBOX[s[SR_IDX[i]]];
    if(!last){
        for(int c=0;c<4;c++){uint8_t a=t[4*c],b=t[4*c+1],d=t[4*c+2],e=t[4*c+3];
            uint8_t ab=a^b^d^e;
            s[4*c+0]=xt(a^b)^a^ab; s[4*c+1]=xt(b^d)^b^ab;
            s[4*c+2]=xt(d^e)^d^ab; s[4*c+3]=xt(e^a)^e^ab; }
    } else memcpy(s,t,16);
    for(int i=0;i<16;i++) s[i]^=rk[i];
}

// =============================================== xoshiro256** RNG =========

static uint64_t S_[4];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t rng64(void){uint64_t r=rotl(S_[1]*5,7)*9,t=S_[1]<<17;
    S_[2]^=S_[0];S_[3]^=S_[1];S_[1]^=S_[2];S_[0]^=S_[3];S_[2]^=t;S_[3]=rotl(S_[3],45);return r;}
static void rng_seed(uint64_t s){for(int i=0;i<4;i++){s+=0x9E3779B97F4A7C15ULL;
    uint64_t z=s;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;S_[i]=z^(z>>31);}}

// =============================================== main =====================

int main(int argc,char**argv){
    int N_keys=(argc>1)?atoi(argv[1]):1;
    uint64_t seed=(argc>2)?strtoull(argv[2],0,10):0xC0FFEE;
    gf_init(); tables_init(); rng_seed(seed);

    // per-byte MC^{-1} row-0 coefficients for v = 14,11,13,9 . x6[0..3]
    static const uint8_t COEF[4]={14,11,13,9};
    static const int K6POS[4]={0,13,10,7};   // z_6 positions -> y_6[0..3]

    for(int K=0;K<N_keys;K++){
        // --- random key, δ-set, encrypt rounds 1..6 ----------------------
        uint8_t key[16],base[16],rk[11][16];
        for(int i=0;i<16;i++){key[i]=(uint8_t)rng64(); base[i]=(uint8_t)rng64();}
        key_schedule(key,rk);
        uint8_t C[256][16], a[256];
        for(int w=0;w<256;w++){
            uint8_t s[16]; memcpy(s,base,16); s[0]=(uint8_t)w;
            for(int rnd=1;rnd<=6;rnd++){
                if(rnd==5) a[w]=s[0];               // x_5[0]
                aes_round(s,rk[rnd+1],rnd==6);
            }
            memcpy(C[w],s,16);
        }
        if(a[0]==0){fprintf(stderr,"  key %d: s=0, skip\n",K);continue;}

        // --- offline fingerprint (both parities) -------------------------
        uint8_t dinv[256];
        for(int w=1;w<256;w++) dinv[w-1]=FIELD_INV[a[0]^a[w]];
        uint8_t Ioff[2][12];
        fp_of(dinv,255,Ioff[0]);
        dinv[255]=0; fp_of(dinv,256,Ioff[1]);

        // --- cross-check: correct k_6 must match -------------------------
        // (also serves as harness self-test for this key)

        // precompute T[j][g][w] = COEF[j] * SBOX_INV[C[w][K6POS[j]] ^ g]
        static uint8_t T[4][256][256] __attribute__((aligned(64)));
        for(int j=0;j<4;j++)
            for(int g=0;g<256;g++)
                for(int w=0;w<256;w++)
                    T[j][g][w]=MUL[COEF[j]][SBOX_INV[C[w][K6POS[j]]^g]];

        uint8_t k6true[4]={rk[7][K6POS[0]],rk[7][K6POS[1]],rk[7][K6POS[2]],rk[7][K6POS[3]]};

        // --- exhaustive sweep over 2^32 guesses --------------------------
        uint64_t hd_hist[13]={0};        // Hamming-dist (bytes) to I_off[p*]
        uint64_t ghost=0;                // matches to I_off under either parity
        uint64_t sentinel=0;
        // prefix histograms on raw fingerprint
        static uint32_t H1[256], H2a[65536];
        static uint32_t *H3; if(!H3)H3=(uint32_t*)calloc(1u<<24,4);
        memset(H1,0,sizeof H1); memset(H2a,0,sizeof H2a);
        memset(H3,0,(size_t)4<<24);

        // determine correct parity once
        int p_ok=-1;
        {   uint8_t v[256],g[256],I[12];
            for(int w=0;w<256;w++)
                v[w]=T[0][k6true[0]][w]^T[1][k6true[1]][w]^T[2][k6true[2]][w]^T[3][k6true[3]][w];
            for(int w=1;w<256;w++) g[w-1]=GINV[v[0]^v[w]];
            fp_of(g,255,I);
            if(!memcmp(I,Ioff[0],12)) p_ok=0;
            else { g[255]=0; fp_of(g,256,I); if(!memcmp(I,Ioff[1],12)) p_ok=1; }
            if(p_ok<0){fprintf(stderr,"  key %d: correct k_6 FAILED to match — harness bug\n",K);return 1;}
        }

        fprintf(stderr,"  key %d: sweeping 2^32 k_6 guesses (parity=%d) ...\n",K,p_ok);
        clock_t t0=clock();
        uint8_t vpart012[256] __attribute__((aligned(64)));
        for(unsigned g0=0;g0<256;g0++){
         for(unsigned g1=0;g1<256;g1++){
          for(unsigned g2=0;g2<256;g2++){
            // partial v without the g3 term
            for(int w=0;w<256;w++)
                vpart012[w]=T[0][g0][w]^T[1][g1][w]^T[2][g2][w];
            for(unsigned g3=0;g3<256;g3++){
                uint8_t v0=vpart012[0]^T[3][g3][0];
                uint8_t g[256],I0[12],I1[12];
                for(int w=1;w<256;w++)
                    g[w-1]=GINV[v0 ^ vpart012[w] ^ T[3][g3][w]];
                int m=fp_both(g,I0,I1); int ok0=m&1, ok1=(m>>1)&1;
                // ghost: either parity matches its I_off
                if((ok0 && !memcmp(I0,Ioff[0],12)) ||
                   (ok1 && !memcmp(I1,Ioff[1],12))) ghost++;
                // histograms on the parity the correct key uses
                const uint8_t *Ip = p_ok? I1:I0;
                int okp = p_ok? ok1:ok0;
                if(!okp){sentinel++;continue;}
                int hd=0; for(int b=0;b<12;b++) hd+=(Ip[b]!=Ioff[p_ok][b]);
                hd_hist[hd]++;
                H1[Ip[0]]++;
                H2a[((unsigned)Ip[0]<<8)|Ip[1]]++;
                H3[((unsigned)Ip[0]<<16)|((unsigned)Ip[1]<<8)|Ip[2]]++;
            }
          }
          if(((g0<<8|g1)&0x3FF)==0){
            double dt=(double)(clock()-t0)/CLOCKS_PER_SEC, frac=((g0*256.0+g1+1)/65536.0);
            fprintf(stderr,"    %5.1f%%  %.2f Mguess/s  ETA %.1f min\r",
                    100*frac, frac*4294967296.0/dt/1e6, dt*(1/frac-1)/60);
          }
         }
        }
        double dt=(double)(clock()-t0)/CLOCKS_PER_SEC;
        fprintf(stderr,"\n    done in %.1fs (%.2f Mguess/s)\n",dt,4294967296.0/dt/1e6);

        // --- report ------------------------------------------------------
        printf("\n=== key %d ===\n",K);
        printf("ghost matches to I_off (incl. correct): %llu  (expect 1)\n",
               (unsigned long long)ghost);
        printf("sentinels: %llu\n",(unsigned long long)sentinel);
        printf("Hamming-distance-to-I_off histogram (bytes differing, 0..12):\n  ");
        double mean=0,tot=0;
        for(int h=0;h<=12;h++){printf("%llu ",(unsigned long long)hd_hist[h]);
            mean+=h*(double)hd_hist[h]; tot+=hd_hist[h];}
        printf("\n  mean = %.4f (uniform-random expects %.4f)\n",mean/tot,12*(1-1.0/255));
        // H2(k) from histograms
        uint64_t M=(uint64_t)tot;
        for(int kk=1;kk<=3;kk++){
            uint64_t Ck=0; uint32_t*Hh=(kk==1)?H1:(kk==2)?H2a:H3;
            size_t B=(kk==1)?256:(kk==2)?65536:(1u<<24);
            for(size_t i=0;i<B;i++){uint64_t v=Hh[i]; Ck+=v*(v-1)/2;}
            double h2=log2((double)M)+log2((double)(M-1))-log2(2.0*Ck);
            printf("H2(k=%d) = %.6f   ideal %.6f   (C=%llu)\n",
                   kk,h2,kk*log2(255.0),(unsigned long long)Ck);
        }
        // TSV for plotting
        FILE*f=fopen("wrongkey_results.tsv","a");
        if(f){if(ftell(f)==0)fprintf(f,"seed\tkey_idx\tghost\tsentinels\tM\t"
            "hd0\thd1\thd2\thd3\thd4\thd5\thd6\thd7\thd8\thd9\thd10\thd11\thd12\n");
            fprintf(f,"%llu\t%d\t%llu\t%llu\t%llu",(unsigned long long)seed,K,
                    (unsigned long long)ghost,(unsigned long long)sentinel,
                    (unsigned long long)M);
            for(int h=0;h<=12;h++)fprintf(f,"\t%llu",(unsigned long long)hd_hist[h]);
            fprintf(f,"\n");fclose(f);}
        fflush(stdout);
    }
    return 0;
}
