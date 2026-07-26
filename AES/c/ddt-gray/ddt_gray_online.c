// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// ddt_gray_online.c — REPRODUCE DDT-GRAY ONLINE (k_6-gray), r2
//
// Structure: per (pair, k_{-1}[diag]) instance, enumerate ~256 k_6[idiag]
// candidates as (~16 valid γ) × (2^4 DDT-branches). Per candidate, compute
// g_v = Δz_5[0]^{(v)} for v=1..255 via partial decrypt. Measure iSB/v.
//
// γ-step (cold): all 4 k_6 bytes change → 4 iSB/v.
// DDT-flip: k_6[pos_r] ↔ k_6[pos_r]⊕ΔC[pos_r] → 1 iSB/v.
//
// Claim (paper, DDT-Gray appendix): amortized 1.1875 iSB/v, c_on=2.1875, T_on=2^89.80.
//
// Honest tie-in: real 7R AES encryption, real ciphertexts, true k_6 among
// candidates, true-k_6 g_v == trace Δz_5[0].
//
// Build: gcc -O3 -march=native -o ddt_gray_online ddt_gray_online.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// primitives
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
static uint8_t GMUL[256][256];
static const uint8_t MC[4][4]={{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
static const uint8_t iMC[4][4]={{14,11,13,9},{9,14,11,13},{13,9,14,11},{11,13,9,14}};
static uint8_t DDTi[256][256],DDTi_sol[256][256]; // iSB DDT

static uint8_t gf_raw(uint8_t a,uint8_t b){uint8_t p=0;for(int i=0;i<8;i++){if(b&1)p^=a;uint8_t h=a&0x80;a<<=1;if(h)a^=0x1b;b>>=1;}return p;}
static void init_tables(void){
  for(int i=0;i<256;i++)iSB[SB[i]]=i;
  for(int a=0;a<256;a++)for(int b=0;b<256;b++)GMUL[a][b]=gf_raw(a,b);
  memset(DDTi,0,sizeof(DDTi));
  for(int d=0;d<256;d++)for(int x=0;x<256;x++){uint8_t o=iSB[x]^iSB[x^d];if(!DDTi[d][o])DDTi_sol[d][o]=x;DDTi[d][o]++;}
}
#define gm(a,b) GMUL[(uint8_t)(a)][(uint8_t)(b)]
static uint64_t rng_s;
static uint64_t xrand(void){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}
static uint8_t rb(void){return(uint8_t)(xrand()>>33);}

// byte-state AES
static void subB(uint8_t s[4][4]){for(int r=0;r<4;r++)for(int c=0;c<4;c++)s[r][c]=SB[s[r][c]];}
static void shrR(uint8_t s[4][4]){uint8_t t[4][4];for(int r=0;r<4;r++)for(int c=0;c<4;c++)t[r][c]=s[r][(c+r)&3];memcpy(s,t,16);}
static void mixC(uint8_t s[4][4]){uint8_t t[4][4];for(int c=0;c<4;c++)for(int r=0;r<4;r++){uint8_t v=0;for(int k=0;k<4;k++)v^=gm(MC[r][k],s[k][c]);t[r][c]=v;}memcpy(s,t,16);}
static void ark(uint8_t s[4][4],const uint8_t rk[4][4]){for(int r=0;r<4;r++)for(int c=0;c<4;c++)s[r][c]^=rk[r][c];}

static void aes_ks(const uint8_t key[16],uint8_t rk[8][4][4]){
  static const uint8_t Rc[10]={1,2,4,8,16,32,64,128,0x1b,0x36};
  for(int c=0;c<4;c++)for(int r=0;r<4;r++)rk[0][r][c]=key[4*c+r];
  for(int i=1;i<8;i++){
    uint8_t t[4]={SB[rk[i-1][1][3]],SB[rk[i-1][2][3]],SB[rk[i-1][3][3]],SB[rk[i-1][0][3]]};
    t[0]^=Rc[i-1];
    for(int r=0;r<4;r++)rk[i][r][0]=rk[i-1][r][0]^t[r];
    for(int c=1;c<4;c++)for(int r=0;r<4;r++)rk[i][r][c]=rk[i-1][r][c]^rk[i][r][c-1];
  }
}
// 7-round encrypt with trace of z_5[0]
static void enc7_trace(const uint8_t P[4][4],const uint8_t rk[8][4][4],
                       uint8_t C[4][4],uint8_t*z5_0){
  uint8_t s[4][4];memcpy(s,P,16);
  ark(s,rk[0]);
  for(int i=1;i<=5;i++){subB(s);shrR(s);mixC(s);ark(s,rk[i]);}
  // round 5 done; s = x_5 after ARK? No: after ark(rk[5]) s=w_4^k_4=x_5... wait.
  // Convention: x_0=P^rk[0]. round i (i=1..6): x_{i-1}→y→z→w; x_i=w^rk[i].
  // After i=5 loop iter: s=x_5. Need z_5: continue SB,SR.
  subB(s);shrR(s); // s=z_5
  *z5_0=s[0][0];
  mixC(s);ark(s,rk[6]); // s=x_6
  subB(s);shrR(s);ark(s,rk[7]); // round 7 (final, no MC): s=C
  memcpy(C,s,16);
}

// ─── iSB instrumentation ─────────────────────────────────────────────
static uint64_t g_isb=0;
static inline uint8_t iSBi(uint8_t x){g_isb++;return iSB[x];}

// ─── online state ────────────────────────────────────────────────────
static const int POS[4]={0,13,10,7}; // byte indices of y_6[r,0] after SR
// but state is [r][c]; byte b=4c+r → (r,c)=(b%4,b/4). pos_r=(r,(4-r)%4).
static inline uint8_t Cget(const uint8_t C[4][4],int r){return C[r][(4-r)&3];}

// Cache per v
static uint8_t pd[4][256];    // pd[r][v]=iSB(C_v[pos_r]^k_6[pos_r])
static uint8_t gdiff[256];    // g_v = Δz_5[0]^{(v)} vs ref (v=0)
static uint8_t k6cur[4];

// COLD: compute g_v for all v given k_6[4]. 4 iSB/v.
static void online_cold(const uint8_t CT[256][4][4],const uint8_t k6[4]){
  for(int r=0;r<4;r++)k6cur[r]=k6[r];
  for(int v=0;v<256;v++)for(int r=0;r<4;r++)pd[r][v]=iSBi(Cget(CT[v],r)^k6[r]);
  gdiff[0]=0;
  for(int v=1;v<256;v++){
    uint8_t g=0;for(int r=0;r<4;r++)g^=gm(iMC[0][r],pd[r][v]^pd[r][0]);gdiff[v]=g;
  }
}

// DDT-flip on byte r: k_6[r] ^= dC[r]. 1 iSB/v.
static void online_flip(const uint8_t CT[256][4][4],int r,uint8_t dCr){
  k6cur[r]^=dCr;
  uint8_t mcr=iMC[0][r];
  uint8_t new_ref=iSBi(Cget(CT[0],r)^k6cur[r]);
  uint8_t dref=new_ref^pd[r][0];pd[r][0]=new_ref;
  for(int v=1;v<256;v++){
    uint8_t np=iSBi(Cget(CT[v],r)^k6cur[r]);
    gdiff[v]^=gm(mcr,(np^pd[r][v])^dref);
    pd[r][v]=np;
  }
}

// cold reference (uninstrumented)
static void online_cold_ref(const uint8_t CT[256][4][4],const uint8_t k6[4],uint8_t g[256]){
  uint8_t p0[4];for(int r=0;r<4;r++)p0[r]=iSB[Cget(CT[0],r)^k6[r]];
  g[0]=0;
  for(int v=1;v<256;v++){uint8_t x=0;for(int r=0;r<4;r++)x^=gm(iMC[0][r],iSB[Cget(CT[v],r)^k6[r]]^p0[r]);g[v]=x;}
}

int main(int argc,char**argv){
  uint64_t seed=(argc>1)?strtoull(argv[1],0,0):0xC0FFEE;
  int nruns=(argc>2)?atoi(argv[2]):200;
  rng_s=seed?seed:1;init_tables();

  printf("=== DDT-GRAY r2 ONLINE (k_6-gray), seed=0x%lx, runs=%d ===\n\n",seed,nruns);

  uint64_t tot_cand=0,tot_cold=0,tot_flip=0,tot_isb=0;
  uint64_t ver_ok=0,ver_fail=0;
  int honest_ok=0, true_in_cands=0;
  double sum_ngamma=0,sum_ncand=0;

  for(int run=0;run<nruns;run++){
    // ── setup: random key, random δ-set, encrypt ───────────────────
    uint8_t key[16],rk[8][4][4];
    for(int i=0;i<16;i++)key[i]=rb();
    aes_ks(key,rk);
    // δ-set: reference x_1 state random; vary x_1[0,0] over 256 values.
    // P = iARK_{rk0}(iMC(iSR(iSB(x_1))))? We don't need P; we start from x_1.
    // Simpler: random y_1 ref (post-SB round 1), vary y_1[0,0] by ω — same as
    // offline convention. But online uses v=Δx_1[0]; here we bypass and just
    // iterate the δ-set directly from y_1 (the per-v vs per-ω distinction is
    // about INDEXING; for cost measurement the set is the same 256 elements).
    uint8_t y1[4][4];for(int r=0;r<4;r++)for(int c=0;c<4;c++)y1[r][c]=rb();
    uint8_t CT[256][4][4],z5true[256];
    for(int w=0;w<256;w++){
      uint8_t s[4][4];memcpy(s,y1,16);s[0][0]^=w;
      // forward from y_1: SR MC ARK(rk1) then rounds 2..7
      shrR(s);mixC(s);ark(s,rk[1]);
      for(int i=2;i<=5;i++){subB(s);shrR(s);mixC(s);ark(s,rk[i]);}
      subB(s);shrR(s);z5true[w]=s[0][0];mixC(s);ark(s,rk[6]);
      subB(s);shrR(s);ark(s,rk[7]);memcpy(CT[w],s,16);
    }
    // True k_6 = rk[7] at positions pos_r
    uint8_t k6true[4];for(int r=0;r<4;r++)k6true[r]=rk[7][r][(4-r)&3];

    // ── HONEST: verify cold(true k_6) g_v == trace Δz_5[0] ──────────
    {
      uint8_t gref[256];online_cold_ref(CT,k6true,gref);
      int ok=1;for(int v=1;v<256;v++)if(gref[v]!=(z5true[v]^z5true[0])){ok=0;break;}
      honest_ok+=ok;
    }

    // ── pair ΔC: synthesize from a real right-pair-like companion ───
    // We don't have a real 2^{105}-data right pair; use the TRUE γ=Δy_5[0] of
    // a random companion in the δ-set (w=1) propagated correctly. Actually we
    // just need ΔC[pos_r]≠0 for cost-measurement. Use random nonzero ΔC.
    uint8_t dC[4];for(int r=0;r<4;r++)dC[r]=rb()|1;

    // Check: is true k_6 a candidate for SOME γ? (it should be iff ΔC is from
    // a real pair; with random ΔC, generically no — so skip this check unless
    // we synthesize ΔC properly.)
    // Proper synthesis: pick random γ_true, compute Δx_6[r]=MC[r][0]·γ_true,
    // set ΔC[pos_r]=SB(x_6[r,0])^SB(x_6[r,0]^Δx_6[r]) where x_6[r,0]=iSB(y_6true[r,0]).
    // y_6true[r,0]=C_ref[pos_r]^k6true[r].
    uint8_t gamma_true=rb()|1;
    for(int r=0;r<4;r++){
      uint8_t y6=Cget(CT[0],r)^k6true[r];
      uint8_t x6=iSB[y6];
      uint8_t dx6=gm(MC[r][0],gamma_true);
      dC[r]=SB[x6]^SB[x6^dx6];
    }

    // ── enumerate k_6 candidates: γ∈[1,255], keep if all 4 DDTi(dC[r],MC[r][0]·γ)>0 ─
    // DDT=4 handling: each byte r may have 1 or 2 solution-PAIRS; iterate
    // Cartesian product of pair-representatives (outer), 4-bit BRGC (inner).
    int ngamma=0;int true_found=0;int n_ddt4=0;
    uint64_t ncand_run=0;
    for(int gmm=1;gmm<256;gmm++){
      uint8_t sols[4][2];int nsol[4];int valid=1;
      for(int r=0;r<4;r++){
        uint8_t dx6=gm(MC[r][0],gmm);
        uint8_t n=DDTi[dC[r]][dx6];
        if(!n){valid=0;break;}
        // collect pair-representatives: scan for solutions, group into {a,a^dC[r]} pairs
        nsol[r]=0;
        for(int a=0;a<256;a++){
          if((iSB[a]^iSB[a^dC[r]])==dx6 && a<(a^dC[r])){
            sols[r][nsol[r]++]=a;if(nsol[r]==2)break;
          }
        }
        if(nsol[r]==2)n_ddt4++;
      }
      if(!valid)continue;
      ngamma++;
      int nprod=nsol[0]*nsol[1]*nsol[2]*nsol[3];
      for(int pp=0;pp<nprod;pp++){
        int idx[4]={pp%nsol[0],(pp/nsol[0])%nsol[1],(pp/nsol[0]/nsol[1])%nsol[2],(pp/nsol[0]/nsol[1]/nsol[2])%nsol[3]};
        uint8_t k6b[4];for(int r=0;r<4;r++)k6b[r]=Cget(CT[0],r)^sols[r][idx[r]];

        // 4-bit BRGC over DDT-flips
        uint64_t isb0=g_isb;
        online_cold(CT,k6b);tot_cold++;
        {uint8_t gref[256];online_cold_ref(CT,k6cur,gref);
         int ok=1;for(int v=1;v<256;v++)if(gdiff[v]!=gref[v]){ok=0;break;}
         if(ok)ver_ok++;else ver_fail++;}
        {int eq=1;for(int r=0;r<4;r++)if(k6cur[r]!=k6true[r])eq=0;if(eq)true_found=1;}
        tot_cand++;ncand_run++;
        for(int i=1;i<16;i++){
          int j=__builtin_ctz(i);
          online_flip(CT,j,dC[j]);tot_flip++;tot_cand++;ncand_run++;
          {uint8_t gref[256];online_cold_ref(CT,k6cur,gref);
           int ok=1;for(int v=1;v<256;v++)if(gdiff[v]!=gref[v]){ok=0;break;}
           if(ok)ver_ok++;else ver_fail++;}
          {int eq=1;for(int r=0;r<4;r++)if(k6cur[r]!=k6true[r])eq=0;if(eq)true_found=1;}
        }
        tot_isb += g_isb-isb0;
      }
    }
    sum_ngamma+=ngamma;sum_ncand+=ncand_run;
    true_in_cands+=true_found;
  }

  double avg_ngamma=sum_ngamma/nruns;
  double avg_ncand=sum_ncand/nruns;
  double isb_per_cand=(double)tot_isb/tot_cand;
  double isb_per_v=isb_per_cand/255.0;
  // Actually cold does 256×4 iSB (includes ref v=0). Let me recompute honestly:
  // cold: 256×4=1024 iSB → 1024/255 per v. flip: 1+255×1=256 iSB → 256/255 per v.
  // amortized over 16: (1024+15×256)/16/255 = (1024+3840)/16/255 = 4864/4080 = 1.192
  double c_on=isb_per_v+1.0;
  double T_on=88+log2(255.0*c_on/160.0);

  printf(" HONEST: cold(true k_6) g_v == trace Δz_5[0]:  %d / %d %s\n",
         honest_ok,nruns,honest_ok==nruns?"✓":"✗");
  printf(" HONEST: true k_6 ∈ enumerated candidates:     %d / %d %s\n",
         true_in_cands,nruns,true_in_cands==nruns?"✓":"✗");
  printf(" Verification (inc==cold):                     %lu / %lu  %s\n",
         ver_ok,ver_ok+ver_fail,ver_fail==0?"✓":"✗ FAIL");
  printf("\n");
  printf(" avg #valid-γ per run:     %.2f  (expect ~255/16≈15.9)\n",avg_ngamma);
  printf(" avg #k_6-cand per run:    %.2f  (expect ~256)\n",avg_ncand);
  printf(" total candidates:         %lu  (cold=%lu flip=%lu)\n",tot_cand,tot_cold,tot_flip);
  printf(" total iSB:                %lu\n",tot_isb);
  printf(" ── iSB / candidate:       %.3f\n",isb_per_cand);
  printf(" ── iSB / (cand,v):        %.5f  [claim: 1.1875]\n",isb_per_v);
  printf(" ── c_on:                  %.4f  [claim: 2.1875]\n",c_on);
  printf(" ── T_on @ N=2^88:         2^%.3f  [paper claim: 2^89.80]\n",T_on);

  return (honest_ok==nruns && true_in_cands==nruns && ver_fail==0)?0:1;
}
