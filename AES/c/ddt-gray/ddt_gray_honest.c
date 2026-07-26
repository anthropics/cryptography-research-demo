// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// ddt_gray_honest.c — tie DDT-GRAY to the mobius-honest.0 pipeline:
//   (A) Verify the 24-ref cold_E formula against REAL AES encryption trace
//       (50 random keys, δ-set at y_1[0], assert E[ω]==Δx_5[0]^{(ω)} ∀ω).
//   (B) Distinctness: are all 2^20 gray-code entries DISTINCT E' multisets?
//   (C) XOR/gfmul accounting alongside SB for full cost-model table.
//
// Build: gcc -O3 -march=native -o ddt_gray_honest ddt_gray_honest.c -lm
// Run:   ./ddt_gray_honest [seed] [nkeys]

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// ─── primitives (same as ddt_gray_e2e.c) ─────────────────────────────
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
static const uint8_t MC[4][4]  = {{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
static const uint8_t iMC[4][4] = {{14,11,13,9},{9,14,11,13},{13,9,14,11},{11,13,9,14}};
static uint8_t DDT[256][256], DDT_sol[256][256];

static uint8_t gf_raw(uint8_t a,uint8_t b){uint8_t p=0;for(int i=0;i<8;i++){if(b&1)p^=a;uint8_t h=a&0x80;a<<=1;if(h)a^=0x1b;b>>=1;}return p;}
static void init_tables(void){
  for(int i=0;i<256;i++)iSB[SB[i]]=i;
  for(int a=0;a<256;a++)for(int b=0;b<256;b++)GMUL[a][b]=gf_raw(a,b);
  memset(DDT,0,sizeof(DDT));
  for(int d=0;d<256;d++)for(int x=0;x<256;x++){uint8_t o=SB[x]^SB[x^d];if(!DDT[d][o])DDT_sol[d][o]=x;DDT[d][o]++;}
}
#define gm(a,b) GMUL[(uint8_t)(a)][(uint8_t)(b)]
static uint64_t rng_s;
static uint64_t xrand(void){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}
static uint8_t rb(void){return(uint8_t)(xrand()>>33);}

// ─── byte-state AES round (for trace) ────────────────────────────────
// state[r][c]; sub, shiftrows, mixcols, ark
static void subB(uint8_t s[4][4]){for(int r=0;r<4;r++)for(int c=0;c<4;c++)s[r][c]=SB[s[r][c]];}
static void shrR(uint8_t s[4][4]){uint8_t t[4][4];for(int r=0;r<4;r++)for(int c=0;c<4;c++)t[r][c]=s[r][(c+r)&3];memcpy(s,t,16);}
static void mixC(uint8_t s[4][4]){uint8_t t[4][4];for(int c=0;c<4;c++)for(int r=0;r<4;r++){uint8_t v=0;for(int k=0;k<4;k++)v^=gm(MC[r][k],s[k][c]);t[r][c]=v;}memcpy(s,t,16);}
static void ark(uint8_t s[4][4],const uint8_t rk[4][4]){for(int r=0;r<4;r++)for(int c=0;c<4;c++)s[r][c]^=rk[r][c];}

// ─── (A) HONEST tie-in: cold_E formula vs real AES trace ────────────
typedef struct{uint8_t x2[4],x3[4][4],x4[4];}Ref24;

static void cold_E(const Ref24*R,uint8_t E[256]){
  uint8_t sb2[4],sb3[4][4],sb4[4];
  for(int r=0;r<4;r++)sb2[r]=SB[R->x2[r]];
  for(int r=0;r<4;r++)for(int c=0;c<4;c++)sb3[r][c]=SB[R->x3[r][c]];
  for(int r=0;r<4;r++)sb4[r]=SB[R->x4[r]];
  E[0]=0;
  for(int w=1;w<256;w++){
    uint8_t dy2[4];for(int r=0;r<4;r++)dy2[r]=SB[R->x2[r]^gm(MC[r][0],w)]^sb2[r];
    uint8_t dy3[4][4];
    for(int r=0;r<4;r++)for(int c=0;c<4;c++){int rc=(4-c)&3;dy3[r][c]=SB[R->x3[r][c]^gm(MC[r][rc],dy2[rc])]^sb3[r][c];}
    uint8_t e=0;
    for(int r=0;r<4;r++){uint8_t dx4=0;for(int rp=0;rp<4;rp++)dx4^=gm(MC[r][rp],dy3[rp][(r+rp)&3]);e^=gm(MC[0][r],SB[R->x4[r]^dx4]^sb4[r]);}
    E[w]=e;
  }
}

static int test_honest_trace(int nkeys){
  // For each key: random round keys rk[0..4] (independent-subkey model, as
  // honest_rebound.c — the formula is key-independent so KS doesn't matter).
  // Build δ-set by varying y_1[0] over all 256 values; extract x_2,x_3,x_4,x_5
  // from trace; assert cold_E(24-refs)==Δx_5[0] for all ω.
  int pass=0;
  for(int K=0;K<nkeys;K++){
    uint8_t rk[5][4][4];
    for(int i=0;i<5;i++)for(int r=0;r<4;r++)for(int c=0;c<4;c++)rk[i][r][c]=rb();
    // Reference state at y_1: random except y_1[0,0] will iterate.
    uint8_t y1ref[4][4];for(int r=0;r<4;r++)for(int c=0;c<4;c++)y1ref[r][c]=rb();
    // Trace all 256 δ-set members to x_5[0]; record refs at ω=0.
    Ref24 R; uint8_t x5_true[256];
    for(int w=0;w<256;w++){
      uint8_t s[4][4];memcpy(s,y1ref,16);s[0][0]=y1ref[0][0]^w; // Δy_1[0]=w
      // y_1→z_1→w_1→x_2
      shrR(s);mixC(s);ark(s,rk[1]);
      if(w==0)for(int r=0;r<4;r++)R.x2[r]=s[r][0];
      // x_2→y_2→z_2→w_2→x_3
      subB(s);shrR(s);mixC(s);ark(s,rk[2]);
      if(w==0)for(int r=0;r<4;r++)for(int c=0;c<4;c++)R.x3[r][c]=s[r][c];
      // x_3→…→x_4
      subB(s);shrR(s);mixC(s);ark(s,rk[3]);
      if(w==0)for(int r=0;r<4;r++)R.x4[r]=s[r][r];
      // x_4→…→x_5
      subB(s);shrR(s);mixC(s);ark(s,rk[4]);
      x5_true[w]=s[0][0];
    }
    uint8_t E[256];cold_E(&R,E);
    int ok=1;
    for(int w=1;w<256;w++)if(E[w]!=(x5_true[w]^x5_true[0])){ok=0;break;}
    pass+=ok;
  }
  return pass;
}

// ─── (B) distinctness of 2^20 entries ───────────────────────────────
// Reuse rebound + incremental from e2e; compute 64-bit hash of E[1..255]
// per entry, store in sorted array, count collisions.
typedef struct{uint8_t Din,x2[4],Dout,z4[4],Dx3[4][4],Dy3[4][4],Dx4[4],Dz4[4];}Params10;
static int rebound(Params10*P,Ref24*R){
  uint8_t Dy2[4];
  for(int r=0;r<4;r++){uint8_t d=gm(MC[r][0],P->Din);Dy2[r]=SB[P->x2[r]^d]^SB[P->x2[r]];if(!Dy2[r])return 0;}
  for(int r=0;r<4;r++)for(int c=0;c<4;c++){int rc=(4-c)&3;P->Dx3[r][c]=gm(MC[r][rc],Dy2[rc]);}
  for(int r=0;r<4;r++){P->Dz4[r]=gm(iMC[r][0],P->Dout);R->x4[r]=iSB[P->z4[r]];P->Dx4[r]=iSB[P->z4[r]^P->Dz4[r]]^R->x4[r];if(!P->Dx4[r])return 0;}
  for(int r=0;r<4;r++)for(int c=0;c<4;c++){int cc=((c-r)+4)&3;P->Dy3[r][c]=gm(iMC[r][cc],P->Dx4[cc]);}
  for(int r=0;r<4;r++)for(int c=0;c<4;c++){if(!DDT[P->Dx3[r][c]][P->Dy3[r][c]])return 0;R->x3[r][c]=DDT_sol[P->Dx3[r][c]][P->Dy3[r][c]];}
  for(int r=0;r<4;r++)R->x2[r]=P->x2[r];
  return 1;
}

// minimal cache + flips (no instrumentation; fast)
static struct{uint8_t dx3[4][4],dy3[4][4],dx4[4],dy4[4],E;}C[256];
static uint8_t sb2c[4],sb3c[4][4],sb4c[4];
static uint64_t g_sb=0,g_gm=0,g_xor=0;
#define SBc(x) (g_sb++,SB[(uint8_t)(x)])
#define GMc(a,b) (g_gm++,GMUL[(uint8_t)(a)][(uint8_t)(b)])
#define XRc(a,b) (g_xor++,(uint8_t)((a)^(b)))

static void cache_init(const Ref24*R){
  for(int r=0;r<4;r++)sb2c[r]=SBc(R->x2[r]);
  for(int r=0;r<4;r++)for(int c=0;c<4;c++)sb3c[r][c]=SBc(R->x3[r][c]);
  for(int r=0;r<4;r++)sb4c[r]=SBc(R->x4[r]);
  for(int w=1;w<256;w++){
    uint8_t dy2[4];for(int r=0;r<4;r++)dy2[r]=XRc(SBc(XRc(R->x2[r],GMc(MC[r][0],w))),sb2c[r]);
    for(int r=0;r<4;r++)for(int c=0;c<4;c++){int rc=(4-c)&3;C[w].dx3[r][c]=GMc(MC[r][rc],dy2[rc]);C[w].dy3[r][c]=XRc(SBc(XRc(R->x3[r][c],C[w].dx3[r][c])),sb3c[r][c]);}
    for(int r=0;r<4;r++){uint8_t d=0;for(int rp=0;rp<4;rp++)d=XRc(d,GMc(MC[r][rp],C[w].dy3[rp][(r+rp)&3]));C[w].dx4[r]=d;C[w].dy4[r]=XRc(SBc(XRc(R->x4[r],d)),sb4c[r]);}
    C[w].E=XRc(XRc(GMc(2,C[w].dy4[0]),GMc(3,C[w].dy4[1])),XRc(C[w].dy4[2],C[w].dy4[3]));
  }
}
static void flip_z4(Ref24*R,int r,uint8_t d){
  R->x4[r]^=d;sb4c[r]=SBc(R->x4[r]);uint8_t m=MC[0][r];
  for(int w=1;w<256;w++){uint8_t n=XRc(SBc(XRc(R->x4[r],C[w].dx4[r])),sb4c[r]);C[w].E=XRc(C[w].E,GMc(m,XRc(n,C[w].dy4[r])));C[w].dy4[r]=n;}
}
static void flip_x3(Ref24*R,int rp,int cp,uint8_t d){
  R->x3[rp][cp]^=d;sb3c[rp][cp]=SBc(R->x3[rp][cp]);int rs=((cp-rp)+4)&3;uint8_t m3=MC[rs][rp],m4=MC[0][rs];
  for(int w=1;w<256;w++){
    uint8_t n3=XRc(SBc(XRc(R->x3[rp][cp],C[w].dx3[rp][cp])),sb3c[rp][cp]);
    C[w].dx4[rs]=XRc(C[w].dx4[rs],GMc(m3,XRc(n3,C[w].dy3[rp][cp])));C[w].dy3[rp][cp]=n3;
    uint8_t n4=XRc(SBc(XRc(R->x4[rs],C[w].dx4[rs])),sb4c[rs]);
    C[w].E=XRc(C[w].E,GMc(m4,XRc(n4,C[w].dy4[rs])));C[w].dy4[rs]=n4;
  }
}

static uint64_t hashE(void){
  uint64_t h=0xcbf29ce484222325ULL;
  for(int w=1;w<256;w++){h^=C[w].E;h*=0x100000001b3ULL;}
  return h;
}
static int cmp64(const void*a,const void*b){
  uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return (x>y)-(x<y);
}

int main(int argc,char**argv){
  uint64_t seed=(argc>1)?strtoull(argv[1],0,0):0xC0FFEE;
  int nkeys=(argc>2)?atoi(argv[2]):50;
  rng_s=seed?seed:1;init_tables();

  printf("=== DDT-GRAY r2 HONEST tie-in, seed=0x%lx ===\n\n",seed);

  // (A) honest trace
  printf("(A) cold_E vs real-AES trace (%d random keys):\n",nkeys);
  int pass=test_honest_trace(nkeys);
  printf("    cold_E(24-refs) == Δx_5[0]^{(ω)} ∀ω:  %d / %d %s\n\n",
         pass,nkeys,pass==nkeys?"✓":"✗ FAIL");

  // (B) distinctness + (C) op-accounting over full 2^20
  printf("(B) 2^20 BRGC distinctness + (C) op accounting:\n");
  Params10 P;Ref24 R;int tries=0;
  for(;;){tries++;P.Din=rb()|1;P.Dout=rb()|1;for(int r=0;r<4;r++){P.x2[r]=rb();P.z4[r]=rb();}if(rebound(&P,&R))break;}
  printf("    base after %d tries\n",tries);

  uint64_t N=1ULL<<20;
  uint64_t*H=malloc(N*sizeof(uint64_t));
  g_sb=g_gm=g_xor=0;
  cache_init(&R);
  uint64_t sb0=g_sb,gm0=g_gm,xr0=g_xor;
  H[0]=hashE();
  for(uint64_t i=1;i<N;i++){
    int j=__builtin_ctzll(i);
    if(j<4)flip_z4(&R,j,P.Dx4[j]);
    else{int p=j-4;flip_x3(&R,p&3,p>>2,P.Dx3[p&3][p>>2]);}
    H[i]=hashE();
  }
  uint64_t sb1=g_sb-sb0,gm1=g_gm-gm0,xr1=g_xor-xr0;

  qsort(H,N,sizeof(uint64_t),cmp64);
  uint64_t coll=0;for(uint64_t i=1;i<N;i++)if(H[i]==H[i-1])coll++;
  free(H);

  double denom=(double)(N-1)*255.0;
  printf("    distinct E-hashes: %lu / %lu  (collisions=%lu)\n",N-coll,N,coll);
  printf("    per-(entry,ω) step-cost:\n");
  printf("      SB    : %.5f\n",sb1/denom);
  printf("      gfmul : %.5f\n",gm1/denom);
  printf("      XOR   : %.5f\n",xr1/denom);
  printf("      [POWINV: 1.00000 — added separately]\n");
  printf("    c_off (SB+POWINV, XOR=free model): %.4f → T_off=2^%.3f @N=2^{88}\n",
         sb1/denom+1.0, 88+log2(255.0*(sb1/denom+1.0)/160.0));
  printf("    c_off (SB+gfmul+POWINV): %.4f → T_off=2^%.3f\n",
         (sb1+gm1)/denom+1.0, 88+log2(255.0*((sb1+gm1)/denom+1.0)/160.0));

  return (pass==nkeys && coll==0)?0:1;
}
