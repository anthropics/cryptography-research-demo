// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// ddt_gray_timing.c — wall-clock timing of pure DDT-GRAY inner loop
// (no instrumentation, no verification) to validate the "1 lookup = 1 SB"
// cost model against real hardware.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

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
static uint8_t iSB[256],GMUL[256][256],DDT[256][256],DDT_sol[256][256];
static const uint8_t MC[4][4]={{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
static const uint8_t iMC[4][4]={{14,11,13,9},{9,14,11,13},{13,9,14,11},{11,13,9,14}};
static uint8_t gf_raw(uint8_t a,uint8_t b){uint8_t p=0;for(int i=0;i<8;i++){if(b&1)p^=a;uint8_t h=a&0x80;a<<=1;if(h)a^=0x1b;b>>=1;}return p;}
#define gm(a,b) GMUL[(uint8_t)(a)][(uint8_t)(b)]
static uint64_t rng_s;
static uint64_t xrand(void){rng_s^=rng_s<<13;rng_s^=rng_s>>7;rng_s^=rng_s<<17;return rng_s;}
static uint8_t rb(void){return(uint8_t)(xrand()>>33);}

typedef struct{uint8_t x2[4],x3[4][4],x4[4];}Ref24;
typedef struct{uint8_t Din,x2[4],Dout,z4[4],Dx3[4][4],Dy3[4][4],Dx4[4],Dz4[4];}P10;

static int rebound(P10*P,Ref24*R){
  uint8_t Dy2[4];
  for(int r=0;r<4;r++){uint8_t d=gm(MC[r][0],P->Din);Dy2[r]=SB[P->x2[r]^d]^SB[P->x2[r]];if(!Dy2[r])return 0;}
  for(int r=0;r<4;r++)for(int c=0;c<4;c++){int rc=(4-c)&3;P->Dx3[r][c]=gm(MC[r][rc],Dy2[rc]);}
  for(int r=0;r<4;r++){P->Dz4[r]=gm(iMC[r][0],P->Dout);R->x4[r]=iSB[P->z4[r]];P->Dx4[r]=iSB[P->z4[r]^P->Dz4[r]]^R->x4[r];if(!P->Dx4[r])return 0;}
  for(int r=0;r<4;r++)for(int c=0;c<4;c++){int cc=((c-r)+4)&3;P->Dy3[r][c]=gm(iMC[r][cc],P->Dx4[cc]);}
  for(int r=0;r<4;r++)for(int c=0;c<4;c++){if(!DDT[P->Dx3[r][c]][P->Dy3[r][c]])return 0;R->x3[r][c]=DDT_sol[P->Dx3[r][c]][P->Dy3[r][c]];}
  for(int r=0;r<4;r++)R->x2[r]=P->x2[r];return 1;
}

static struct{uint8_t dx3[4][4],dy3[4][4],dx4[4],dy4[4],E;}C[256];
static uint8_t sb2c[4],sb3c[4][4],sb4c[4];
// xtime tables for ×2,×3 (only MC coeffs used in inner loop)
static uint8_t X2[256],X3[256];
static uint8_t POWINV[256][16];

static void cache_init(const Ref24*R){
  for(int r=0;r<4;r++)sb2c[r]=SB[R->x2[r]];
  for(int r=0;r<4;r++)for(int c=0;c<4;c++)sb3c[r][c]=SB[R->x3[r][c]];
  for(int r=0;r<4;r++)sb4c[r]=SB[R->x4[r]];
  for(int w=1;w<256;w++){
    uint8_t dy2[4];for(int r=0;r<4;r++)dy2[r]=SB[R->x2[r]^gm(MC[r][0],w)]^sb2c[r];
    for(int r=0;r<4;r++)for(int c=0;c<4;c++){int rc=(4-c)&3;C[w].dx3[r][c]=gm(MC[r][rc],dy2[rc]);C[w].dy3[r][c]=SB[R->x3[r][c]^C[w].dx3[r][c]]^sb3c[r][c];}
    for(int r=0;r<4;r++){uint8_t d=0;for(int rp=0;rp<4;rp++)d^=gm(MC[r][rp],C[w].dy3[rp][(r+rp)&3]);C[w].dx4[r]=d;C[w].dy4[r]=SB[R->x4[r]^d]^sb4c[r];}
    C[w].E=X2[C[w].dy4[0]]^X3[C[w].dy4[1]]^C[w].dy4[2]^C[w].dy4[3];
  }
}
// Inlined flips using xtime (×1,×2,×3 only)
static inline uint8_t mcmul(int a,uint8_t v){return a==1?v:(a==2?X2[v]:X3[v]);}
static inline void flip_z4(Ref24*R,int r,uint8_t d,uint8_t Sk[16]){
  R->x4[r]^=d;sb4c[r]=SB[R->x4[r]];int mc=MC[0][r];
  for(int w=1;w<256;w++){
    uint8_t n=SB[R->x4[r]^C[w].dx4[r]]^sb4c[r];
    C[w].E^=mcmul(mc,n^C[w].dy4[r]);C[w].dy4[r]=n;
    const uint8_t*pw=POWINV[C[w].E];
    for(int i=0;i<16;i++)Sk[i]^=pw[i]; // fingerprint accumulate (1 SIMD-XOR model)
  }
}
static inline void flip_x3(Ref24*R,int rp,int cp,uint8_t d,uint8_t Sk[16]){
  R->x3[rp][cp]^=d;sb3c[rp][cp]=SB[R->x3[rp][cp]];int rs=((cp-rp)+4)&3;int m3=MC[rs][rp],m4=MC[0][rs];
  for(int w=1;w<256;w++){
    uint8_t n3=SB[R->x3[rp][cp]^C[w].dx3[rp][cp]]^sb3c[rp][cp];
    C[w].dx4[rs]^=mcmul(m3,n3^C[w].dy3[rp][cp]);C[w].dy3[rp][cp]=n3;
    uint8_t n4=SB[R->x4[rs]^C[w].dx4[rs]]^sb4c[rs];
    C[w].E^=mcmul(m4,n4^C[w].dy4[rs]);C[w].dy4[rs]=n4;
    const uint8_t*pw=POWINV[C[w].E];
    for(int i=0;i<16;i++)Sk[i]^=pw[i];
  }
}

int main(int argc,char**argv){
  uint64_t seed=(argc>1)?strtoull(argv[1],0,0):0xC0FFEE;
  int nrep=(argc>2)?atoi(argv[2]):4;
  rng_s=seed?seed:1;
  for(int i=0;i<256;i++)iSB[SB[i]]=i;
  for(int a=0;a<256;a++)for(int b=0;b<256;b++)GMUL[a][b]=gf_raw(a,b);
  for(int i=0;i<256;i++){X2[i]=gf_raw(2,i);X3[i]=gf_raw(3,i);}
  memset(DDT,0,sizeof(DDT));
  for(int d=0;d<256;d++)for(int x=0;x<256;x++){uint8_t o=SB[x]^SB[x^d];if(!DDT[d][o])DDT_sol[d][o]=x;DDT[d][o]++;}
  for(int d=0;d<256;d++)for(int i=0;i<16;i++)POWINV[d][i]=rb(); // dummy content for timing

  printf("=== DDT-GRAY r2 wall-clock timing (pure inner loop, no verify) ===\n");

  P10 P;Ref24 R;
  for(;;){P.Din=rb()|1;P.Dout=rb()|1;for(int r=0;r<4;r++){P.x2[r]=rb();P.z4[r]=rb();}if(rebound(&P,&R))break;}

  uint64_t N=1ULL<<20;
  uint8_t Sk[16]={0};
  volatile uint64_t sink=0;

  for(int rep=0;rep<nrep;rep++){
    Ref24 Rr=R;
    cache_init(&Rr);
    struct timespec t0,t1;clock_gettime(CLOCK_MONOTONIC,&t0);
    for(uint64_t i=1;i<N;i++){
      int j=__builtin_ctzll(i);
      if(j<4)flip_z4(&Rr,j,P.Dx4[j],Sk);
      else{int p=j-4;flip_x3(&Rr,p&3,p>>2,P.Dx3[p&3][p>>2],Sk);}
    }
    clock_gettime(CLOCK_MONOTONIC,&t1);
    double dt=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
    double ns_per_entry=dt*1e9/N;
    double ns_per_lookup=ns_per_entry/(255*2.0664); // 2.0664 lookups/(entry,ω)
    sink+=Sk[0];
    printf(" rep %d: 2^20 entries in %.3fs → %.1f ns/entry → %.2f ns/(entry,ω) → %.3f ns/lookup-eq\n",
           rep,dt,ns_per_entry,ns_per_entry/255,ns_per_lookup);
  }
  // Reference: time a single 7R AES encryption (T-table, ~160 SB-eq)
  // for "enc-equivalent" calibration would go here; omitted (use 160 SB convention).
  printf(" (sink=%lu)\n",sink);
  return 0;
}
