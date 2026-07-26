// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// ddt_gray_e2e.c — REPRODUCE DDT-GRAY (r2): full 2^20 branch enumeration
// with instrumented SB counting and cold-vs-incremental verification.
//
// Approach (different from r1's micro-experiments):
//   1. Start from mobius-honest.0 pipeline: 10 params → rebound → 24 refs
//      → per-ω forward-prop → E[ω]=Δx_5[0]^{(ω)}.
//   2. Add DDT-GRAY: enumerate all 2^20 branches (4 z_4-flips innermost,
//      16 x_3-flips outer) via BRGC, incrementally updating E[ω].
//   3. Instrument: count every SB lookup in the incremental path.
//   4. Verify: for each branch (or dense sample), recompute E[ω] cold and
//      assert byte-equal.
//   5. Report: amortized SB/(entry,ω), does T_off=2^{89.72} hold?
//
// Build: gcc -O3 -march=native -o ddt_gray_e2e ddt_gray_e2e.c -lm
// Run:   ./ddt_gray_e2e [seed] [verify_mode]
//        verify_mode: 0=sample(fast), 1=full(all 2^20 entries cold-checked)

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// ─────────────────────────────────────────────────────────────────────
// GF(256) + AES primitives
// ─────────────────────────────────────────────────────────────────────
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
static uint8_t GINV[256];
// DDT[din][dout] = #{x : SB(x)^SB(x^din)=dout}; DDT_sol[din][dout] = smallest such x
static uint8_t DDT[256][256];
static uint8_t DDT_sol[256][256];

// MC matrix (row r, col c): circ(2,3,1,1)
static const uint8_t MC[4][4]  = {{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
static const uint8_t iMC[4][4] = {{14,11,13,9},{9,14,11,13},{13,9,14,11},{11,13,9,14}};

static uint8_t gf_raw(uint8_t a, uint8_t b){
  uint8_t p=0;
  for(int i=0;i<8;i++){ if(b&1)p^=a; uint8_t h=a&0x80; a<<=1; if(h)a^=0x1b; b>>=1; }
  return p;
}
static void init_tables(void){
  for(int i=0;i<256;i++) iSB[SB[i]]=i;
  for(int a=0;a<256;a++) for(int b=0;b<256;b++) GMUL[a][b]=gf_raw(a,b);
  GINV[0]=0;
  for(int a=1;a<256;a++) for(int b=1;b<256;b++) if(GMUL[a][b]==1){GINV[a]=b;break;}
  memset(DDT,0,sizeof(DDT));
  for(int din=0;din<256;din++) for(int x=0;x<256;x++){
    uint8_t dout=SB[x]^SB[x^din];
    if(DDT[din][dout]==0) DDT_sol[din][dout]=x;
    DDT[din][dout]++;
  }
}
#define gm(a,b) GMUL[(uint8_t)(a)][(uint8_t)(b)]

// ─────────────────────────────────────────────────────────────────────
// SB instrumentation (for incremental path ONLY)
// ─────────────────────────────────────────────────────────────────────
static uint64_t g_sb_count = 0;
static inline uint8_t SBi(uint8_t x){ g_sb_count++; return SB[x]; }

// ─────────────────────────────────────────────────────────────────────
// xorshift PRNG
// ─────────────────────────────────────────────────────────────────────
static uint64_t rng_s;
static uint64_t xrand(void){
  rng_s^=rng_s<<13; rng_s^=rng_s>>7; rng_s^=rng_s<<17; return rng_s;
}
static uint8_t rb(void){ return (uint8_t)(xrand()>>33); }

// ─────────────────────────────────────────────────────────────────────
// 24-ref state
// ─────────────────────────────────────────────────────────────────────
typedef struct {
  uint8_t x2[4];      // x_2[r,0], r=0..3
  uint8_t x3[4][4];   // x_3[r][c]
  uint8_t x4[4];      // x_4[r,r], r=0..3
} Ref24;

// Right-pair params (10 bytes) + derived DDT-pair flips
typedef struct {
  uint8_t Din;          // Δy_1[0] of right pair (= "Δz_1[0]")
  uint8_t x2[4];
  uint8_t Dout;         // Δw_4[0] of right pair
  uint8_t z4[4];
  // derived:
  uint8_t Dx3[4][4];    // right-pair Δx_3[r][c] (forward)
  uint8_t Dy3[4][4];    // right-pair Δy_3[r][c] (backward)
  uint8_t Dx4[4];       // right-pair Δx_4[r,r]
  uint8_t Dz4[4];       // right-pair Δz_4[r,0] = iMC[r][0]·Dout
} Params10;

// ─────────────────────────────────────────────────────────────────────
// COLD path: 24 refs → E[ω] for all ω. Ground truth (uninstrumented).
// ─────────────────────────────────────────────────────────────────────
static void cold_E(const Ref24 *R, uint8_t E[256]){
  uint8_t sb_x2[4], sb_x3[4][4], sb_x4[4];
  for(int r=0;r<4;r++) sb_x2[r]=SB[R->x2[r]];
  for(int r=0;r<4;r++) for(int c=0;c<4;c++) sb_x3[r][c]=SB[R->x3[r][c]];
  for(int r=0;r<4;r++) sb_x4[r]=SB[R->x4[r]];
  E[0]=0;
  for(int w=1;w<256;w++){
    uint8_t dy2[4];
    for(int r=0;r<4;r++) dy2[r]=SB[R->x2[r]^gm(MC[r][0],w)]^sb_x2[r];
    uint8_t dy3[4][4];
    for(int r=0;r<4;r++) for(int c=0;c<4;c++){
      int rc=(4-c)&3;
      uint8_t dx3=gm(MC[r][rc],dy2[rc]);
      dy3[r][c]=SB[R->x3[r][c]^dx3]^sb_x3[r][c];
    }
    uint8_t dy4[4];
    for(int r=0;r<4;r++){
      uint8_t dx4=0;
      for(int rp=0;rp<4;rp++) dx4^=gm(MC[r][rp],dy3[rp][(r+rp)&3]);
      dy4[r]=SB[R->x4[r]^dx4]^sb_x4[r];
    }
    E[w]=gm(2,dy4[0])^gm(3,dy4[1])^dy4[2]^dy4[3];
  }
}

// ─────────────────────────────────────────────────────────────────────
// Rebound: 10 params → (Δx_3, Δy_3, Δx_4) and a base Ref24. Returns 0
// if any DDT entry is zero (invalid base). Fills P->Dx3,Dy3,Dx4,Dz4 and R.
// ─────────────────────────────────────────────────────────────────────
static int rebound(Params10 *P, Ref24 *R){
  // Forward: Δx_2 → Δy_2 → Δx_3
  uint8_t Dy2[4];
  for(int r=0;r<4;r++){
    uint8_t Dx2=gm(MC[r][0],P->Din);
    Dy2[r]=SB[P->x2[r]^Dx2]^SB[P->x2[r]];
    if(Dy2[r]==0) return 0; // degenerate (shouldn't happen for Din≠0)
  }
  for(int r=0;r<4;r++) for(int c=0;c<4;c++){
    int rc=(4-c)&3;
    P->Dx3[r][c]=gm(MC[r][rc],Dy2[rc]);
  }
  // Backward: Δz_4 → Δx_4 → Δw_3(diag) → Δz_3 → Δy_3
  for(int r=0;r<4;r++){
    P->Dz4[r]=gm(iMC[r][0],P->Dout);
    uint8_t x4=iSB[P->z4[r]];
    P->Dx4[r]=iSB[P->z4[r]^P->Dz4[r]]^x4;
    R->x4[r]=x4;
    if(P->Dx4[r]==0) return 0;
  }
  for(int r=0;r<4;r++) for(int c=0;c<4;c++){
    int cc=((c-r)+4)&3; // column of z_3
    P->Dy3[r][c]=gm(iMC[r][cc],P->Dx4[cc]);
  }
  // Rebound middle: DDT(Δx_3,Δy_3) must all be nonzero; take base solution.
  for(int r=0;r<4;r++) for(int c=0;c<4;c++){
    uint8_t n=DDT[P->Dx3[r][c]][P->Dy3[r][c]];
    if(n==0) return 0;
    R->x3[r][c]=DDT_sol[P->Dx3[r][c]][P->Dy3[r][c]];
  }
  for(int r=0;r<4;r++) R->x2[r]=P->x2[r];
  return 1;
}

// Sanity: verify that base R actually satisfies SB(x3)^SB(x3^Dx3)=Dy3
static int check_rebound(const Params10 *P, const Ref24 *R){
  for(int r=0;r<4;r++) for(int c=0;c<4;c++){
    if((SB[R->x3[r][c]]^SB[R->x3[r][c]^P->Dx3[r][c]])!=P->Dy3[r][c]) return 0;
  }
  for(int r=0;r<4;r++){
    if((SB[R->x4[r]]^SB[R->x4[r]^P->Dx4[r]])!=P->Dz4[r]) return 0;
  }
  return 1;
}

// ─────────────────────────────────────────────────────────────────────
// Per-ω cache for incremental updates
// ─────────────────────────────────────────────────────────────────────
typedef struct {
  uint8_t dx3[4][4];
  uint8_t dy3[4][4];
  uint8_t dx4[4];
  uint8_t dy4[4];
  uint8_t E;
} OmegaCache;

static OmegaCache C[256];   // C[0] unused
static uint8_t sb_x2c[4], sb_x3c[4][4], sb_x4c[4]; // cached SB(refs)

// Initialize cache from Ref24 (this IS cold; instrumented so initial cost
// is accounted honestly, but amortizes to ~0 over 2^20).
static void cache_init(const Ref24 *R){
  for(int r=0;r<4;r++) sb_x2c[r]=SBi(R->x2[r]);
  for(int r=0;r<4;r++) for(int c=0;c<4;c++) sb_x3c[r][c]=SBi(R->x3[r][c]);
  for(int r=0;r<4;r++) sb_x4c[r]=SBi(R->x4[r]);
  for(int w=1;w<256;w++){
    uint8_t dy2[4];
    for(int r=0;r<4;r++) dy2[r]=SBi(R->x2[r]^gm(MC[r][0],w))^sb_x2c[r];
    for(int r=0;r<4;r++) for(int c=0;c<4;c++){
      int rc=(4-c)&3;
      C[w].dx3[r][c]=gm(MC[r][rc],dy2[rc]);
      C[w].dy3[r][c]=SBi(R->x3[r][c]^C[w].dx3[r][c])^sb_x3c[r][c];
    }
    for(int r=0;r<4;r++){
      uint8_t dx4=0;
      for(int rp=0;rp<4;rp++) dx4^=gm(MC[r][rp],C[w].dy3[rp][(r+rp)&3]);
      C[w].dx4[r]=dx4;
      C[w].dy4[r]=SBi(R->x4[r]^dx4)^sb_x4c[r];
    }
    C[w].E=gm(2,C[w].dy4[0])^gm(3,C[w].dy4[1])^C[w].dy4[2]^C[w].dy4[3];
  }
}

// z_4[r]-DDT-flip: toggle x_4[r] by Dx4_pair[r]. Per-ω: 1 SB.
static void flip_z4(Ref24 *R, int r, uint8_t Dx4_pair_r){
  R->x4[r]^=Dx4_pair_r;
  sb_x4c[r]=SBi(R->x4[r]);   // 1 SB (not per-ω; amortizes to 1/255)
  uint8_t mcr=MC[0][r];
  for(int w=1;w<256;w++){
    uint8_t nd=SBi(R->x4[r]^C[w].dx4[r])^sb_x4c[r];   // 1 SB/ω
    C[w].E^=gm(mcr,nd^C[w].dy4[r]);
    C[w].dy4[r]=nd;
  }
}

// x_3[p]-DDT-flip: toggle x_3[rp][cp] by Dx3_pair[rp][cp]. Per-ω: 2 SB.
static void flip_x3(Ref24 *R, int rp, int cp, uint8_t Dx3_pair_p){
  R->x3[rp][cp]^=Dx3_pair_p;
  sb_x3c[rp][cp]=SBi(R->x3[rp][cp]);   // 1 SB (not per-ω)
  int rs=((cp-rp)+4)&3; // which x_4 diag byte is affected
  uint8_t mc3=MC[rs][rp], mc4=MC[0][rs];
  for(int w=1;w<256;w++){
    uint8_t nd3=SBi(R->x3[rp][cp]^C[w].dx3[rp][cp])^sb_x3c[rp][cp]; // 1 SB/ω
    C[w].dx4[rs]^=gm(mc3,nd3^C[w].dy3[rp][cp]);
    C[w].dy3[rp][cp]=nd3;
    uint8_t nd4=SBi(R->x4[rs]^C[w].dx4[rs])^sb_x4c[rs];             // 1 SB/ω
    C[w].E^=gm(mc4,nd4^C[w].dy4[rs]);
    C[w].dy4[rs]=nd4;
  }
}

// ─────────────────────────────────────────────────────────────────────
// I_{m,n} fingerprint (POWINV-128b model: 1 lookup/ω → 13 power sums)
// ─────────────────────────────────────────────────────────────────────
static const int MSET[12]={3,5,7,11,13,19,21,25,37,41,49,81};
static uint8_t POWINV[256][16]; // POWINV[d][i] = (1/d)^k_i, k_i odd 1..25; 16-wide SIMD model

static void init_powinv(void){
  // collect odd k needed for MSET via Lucas: k submask of m, but we just
  // precompute 13 odd k's: 1,3,5,7,9,11,13,17,19,21,25,33,... — enough for MSET.
  // For SB-count purposes the content doesn't matter; 1 lookup/ω is the model.
  // But compute correctly so fingerprint check is meaningful.
  static const int Ks[13]={1,3,5,7,9,11,13,17,19,21,25,33,49};
  for(int d=0;d<256;d++){
    uint8_t inv=GINV[d];
    for(int i=0;i<13;i++){
      uint8_t p=1,b=inv; int e=Ks[i];
      while(e){ if(e&1)p=gm(p,b); b=gm(b,b); e>>=1; }
      POWINV[d][i]=p;
    }
    for(int i=13;i<16;i++) POWINV[d][i]=0;
  }
}

// Compute S_k accumulators from E[1..255] via 1 lookup/ω model.
// Returns a 128-bit (as 2x u64) XOR-fold fingerprint; NOT the full I_{m,n}
// (that's post-processing, negligible). Also counts lookups in g_fp_count.
static uint64_t g_fp_count=0;
static void fingerprint_Sk(const uint8_t E[256], uint8_t Sk[16]){
  memset(Sk,0,16);
  for(int w=1;w<256;w++){
    g_fp_count++;
    const uint8_t *row=POWINV[E[w]];
    // SIMD-128b XOR model: 1 wide-XOR. We do it byte-wise; counted as 1 lookup.
    for(int i=0;i<16;i++) Sk[i]^=row[i];
  }
}

// ─────────────────────────────────────────────────────────────────────
// Main: full 2^20 BRGC enumeration
// ─────────────────────────────────────────────────────────────────────
int main(int argc, char**argv){
  uint64_t seed = (argc>1)?strtoull(argv[1],0,0):0xC0FFEE;
  int verify_full = (argc>2)?atoi(argv[2]):0;
  int nbits = (argc>3)?atoi(argv[3]):20;
  if(nbits<4||nbits>20){fprintf(stderr,"nbits must be 4..20\n");return 1;}
  rng_s=seed?seed:1;
  init_tables();
  init_powinv();

  printf("=== DDT-GRAY r2: full 2^%d BRGC enumeration, seed=0x%lx, verify=%s ===\n",
         nbits,seed,verify_full?"FULL":"SAMPLE");

  // ── Find a valid 10-param base via rejection sampling ───────────────
  Params10 P; Ref24 R;
  int tries=0, ddt4_ct=0;
  for(;;){
    tries++;
    P.Din=rb()|1; P.Dout=rb()|1;
    for(int r=0;r<4;r++){P.x2[r]=rb(); P.z4[r]=rb();}
    if(!rebound(&P,&R)) continue;
    // count DDT=4 bytes (informational; gray covers one pair each)
    ddt4_ct=0;
    for(int r=0;r<4;r++)for(int c=0;c<4;c++)
      if(DDT[P.Dx3[r][c]][P.Dy3[r][c]]==4) ddt4_ct++;
    break;
  }
  assert(check_rebound(&P,&R));
  printf("[base] found valid 10-tuple after %d tries; x_3 DDT=4 bytes: %d/16\n",
         tries, ddt4_ct);
  printf("       Din=%02x Dout=%02x x2=[%02x,%02x,%02x,%02x] z4=[%02x,%02x,%02x,%02x]\n",
         P.Din,P.Dout,P.x2[0],P.x2[1],P.x2[2],P.x2[3],P.z4[0],P.z4[1],P.z4[2],P.z4[3]);

  // ── Self-test: cold_E vs mobius-honest.0-style 24-ref formula ───────
  {
    uint8_t E0[256]; cold_E(&R,E0);
    int nz=0; for(int w=1;w<256;w++) if(E0[w]) nz++;
    printf("[self] cold E' computed: %d/255 nonzero Δx_5[0] values\n",nz);
  }

  // ── Initialize incremental cache (instrumented) ────────────────────
  g_sb_count=0; g_fp_count=0;
  uint64_t sb_init_start=g_sb_count;
  cache_init(&R);
  uint64_t sb_init=g_sb_count-sb_init_start;
  printf("[init] cache_init SB ops: %lu (= %d + 255*%d expected %d)\n",
         sb_init, 24, 24, 24+255*24);

  // Verify cache matches cold at base
  {
    uint8_t E0[256]; cold_E(&R,E0);
    for(int w=1;w<256;w++) assert(C[w].E==E0[w]);
    printf("[self] incremental cache @ base == cold: 255/255 ✓\n");
  }

  // ── BRGC enumeration over nbits (bits 0..3 = z_4, bits 4..19 = x_3) ─
  uint64_t N=(uint64_t)1<<nbits;
  uint64_t sb_steps=0, n_z4=0, n_x3=0;
  uint64_t ver_ok=0, ver_fail=0, ver_done=0;
  uint64_t fp_hash=0; // running xor of all fingerprints (distinctness proxy)
  uint64_t checksum=0; // running checksum of all E arrays (detects any drift)
  uint8_t Sk[16];

  // checksum/fingerprint for base (entry 0)
  for(int w=1;w<256;w++) checksum += (uint64_t)C[w].E * (w*1099511628211ULL);
  fingerprint_Sk((uint8_t[256]){0},Sk); // dummy to keep g_fp_count accounting clean
  g_fp_count=0;
  {
    uint8_t Ecur[256]; Ecur[0]=0; for(int w=1;w<256;w++)Ecur[w]=C[w].E;
    fingerprint_Sk(Ecur,Sk);
    uint64_t h=0; for(int i=0;i<16;i++) h=h*131+Sk[i];
    fp_hash^=h;
  }

  // save base refs for final-state check
  Ref24 R0=R;

  uint64_t sb_before_loop=g_sb_count;

  for(uint64_t i=1;i<N;i++){
    int j=__builtin_ctzll(i);
    uint64_t sb0=g_sb_count;
    if(j<4){
      flip_z4(&R,j,P.Dx4[j]);
      n_z4++;
    } else {
      int p=j-4, rp=p&3, cp=p>>2;
      flip_x3(&R,rp,cp,P.Dx3[rp][cp]);
      n_x3++;
    }
    sb_steps += g_sb_count-sb0;

    // checksum (cheap, every entry)
    uint64_t cs=0;
    for(int w=1;w<256;w++) cs += (uint64_t)C[w].E * (w*1099511628211ULL);
    checksum ^= cs*(i|1);

    // fingerprint (1 lookup/ω model)
    {
      uint8_t Ecur[256]; Ecur[0]=0; for(int w=1;w<256;w++)Ecur[w]=C[w].E;
      fingerprint_Sk(Ecur,Sk);
      uint64_t h=0; for(int k=0;k<16;k++) h=h*131+Sk[k];
      fp_hash^=h;
    }

    // cold verification
    int do_verify = verify_full || (i<1000) || ((i&0xFFF)==0) || (i>N-1000);
    if(do_verify){
      uint8_t Ecold[256]; cold_E(&R,Ecold);
      int ok=1;
      for(int w=1;w<256;w++) if(C[w].E!=Ecold[w]){ok=0;break;}
      if(ok) ver_ok++; else ver_fail++;
      ver_done++;
    }

    if((i&((N>>4)-1))==0){
      printf("  [%5.1f%%] i=%lu  z4_flips=%lu x3_flips=%lu  SB/step=%.4f  ver=%lu/%lu\n",
             100.0*i/N, i, n_z4, n_x3, (double)sb_steps/i, ver_ok, ver_done);
    }
  }

  // BRGC over N=2^n visits every n-bit string once; at i=N-1 the gray code
  // is at position gray(N-1)=N/2 (binary 1000..0). To return to all-zeros
  // would need one more flip of bit n-1. We don't do it; the 24-refs at end
  // differ from base in exactly bit (n-1). Let's verify that.
  {
    uint64_t g=(N-1)^((N-1)>>1); // gray code of last index
    printf("[end] final gray state = 0x%05lx (expect single bit %d set: 0x%05lx)\n",
           g, nbits-1, (uint64_t)1<<(nbits-1));
    // Verify ref diff matches gray bits
    int ok=1;
    for(int r=0;r<4;r++){
      uint8_t exp=((g>>r)&1)?P.Dx4[r]:0;
      if((R.x4[r]^R0.x4[r])!=exp){ok=0;printf("  x4[%d] mismatch\n",r);}
    }
    for(int p=0;p<nbits-4;p++){
      int rp=p&3,cp=p>>2;
      uint8_t exp=((g>>(p+4))&1)?P.Dx3[rp][cp]:0;
      if((R.x3[rp][cp]^R0.x3[rp][cp])!=exp){ok=0;printf("  x3[%d][%d] mismatch\n",rp,cp);}
    }
    printf("[end] ref-state consistent with gray bits: %s\n",ok?"✓":"✗ FAIL");
  }

  // ── RESULTS ─────────────────────────────────────────────────────────
  uint64_t sb_loop=g_sb_count-sb_before_loop;
  uint64_t steps=N-1;
  double sb_per_step=(double)sb_steps/steps;
  double sb_per_omega=sb_per_step/255.0;
  double sb_per_entry_total=(double)(sb_init+sb_loop)/N; // incl. init, over all N entries
  double sb_per_omega_total=sb_per_entry_total/255.0;

  double fp_per_entry=(double)g_fp_count/N;
  double fp_per_omega=fp_per_entry/255.0;

  double c_off = sb_per_omega_total + fp_per_omega;
  // T_off = N_off × 255 × c_off / 160, N_off=2^{80} (DFJ) or 2^{88} (this paper).
  // The paper uses 2^{88}. Report both.
  double T_off_80 = 80 + log2(255.0*c_off/160.0);
  double T_off_88 = 88 + log2(255.0*c_off/160.0);

  printf("\n════════════════════════════════════════════════════════════\n");
  printf(" RESULTS (seed=0x%lx, N=2^%d=%lu entries)\n",seed,nbits,N);
  printf("════════════════════════════════════════════════════════════\n");
  printf(" Verification:        %lu / %lu entries cold-checked, %lu FAIL\n",
         ver_ok,ver_done,ver_fail);
  printf(" z_4-flips:           %lu (expect %lu)\n", n_z4, N-(N>>4));
  printf(" x_3-flips:           %lu (expect %lu)\n", n_x3, (N>>4)-1);
  printf(" SB ops (init):       %lu\n", sb_init);
  printf(" SB ops (loop):       %lu\n", sb_loop);
  printf(" SB / step:           %.4f  (analytic: %.4f)\n",
         sb_per_step, 255.0*((double)(N-(N>>4))*1+(double)((N>>4)-1)*2)/steps
                      + (double)(n_z4+n_x3)/steps /*ref-SB overhead*/);
  printf(" ── SB / (entry,ω):   %.5f  (step-only)\n", sb_per_omega);
  printf(" ── SB / (entry,ω):   %.5f  (incl. init, amortized over N)\n", sb_per_omega_total);
  printf("    [claim: 1.0625]\n");
  printf(" POWINV / (entry,ω):  %.5f  (fingerprint; claim: 1.000)\n", fp_per_omega);
  printf(" ── c_off:            %.4f  lookups/(entry,ω)  [claim: 2.06]\n", c_off);
  printf(" ── T_off @ N=2^80:   2^%.3f\n", T_off_80);
  printf(" ── T_off @ N=2^88:   2^%.3f  [paper claim: 2^89.72]\n", T_off_88);
  printf(" checksum:            0x%016lx\n", checksum);
  printf(" fp_hash:             0x%016lx\n", fp_hash);
  printf("════════════════════════════════════════════════════════════\n");

  return ver_fail?1:0;
}
