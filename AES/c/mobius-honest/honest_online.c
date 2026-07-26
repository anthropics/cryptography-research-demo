// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* honest_online.c — HONEST-ONLINE verification of the I_{m,n} parity-complete
 * bridge for 7-round AES-128 DS-MITM.
 *
 * For NKEYS random keys:
 *   1. Build δ-set on x_1[0] (256 plaintexts), encrypt to 7r ciphertexts.
 *   2. Record INTERNAL x_5[0], z_5[0] per element (ground truth).
 *   3. HONEST δ-set reconstruction: from P_ref + k_{-1}[0,5,10,15] ONLY,
 *      rebuild the 256 plaintexts → assert identical set.
 *   4. HONEST partial decryption: from C + k_6[0,7,10,13] ONLY, compute
 *      Δz_5[0] → assert byte-equal to internal Δz_5[0].
 *   5. H := {1/M^{-1}(Δz_5[0])} (honest), E' := {1/Δx_5[0]} (internal).
 *      Check I_{m,n} parity: even-match XOR odd-match, compare to k_A parity.
 *   6. FP: NWRONG random wrong k_6[idiag] → H_wrong, check no match.
 *
 * Build:  gcc -O3 -march=native -o honest_online honest_online.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

/* ---------------- GF(256) ---------------- */
static uint8_t gf_log[256], gf_alog[256];
static void gf_init(void){
    uint8_t x=1; for(int i=0;i<255;i++){ gf_alog[i]=x; gf_log[x]=i;
        x = (x<<1)^(x&0x80?0x1b:0)^x; /* ×3 */ }
    gf_alog[255]=1;
}
static inline uint8_t gf_mul(uint8_t a,uint8_t b){
    if(!a||!b) return 0; return gf_alog[(gf_log[a]+gf_log[b])%255];
}
static inline uint8_t gf_inv(uint8_t a){
    if(!a) return 0; return gf_alog[(255-gf_log[a])%255];
}
static inline uint8_t gf_pow(uint8_t a,int e){
    if(!a) return (e==0)?1:0; return gf_alog[((long)gf_log[a]*e)%255];
}

/* ---------------- AES primitives ---------------- */
static const uint8_t Sbox[256] = {
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
static uint8_t iSbox[256];
static uint8_t Minv_tab[256]; /* M^{-1} of the AES affine layer */

static void aes_tabs_init(void){
    for(int i=0;i<256;i++) iSbox[Sbox[i]]=i;
    /* Sbox(x)=M(inv(x))^0x63  ⟹  M^{-1}(w)=inv(iSbox[w^0x63]) */
    for(int w=0;w<256;w++) Minv_tab[w]=gf_inv(iSbox[w^0x63]);
    /* sanity: M^{-1} is GF(2)-linear */
    assert(Minv_tab[0]==0);
    for(int i=0;i<1000;i++){
        uint8_t a=rand()&255,b=rand()&255;
        assert(Minv_tab[a^b]==(Minv_tab[a]^Minv_tab[b]));
    }
    /* sanity: M^{-1}(Sbox(x)^Sbox(0)) = inv(x)^inv(0) = inv(x) */
    for(int x=1;x<256;x++)
        assert(Minv_tab[Sbox[x]^Sbox[0]]==gf_inv(x));
}

/* MixColumns on one column (in place, 4 bytes) */
static inline void mc_col(uint8_t *c){
    uint8_t a0=c[0],a1=c[1],a2=c[2],a3=c[3];
    c[0]=gf_mul(2,a0)^gf_mul(3,a1)^a2^a3;
    c[1]=a0^gf_mul(2,a1)^gf_mul(3,a2)^a3;
    c[2]=a0^a1^gf_mul(2,a2)^gf_mul(3,a3);
    c[3]=gf_mul(3,a0)^a1^a2^gf_mul(2,a3);
}
static inline void imc_col(uint8_t *c){
    uint8_t a0=c[0],a1=c[1],a2=c[2],a3=c[3];
    c[0]=gf_mul(0x0e,a0)^gf_mul(0x0b,a1)^gf_mul(0x0d,a2)^gf_mul(0x09,a3);
    c[1]=gf_mul(0x09,a0)^gf_mul(0x0e,a1)^gf_mul(0x0b,a2)^gf_mul(0x0d,a3);
    c[2]=gf_mul(0x0d,a0)^gf_mul(0x09,a1)^gf_mul(0x0e,a2)^gf_mul(0x0b,a3);
    c[3]=gf_mul(0x0b,a0)^gf_mul(0x0d,a1)^gf_mul(0x09,a2)^gf_mul(0x0e,a3);
}
/* ShiftRows (column-major byte layout: s[4*c+r]) */
static inline void sr(uint8_t *s){
    uint8_t t[16]; for(int c=0;c<4;c++)for(int r=0;r<4;r++) t[4*c+r]=s[4*((c+r)%4)+r];
    memcpy(s,t,16);
}
static inline void isr(uint8_t *s){
    uint8_t t[16]; for(int c=0;c<4;c++)for(int r=0;r<4;r++) t[4*c+r]=s[4*((c+4-r)%4)+r];
    memcpy(s,t,16);
}
static inline void sb(uint8_t *s){for(int i=0;i<16;i++)s[i]=Sbox[s[i]];}
static inline void isb(uint8_t *s){for(int i=0;i<16;i++)s[i]=iSbox[s[i]];}
static inline void ark(uint8_t *s,const uint8_t *k){for(int i=0;i<16;i++)s[i]^=k[i];}
static inline void mc(uint8_t *s){for(int c=0;c<4;c++)mc_col(s+4*c);}
static inline void imc(uint8_t *s){for(int c=0;c<4;c++)imc_col(s+4*c);}

/* Key schedule (bytewise, column-major) */
static const uint8_t Rcon[11]={0,0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1b,0x36};
static void key_expand(const uint8_t key[16], uint8_t rk[11][16]){
    memcpy(rk[0],key,16);
    for(int r=1;r<=10;r++){
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        /* RotWord + SubWord */
        uint8_t tt=t[0]; t[0]=Sbox[t[1]];t[1]=Sbox[t[2]];t[2]=Sbox[t[3]];t[3]=Sbox[tt];
        t[0]^=Rcon[r];
        for(int i=0;i<4;i++) rk[r][i]=rk[r-1][i]^t[i];
        for(int c=1;c<4;c++)for(int i=0;i<4;i++) rk[r][4*c+i]=rk[r-1][4*c+i]^rk[r][4*(c-1)+i];
    }
}

/* 7-round encrypt with full internal trace. trace[i] = x_i (i=0..6), trace[7]=C */
static void enc7_trace(const uint8_t rk[11][16], const uint8_t P[16],
                       uint8_t trace[8][16]){
    uint8_t s[16]; memcpy(s,P,16); ark(s,rk[0]);
    memcpy(trace[0],s,16);
    for(int r=0;r<6;r++){
        sb(s); sr(s); mc(s); ark(s,rk[r+1]);
        memcpy(trace[r+1],s,16);
    }
    sb(s); sr(s); ark(s,rk[7]);
    memcpy(trace[7],s,16);
}

/* ---------------- I_{m,n} fingerprint ---------------- */
/* P_m(S) = XOR_{i<j} (s_i ^ s_j)^m, computed via power-sum identity:
 * In char 2, Σ_{i<j}(a_i+a_j)^m = Σ_{k: C(m,k) odd, 0<k<m} S_k · S_{m-k} / ... no.
 * Direct O(n^2) is fine for n≤256. */
#define NPM 12
static const int Pm_exp[NPM] = {3,5,7,11,13,19,21,25,37,41,49,81}; /* popcount≥2 */

static void compute_Pm(const uint8_t *S, int n, uint8_t Pm[NPM]){
    memset(Pm,0,NPM);
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
        uint8_t d=S[i]^S[j]; if(!d) continue;
        int lg=gf_log[d];
        for(int k=0;k<NPM;k++) Pm[k]^=gf_alog[(lg*Pm_exp[k])%255];
    }
}
/* Fingerprint: I_{m,n} = P_m^n · P_n^{-m} for all ordered (m,n), m≠n.
 * Pack into 64-bit hash. Returns 0 if ALL P_m are 0 (degenerate). */
static uint64_t fingerprint(const uint8_t Pm[NPM]){
    /* find first nonzero as normalizer */
    int n0=-1; for(int k=0;k<NPM;k++) if(Pm[k]){n0=k;break;}
    if(n0<0) return 0;
    uint64_t fp=1;
    int n=Pm_exp[n0]; uint8_t Pn=Pm[n0];
    for(int k=0;k<NPM;k++){
        if(k==n0) continue;
        int m=Pm_exp[k];
        /* I = P_m^n · P_n^{-m} */
        uint8_t I = gf_mul(gf_pow(Pm[k],n), gf_pow(gf_inv(Pn),m));
        fp = fp*257 + I;
    }
    fp = fp*257 + (uint8_t)(n0+1); /* encode which normalizer */
    return fp;
}
/* compute fp(S) and fp(S∪{0}) */
static void fp_pair(const uint8_t *S, int n, uint64_t *fp_even, uint64_t *fp_odd){
    uint8_t Pm[NPM]; compute_Pm(S,n,Pm);
    *fp_even = fingerprint(Pm);
    /* S∪{0}: adds pairs (0,s_i) with diff s_i. Incremental: Pm += Σ s_i^m */
    uint8_t Pm2[NPM]; memcpy(Pm2,Pm,NPM);
    for(int i=0;i<n;i++){ uint8_t s=S[i]; if(!s) continue;
        int lg=gf_log[s];
        for(int k=0;k<NPM;k++) Pm2[k]^=gf_alog[(lg*Pm_exp[k])%255];
    }
    *fp_odd = fingerprint(Pm2);
}

/* ---------------- main test ---------------- */
static uint64_t rng_state=0x243f6a8885a308d3ULL;
static uint32_t rng(void){ rng_state=rng_state*6364136223846793005ULL+1; return rng_state>>32; }

int main(int argc,char**argv){
    int NKEYS=50, NWRONG=20;
    if(argc>1) NKEYS=atoi(argv[1]);
    if(argc>2) rng_state^=strtoull(argv[2],0,0)*0x9e3779b97f4a7c15ULL;
    gf_init(); aes_tabs_init();

    /* self-test encryption vs known vector */
    {
        uint8_t key[16]={0},rk[11][16],tr[8][16],pt[16]={0};
        key_expand(key,rk);
        /* full 10r check: encrypt 0 with key 0 → 66e94bd4ef8a2c3b884cfa59ca342b2e */
        uint8_t s[16]={0}; ark(s,rk[0]);
        for(int r=0;r<9;r++){sb(s);sr(s);mc(s);ark(s,rk[r+1]);}
        sb(s);sr(s);ark(s,rk[10]);
        static const uint8_t kat[16]={0x66,0xe9,0x4b,0xd4,0xef,0x8a,0x2c,0x3b,0x88,0x4c,0xfa,0x59,0xca,0x34,0x2b,0x2e};
        assert(memcmp(s,kat,16)==0);
        (void)tr;(void)pt;
    }

    int n_dset_ok=0, n_pdec_ok=0, n_off24_ok=0;
    int n_even=0,n_odd=0,n_either=0,n_a0=0,n_parity_correct=0;
    int n_fp=0, n_fp_trials=0;
    int kA_hist[8]={0};

    for(int trial=0;trial<NKEYS;trial++){
        uint8_t key[16],rk[11][16];
        for(int i=0;i<16;i++) key[i]=rng();
        key_expand(key,rk);
        /* k_{-1}=rk[0], k_0=rk[1], ..., k_6=rk[7] */
        const uint8_t *km1=rk[0], *k0=rk[1], *k6=rk[7];

        /* ---- Build δ-set via full key: random base x_1, vary x_1[0]=v ---- */
        uint8_t x1_base[16]; for(int i=0;i<16;i++) x1_base[i]=rng();
        uint8_t P[256][16], C[256][16];
        uint8_t x5_0[256], z5_0[256];
        for(int v=0;v<256;v++){
            uint8_t x1[16]; memcpy(x1,x1_base,16); x1[0]=v;
            /* backward to plaintext: x1 → w0=x1^k0 → z0=iMC(w0) → y0=iSR(z0) → x0=iSB(y0) → P=x0^km1 */
            uint8_t s[16]; memcpy(s,x1,16); ark(s,k0);
            imc(s); isr(s); isb(s); ark(s,km1);
            memcpy(P[v],s,16);
            /* forward encrypt with trace */
            uint8_t tr[8][16]; enc7_trace(rk,P[v],tr);
            memcpy(C[v],tr[7],16);
            x5_0[v]=tr[5][0];
            z5_0[v]=Sbox[tr[5][0]]; /* z_5[0]=y_5[0]=SB(x_5[0]) since SR row 0 */
            /* sanity: tr[1] should equal x1 */
            assert(memcmp(tr[1],x1,16)==0);
        }
        /* sanity: all P[v] equal on non-diag bytes */
        for(int v=1;v<256;v++)for(int i=0;i<16;i++)
            if(i!=0&&i!=5&&i!=10&&i!=15) assert(P[v][i]==P[0][i]);

        /* ---- Pick reference: v_ref = x1_base[0] (the "t=0" element). Actually
         * use index 0 of the array, i.e. v=0 ⟹ x_1[0]=0. But any works.
         * In the real attack the ref is the right-pair plaintext; here arbitrary. */
        int vref = x1_base[0]; /* so that "t=0" below corresponds to P[vref] */
        uint8_t a = x5_0[vref]; /* = x_5[0]_ref */

        /* ---- Test 3: HONEST δ-set reconstruction from k_{-1}[diag]+P_ref ---- */
        {
            const int diag[4]={0,5,10,15};
            const uint8_t coef[4]={0x0e,0x09,0x0d,0x0b};
            uint8_t y0ref[4];
            for(int i=0;i<4;i++) y0ref[i]=Sbox[P[vref][diag[i]]^km1[diag[i]]];
            /* mark which P[v] appear */
            int seen[256]={0};
            for(int t=0;t<256;t++){
                uint8_t Pt[16]; memcpy(Pt,P[vref],16);
                for(int i=0;i<4;i++)
                    Pt[diag[i]]=iSbox[y0ref[i]^gf_mul(coef[i],t)]^km1[diag[i]];
                /* find which v this is */
                int found=-1;
                for(int v=0;v<256;v++) if(memcmp(Pt,P[v],16)==0){found=v;break;}
                if(found>=0) seen[found]++;
            }
            int ok=1; for(int v=0;v<256;v++) if(seen[v]!=1){ok=0;break;}
            if(ok) n_dset_ok++;
            else printf("  [key %d] δ-set reconstruction MISMATCH\n",trial);
        }

        /* ---- Test 4: HONEST Δz_5[0] from C + k_6[idiag] ---- */
        /* pos[r] = byte index of z_6 feeding y_6[r,0] */
        const int pos[4]={0,13,10,7};
        uint8_t z5hon[256]; int pdec_ok=1;
        for(int v=0;v<256;v++){
            uint8_t x6c0[4];
            for(int r=0;r<4;r++) x6c0[r]=iSbox[C[v][pos[r]]^k6[pos[r]]];
            uint8_t z=gf_mul(0x0e,x6c0[0])^gf_mul(0x0b,x6c0[1])^
                      gf_mul(0x0d,x6c0[2])^gf_mul(0x09,x6c0[3]);
            z5hon[v]=z; /* = z_5[0] ^ u_5[0]; diff cancels u_5 */
        }
        for(int v=0;v<256;v++){
            uint8_t dhon=z5hon[v]^z5hon[vref];
            uint8_t dint=z5_0[v]^z5_0[vref];
            if(dhon!=dint){pdec_ok=0;break;}
        }
        if(pdec_ok) n_pdec_ok++;
        else printf("  [key %d] partial-decrypt Δz_5[0] MISMATCH\n",trial);

        /* ---- Test 5: OFFLINE E' from 24-byte parametrization ONLY ----
         * Extract x_2[col0], x_3[full], x_4[diag] at ref. Propagate δ-set
         * forward using ONLY these 24 bytes (no round keys). Verify
         * {Δw_4[0]_ω : ω=Δy_1[0]∈GF(256)} == internal {Δx_5[0]}. */
        {
            uint8_t tr_ref[8][16]; enc7_trace(rk,P[vref],tr_ref);
            uint8_t x2c0[4]={tr_ref[2][0],tr_ref[2][1],tr_ref[2][2],tr_ref[2][3]};
            uint8_t x3[16]; memcpy(x3,tr_ref[3],16);
            uint8_t x4d[4]={tr_ref[4][0],tr_ref[4][5],tr_ref[4][10],tr_ref[4][15]};
            /* δ-set indexed by ω=Δy_1[0]: x_2[r,0](ω)=x2c0[r]⊕MC[r,0]·ω */
            uint8_t Dw4off[256]={0};
            static const uint8_t mc0[4]={2,1,1,3};
            for(int om=0;om<256;om++){
                /* Δx_2[col0]=ω·[2,1,1,3] */
                uint8_t x2w[4],dy2[4];
                for(int r=0;r<4;r++){x2w[r]=x2c0[r]^gf_mul(mc0[r],om);
                                     dy2[r]=Sbox[x2w[r]]^Sbox[x2c0[r]];}
                /* Δz_2: SR of col0-only: Δz_2[r,(4-r)%4]=dy2[r] */
                uint8_t dz2[16]={0};for(int r=0;r<4;r++)dz2[4*((4-r)%4)+r]=dy2[r];
                uint8_t dx3[16];memcpy(dx3,dz2,16);mc(dx3);
                uint8_t x3w[16];for(int i=0;i<16;i++)x3w[i]=x3[i]^dx3[i];
                uint8_t dy3[16];for(int i=0;i<16;i++)dy3[i]=Sbox[x3w[i]]^Sbox[x3[i]];
                /* Δz_3=SR(dy3), Δw_3=MC(Δz_3), Δx_4=Δw_3 */
                uint8_t t[16];memcpy(t,dy3,16);sr(t);mc(t);
                /* x_4[diag](ω)=x4d⊕Δx_4[diag]; Δy_4[diag]=SB(x_4[diag](ω))⊕SB(x4d) */
                uint8_t dy4d[4];static const int diag[4]={0,5,10,15};
                for(int r=0;r<4;r++){uint8_t x4w=x4d[r]^t[diag[r]];
                                     dy4d[r]=Sbox[x4w]^Sbox[x4d[r]];}
                /* Δz_4[col0]=dy4d (SR of diag→col0); Δw_4[0]=MC[row0]·Δz_4[col0] */
                Dw4off[om]=gf_mul(2,dy4d[0])^gf_mul(3,dy4d[1])^dy4d[2]^dy4d[3];
            }
            /* Compare as multisets against internal {Δx_5[0]_v:v} */
            int hI[256]={0},hO[256]={0};
            for(int v=0;v<256;v++)hI[x5_0[v]^x5_0[vref]]++;
            for(int om=0;om<256;om++)hO[Dw4off[om]]++;
            int ok=1;for(int i=0;i<256;i++)if(hI[i]!=hO[i]){ok=0;break;}
            if(ok)n_off24_ok++;
            else if(trial<3)printf("  [key %d] OFFLINE-24 MISMATCH\n",trial);
        }

        /* ---- Build H (honest) and E' (internal) ---- */
        uint8_t H[256],E[256]; int nH=0,nE=0;
        for(int v=0;v<256;v++){
            if(v==vref) continue;
            uint8_t dz=z5hon[v]^z5hon[vref];
            uint8_t g=Minv_tab[dz];
            if(g) H[nH++]=gf_inv(g);
            uint8_t dx=x5_0[v]^x5_0[vref];
            if(dx) E[nE++]=gf_inv(dx);
        }
        /* sanity: g=0 ⟺ dz=0 ⟺ dx=0 (same exclusion set) */
        assert(nH==nE);

        int kA=0; for(int v=0;v<256;v++) if(x5_0[v]==0) kA++;
        if(kA<8) kA_hist[kA]++;

        /* ---- I_{m,n} parity match ---- */
        uint64_t fpH_ev,fpH_od,fpE_ev,fpE_od;
        fp_pair(H,nH,&fpH_ev,&fpH_od);
        fp_pair(E,nE,&fpE_ev,&fpE_od);
        int ev=(fpH_ev==fpE_ev); /* 0==0 counts as (degenerate) match */
        int od=(fpH_od==fpE_od);
        if(a==0){ n_a0++; /* edge: expect neither in general */ }
        if(ev) n_even++;
        if(od) n_odd++;
        if(ev||od) n_either++;
        /* parity-lemma check: ev ⟺ kA even, od ⟺ kA odd (when a≠0) */
        if(a!=0){
            int pred_ev=(kA%2==0), pred_od=(kA%2==1);
            if(ev==pred_ev && od==pred_od) n_parity_correct++;
            else {
                printf("  [key %d] PARITY MISMATCH: kA=%d ev=%d od=%d a=%02x nH=%d nE=%d\n",
                        trial,kA,ev,od,a,nH,nE);
                uint8_t PmH[NPM],PmE[NPM];
                compute_Pm(H,nH,PmH); compute_Pm(E,nE,PmE);
                printf("    m:   "); for(int k=0;k<NPM;k++)printf("%3d ",Pm_exp[k]); printf("\n");
                printf("    PmH: "); for(int k=0;k<NPM;k++)printf("%02x  ",PmH[k]); printf("\n");
                printf("    PmE: "); for(int k=0;k<NPM;k++)printf("%02x  ",PmE[k]); printf("\n");
                printf("    a^2m·PmE: ");
                for(int k=0;k<NPM;k++)printf("%02x  ",gf_mul(gf_pow(a,2*Pm_exp[k]),PmE[k]));
                printf("\n");
                printf("    fpH_ev=%016llx fpE_ev=%016llx\n",(unsigned long long)fpH_ev,(unsigned long long)fpE_ev);
                printf("    fpH_od=%016llx fpE_od=%016llx\n",(unsigned long long)fpH_od,(unsigned long long)fpE_od);
                /* check H == T_a(E') elementwise (sorted) */
                uint8_t TE[256]; for(int i=0;i<nE;i++) TE[i]=a^gf_mul(gf_mul(a,a),E[i]);
                /* count mismatches as multisets */
                int histH[256]={0},histT[256]={0};
                for(int i=0;i<nH;i++)histH[H[i]]++;
                for(int i=0;i<nE;i++)histT[TE[i]]++;
                int diff=0; for(int i=0;i<256;i++) diff+=abs(histH[i]-histT[i]);
                printf("    |H Δ T_a(E')| (multiset sym-diff) = %d\n",diff);
            }
        }

        /* ---- FP: wrong k_{-1}[diag] (wrong δ-set) + true k_6 ---- */
        for(int w=0;w<NWRONG;w++){
            const int diag[4]={0,5,10,15};
            const uint8_t coef[4]={0x0e,0x09,0x0d,0x0b};
            uint8_t wkm1[4]={rng(),rng(),rng(),rng()};
            if(wkm1[0]==km1[0]) wkm1[0]^=1;
            uint8_t y0ref[4];
            for(int i=0;i<4;i++) y0ref[i]=Sbox[P[vref][diag[i]]^wkm1[i]];
            /* build wrong δ-set (different plaintexts, but still in structure) */
            uint8_t Hw[256]; int nHw=0;
            uint8_t zref2=z5hon[vref]; /* t=0 ⟹ P_ref itself, same ref */
            for(int t=1;t<256;t++){
                uint8_t Pt[16]; memcpy(Pt,P[vref],16);
                for(int i=0;i<4;i++)
                    Pt[diag[i]]=iSbox[y0ref[i]^gf_mul(coef[i],t)]^wkm1[i];
                uint8_t tr[8][16]; enc7_trace(rk,Pt,tr);
                uint8_t x6c0[4];
                for(int r=0;r<4;r++) x6c0[r]=iSbox[tr[7][pos[r]]^k6[pos[r]]];
                uint8_t z=gf_mul(0x0e,x6c0[0])^gf_mul(0x0b,x6c0[1])^
                          gf_mul(0x0d,x6c0[2])^gf_mul(0x09,x6c0[3]);
                uint8_t g=Minv_tab[z^zref2];
                if(g) Hw[nHw++]=gf_inv(g);
            }
            uint64_t fev,fod; fp_pair(Hw,nHw,&fev,&fod);
            n_fp_trials++;
            if((fev==fpE_ev&&fev)||(fod==fpE_od&&fod)) n_fp++;
        }

        /* ---- FP: wrong k_6[idiag] ---- */
        for(int w=0;w<NWRONG;w++){
            uint8_t wk6[4]={rng(),rng(),rng(),rng()};
            /* ensure != true */
            if(wk6[0]==k6[0]&&wk6[1]==k6[13]&&wk6[2]==k6[10]&&wk6[3]==k6[7]) wk6[0]^=1;
            uint8_t Hw[256]; int nHw=0;
            uint8_t zref=0;
            {
                uint8_t x6c0[4];
                x6c0[0]=iSbox[C[vref][0]^wk6[0]];
                x6c0[1]=iSbox[C[vref][13]^wk6[1]];
                x6c0[2]=iSbox[C[vref][10]^wk6[2]];
                x6c0[3]=iSbox[C[vref][7]^wk6[3]];
                zref=gf_mul(0x0e,x6c0[0])^gf_mul(0x0b,x6c0[1])^
                     gf_mul(0x0d,x6c0[2])^gf_mul(0x09,x6c0[3]);
            }
            for(int v=0;v<256;v++){
                if(v==vref) continue;
                uint8_t x6c0[4];
                x6c0[0]=iSbox[C[v][0]^wk6[0]];
                x6c0[1]=iSbox[C[v][13]^wk6[1]];
                x6c0[2]=iSbox[C[v][10]^wk6[2]];
                x6c0[3]=iSbox[C[v][7]^wk6[3]];
                uint8_t z=gf_mul(0x0e,x6c0[0])^gf_mul(0x0b,x6c0[1])^
                          gf_mul(0x0d,x6c0[2])^gf_mul(0x09,x6c0[3]);
                uint8_t g=Minv_tab[z^zref];
                if(g) Hw[nHw++]=gf_inv(g);
            }
            uint64_t fev,fod; fp_pair(Hw,nHw,&fev,&fod);
            n_fp_trials++;
            if((fev==fpE_ev&&fev)||(fod==fpE_od&&fod)) n_fp++;
        }
    }

    printf("\n========== HONEST-ONLINE RESULTS (NKEYS=%d) ==========\n",NKEYS);
    printf("δ-set reconstruction (k_{-1}[diag] only) : %d/%d\n",n_dset_ok,NKEYS);
    printf("Δz_5[0] partial-decrypt == internal      : %d/%d\n",n_pdec_ok,NKEYS);
    printf("OFFLINE E' from 24-byte params == intern : %d/%d\n",n_off24_ok,NKEYS);
    printf("I_{m,n} even-match                       : %d/%d\n",n_even,NKEYS);
    printf("I_{m,n} odd-match (∪{0})                 : %d/%d\n",n_odd,NKEYS);
    printf("EITHER match                             : %d/%d\n",n_either,NKEYS);
    printf("  (a=0 edge cases: %d)\n",n_a0);
    printf("Parity-Lemma exact (ev⟺kA even,a≠0)     : %d/%d\n",n_parity_correct,NKEYS-n_a0);
    printf("Wrong-k_6 FP                             : %d/%d\n",n_fp,n_fp_trials);
    printf("k_A histogram [0..7]: ");
    for(int i=0;i<8;i++)printf("%d ",kA_hist[i]);printf("\n");
    return 0;
}
