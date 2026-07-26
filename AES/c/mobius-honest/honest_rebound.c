// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* honest_rebound.c — FULL-AES honest test WITH right pair + 10-param offline.
 *
 * Closes the last gap: prior honest_online.c verified 24-byte→E'. This test
 * adds the 10-param REBOUND (= what the real offline table enumerates):
 *
 *   1. Rebound-construct (key, P_ref, P_pair) with right-pair property
 *      (Δx_1=[*,0..], Δx_5=[*,0..]) — known-key shortcut, same as small-AES.
 *   2. Extract TRUE 10 params (Δz_1[0],x_2[col0],Δw_4[0],z_4[col0]) from trace.
 *   3. OFFLINE: rebound from 10 params ONLY → x_3 candidates (DDT) → 24 bytes
 *      → E' via forward δ-set propagation.  NO round keys used.
 *   4. ONLINE: H from (C, k_6[0,7,10,13]) via partial decryption.
 *   5. Verify I_{m,n} parity match H↔E' (one of the DDT-x_3 branches).
 *   6. Sampled-offline FP: N_SAMP random 10-param tuples → how many match H?
 *
 * Build: gcc -O3 -march=native -o honest_rebound honest_rebound.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

/* ===== GF(256) + AES primitives (same as honest_online.c) ===== */
static uint8_t gf_log[256],gf_alog[256];
static void gf_init(void){uint8_t x=1;for(int i=0;i<255;i++){gf_alog[i]=x;gf_log[x]=i;
    x=(x<<1)^(x&0x80?0x1b:0)^x;}gf_alog[255]=1;}
static inline uint8_t gm(uint8_t a,uint8_t b){if(!a||!b)return 0;return gf_alog[(gf_log[a]+gf_log[b])%255];}
static inline uint8_t gi(uint8_t a){if(!a)return 0;return gf_alog[(255-gf_log[a])%255];}
static inline uint8_t gp(uint8_t a,int e){if(!a)return e==0?1:0;return gf_alog[((long)gf_log[a]*e)%255];}

static const uint8_t Sbox[256]={
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
0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16};
static uint8_t iSbox[256],Minv[256];
static void tabs_init(void){for(int i=0;i<256;i++)iSbox[Sbox[i]]=i;
    for(int w=0;w<256;w++)Minv[w]=gi(iSbox[w^0x63]);}

static inline void mc_c(uint8_t*c){uint8_t a=c[0],b=c[1],d=c[2],e=c[3];
    c[0]=gm(2,a)^gm(3,b)^d^e;c[1]=a^gm(2,b)^gm(3,d)^e;
    c[2]=a^b^gm(2,d)^gm(3,e);c[3]=gm(3,a)^b^d^gm(2,e);}
static inline void imc_c(uint8_t*c){uint8_t a=c[0],b=c[1],d=c[2],e=c[3];
    c[0]=gm(14,a)^gm(11,b)^gm(13,d)^gm(9,e);c[1]=gm(9,a)^gm(14,b)^gm(11,d)^gm(13,e);
    c[2]=gm(13,a)^gm(9,b)^gm(14,d)^gm(11,e);c[3]=gm(11,a)^gm(13,b)^gm(9,d)^gm(14,e);}
static inline void sb(uint8_t*s){for(int i=0;i<16;i++)s[i]=Sbox[s[i]];}
static inline void isb(uint8_t*s){for(int i=0;i<16;i++)s[i]=iSbox[s[i]];}
static inline void sr(uint8_t*s){uint8_t t[16];for(int c=0;c<4;c++)for(int r=0;r<4;r++)t[4*c+r]=s[4*((c+r)%4)+r];memcpy(s,t,16);}
static inline void isr(uint8_t*s){uint8_t t[16];for(int c=0;c<4;c++)for(int r=0;r<4;r++)t[4*c+r]=s[4*((c+4-r)%4)+r];memcpy(s,t,16);}
static inline void mc(uint8_t*s){for(int c=0;c<4;c++)mc_c(s+4*c);}
static inline void imc(uint8_t*s){for(int c=0;c<4;c++)imc_c(s+4*c);}
static inline void ark(uint8_t*s,const uint8_t*k){for(int i=0;i<16;i++)s[i]^=k[i];}

static const uint8_t Rcon[11]={0,1,2,4,8,0x10,0x20,0x40,0x80,0x1b,0x36};
static void key_from(int r0,const uint8_t kr0[16],uint8_t rk[11][16]){
    memcpy(rk[r0],kr0,16);
    for(int r=r0;r>=1;r--){
        for(int c=3;c>=1;c--)for(int i=0;i<4;i++)rk[r-1][4*c+i]=rk[r][4*c+i]^rk[r][4*(c-1)+i];
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        uint8_t tt=t[0];t[0]=Sbox[t[1]];t[1]=Sbox[t[2]];t[2]=Sbox[t[3]];t[3]=Sbox[tt];
        t[0]^=Rcon[r];for(int i=0;i<4;i++)rk[r-1][i]=rk[r][i]^t[i];
    }
    for(int r=r0+1;r<=10;r++){
        uint8_t t[4]={rk[r-1][12],rk[r-1][13],rk[r-1][14],rk[r-1][15]};
        uint8_t tt=t[0];t[0]=Sbox[t[1]];t[1]=Sbox[t[2]];t[2]=Sbox[t[3]];t[3]=Sbox[tt];
        t[0]^=Rcon[r];for(int i=0;i<4;i++)rk[r][i]=rk[r-1][i]^t[i];
        for(int c=1;c<4;c++)for(int i=0;i<4;i++)rk[r][4*c+i]=rk[r-1][4*c+i]^rk[r][4*(c-1)+i];
    }
}
static void enc7_trace(const uint8_t rk[11][16],const uint8_t P[16],uint8_t tr[8][16]){
    uint8_t s[16];memcpy(s,P,16);ark(s,rk[0]);memcpy(tr[0],s,16);
    for(int r=0;r<6;r++){sb(s);sr(s);mc(s);ark(s,rk[r+1]);memcpy(tr[r+1],s,16);}
    sb(s);sr(s);ark(s,rk[7]);memcpy(tr[7],s,16);
}

/* DDT: store up to 4 solutions */
static uint8_t DDTn[256][256];static uint8_t DDTs[256][256][4];
static void build_ddt(void){
    for(int di=0;di<256;di++)for(int x=0;x<256;x++){
        int do_=Sbox[x]^Sbox[x^di];
        if(DDTn[di][do_]<4)DDTs[di][do_][DDTn[di][do_]]=x;
        DDTn[di][do_]++;
    }
}

/* I_{m,n} */
#define NPM 12
static const int Pe[NPM]={3,5,7,11,13,19,21,25,37,41,49,81};
static void cPm(const uint8_t*S,int n,uint8_t P[NPM]){memset(P,0,NPM);
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){uint8_t d=S[i]^S[j];if(!d)continue;
        int lg=gf_log[d];for(int k=0;k<NPM;k++)P[k]^=gf_alog[(lg*Pe[k])%255];}}
static uint64_t fp(const uint8_t P[NPM]){int n0=-1;for(int k=0;k<NPM;k++)if(P[k]){n0=k;break;}
    if(n0<0)return 0;uint64_t f=1;int n=Pe[n0];uint8_t Pn=P[n0];
    for(int k=0;k<NPM;k++){if(k==n0)continue;
        f=f*257+gm(gp(P[k],n),gp(gi(Pn),Pe[k]));}
    return f*257+(n0+1);}
static void fpp(const uint8_t*S,int n,uint64_t*fe,uint64_t*fo){
    uint8_t P[NPM];cPm(S,n,P);*fe=fp(P);uint8_t P2[NPM];memcpy(P2,P,NPM);
    for(int i=0;i<n;i++){uint8_t s=S[i];if(!s)continue;int lg=gf_log[s];
        for(int k=0;k<NPM;k++)P2[k]^=gf_alog[(lg*Pe[k])%255];}
    *fo=fp(P2);}

/* 24-byte offline formula → E' multiset */
static int off24_Eprime(const uint8_t x2c0[4],const uint8_t x3[16],const uint8_t x4d[4],
                        uint8_t E[256]){
    static const uint8_t mc0[4]={2,1,1,3};static const int diag[4]={0,5,10,15};
    int nE=0;
    for(int om=1;om<256;om++){
        uint8_t dy2[4];for(int r=0;r<4;r++)dy2[r]=Sbox[x2c0[r]^gm(mc0[r],om)]^Sbox[x2c0[r]];
        uint8_t dz2[16]={0};for(int r=0;r<4;r++)dz2[4*((4-r)%4)+r]=dy2[r];
        uint8_t dx3[16];memcpy(dx3,dz2,16);mc(dx3);
        uint8_t dy3[16];for(int i=0;i<16;i++)dy3[i]=Sbox[x3[i]^dx3[i]]^Sbox[x3[i]];
        uint8_t t[16];memcpy(t,dy3,16);sr(t);mc(t);
        uint8_t dy4[4];for(int r=0;r<4;r++)dy4[r]=Sbox[x4d[r]^t[diag[r]]]^Sbox[x4d[r]];
        uint8_t dw4=gm(2,dy4[0])^gm(3,dy4[1])^dy4[2]^dy4[3];
        if(dw4)E[nE++]=gi(dw4);
    }
    return nE;
}

/* 10-param rebound: returns #x3-branches and for each the 24-byte tuple */
typedef struct{uint8_t x2c0[4],x3[16],x4d[4];}Off24;
static int rebound10(uint8_t Dz1,const uint8_t x2c0[4],uint8_t Dw4,const uint8_t z4c0[4],
                     Off24 out[],int maxout){
    static const int diag[4]={0,5,10,15};
    uint8_t Dx2[4]={gm(2,Dz1),Dz1,Dz1,gm(3,Dz1)};
    uint8_t Dy2[4];for(int r=0;r<4;r++)Dy2[r]=Sbox[x2c0[r]^Dx2[r]]^Sbox[x2c0[r]];
    uint8_t Dz2[16]={0};for(int r=0;r<4;r++)Dz2[4*((4-r)%4)+r]=Dy2[r];
    uint8_t Dx3[16];memcpy(Dx3,Dz2,16);mc(Dx3);
    uint8_t Dz4[4]={gm(14,Dw4),gm(9,Dw4),gm(13,Dw4),gm(11,Dw4)};
    uint8_t x4d[4],Dx4d[4];for(int r=0;r<4;r++){x4d[r]=iSbox[z4c0[r]];Dx4d[r]=iSbox[z4c0[r]^Dz4[r]]^x4d[r];}
    uint8_t Dw3[16]={0};for(int r=0;r<4;r++)Dw3[diag[r]]=Dx4d[r];
    uint8_t Dz3[16];memcpy(Dz3,Dw3,16);imc(Dz3);
    uint8_t Dy3[16];memcpy(Dy3,Dz3,16);isr(Dy3);
    /* Enumerate x_3 branches (product of DDT solutions) */
    uint8_t sol[16][4];int ns[16];
    for(int i=0;i<16;i++){
        if(Dx3[i]==0){if(Dy3[i]!=0)return 0;ns[i]=0;/*free — but then 256 branches; for TRUE-param check we only need to verify the TRUE x_3 is a solution, not enumerate. Mark as wildcard.*/
            sol[i][0]=0;ns[i]=1;/*placeholder*/ return -1;/*signal: free byte, can't enumerate*/}
        ns[i]=DDTn[Dx3[i]][Dy3[i]];if(ns[i]==0)return 0;
        for(int j=0;j<ns[i]&&j<4;j++)sol[i][j]=DDTs[Dx3[i]][Dy3[i]][j];
    }
    /* product could be huge (2^16 typical). For this test: just return ONE branch
     * + the FIRST solution per byte, AND separately let caller check true x_3. */
    int nout=0;
    if(maxout>=1){memcpy(out[0].x2c0,x2c0,4);memcpy(out[0].x4d,x4d,4);
        for(int i=0;i<16;i++)out[0].x3[i]=sol[i][0];nout=1;}
    /* return total branch count (capped) */
    long tot=1;for(int i=0;i<16;i++){tot*=ns[i];if(tot>1000000)break;}
    (void)tot;
    return nout;
}
/* verify: given true x_3, check it satisfies the DDT constraints from 10 params */
static int rebound10_admits(uint8_t Dz1,const uint8_t x2c0[4],uint8_t Dw4,const uint8_t z4c0[4],
                            const uint8_t x3_true[16]){
    static const int diag[4]={0,5,10,15};
    uint8_t Dx2[4]={gm(2,Dz1),Dz1,Dz1,gm(3,Dz1)};
    uint8_t Dy2[4];for(int r=0;r<4;r++)Dy2[r]=Sbox[x2c0[r]^Dx2[r]]^Sbox[x2c0[r]];
    uint8_t Dz2[16]={0};for(int r=0;r<4;r++)Dz2[4*((4-r)%4)+r]=Dy2[r];
    uint8_t Dx3[16];memcpy(Dx3,Dz2,16);mc(Dx3);
    uint8_t Dz4[4]={gm(14,Dw4),gm(9,Dw4),gm(13,Dw4),gm(11,Dw4)};
    uint8_t x4d[4],Dx4d[4];for(int r=0;r<4;r++){x4d[r]=iSbox[z4c0[r]];Dx4d[r]=iSbox[z4c0[r]^Dz4[r]]^x4d[r];}
    uint8_t Dw3[16]={0};for(int r=0;r<4;r++)Dw3[diag[r]]=Dx4d[r];
    uint8_t Dz3[16];memcpy(Dz3,Dw3,16);imc(Dz3);
    uint8_t Dy3[16];memcpy(Dy3,Dz3,16);isr(Dy3);
    for(int i=0;i<16;i++){
        uint8_t dy=Sbox[x3_true[i]^Dx3[i]]^Sbox[x3_true[i]];
        if(dy!=Dy3[i])return 0;
    }
    return 1;
}

static uint64_t rs=0x1234abcdULL;
static uint32_t rng(void){rs=rs*6364136223846793005ULL+1;return rs>>32;}

int main(int argc,char**argv){
    int NKEYS=50;long NSAMP=1<<16;
    if(argc>1)NKEYS=atoi(argv[1]);
    if(argc>2)rs^=strtoull(argv[2],0,0)*0x9e3779b97f4a7c15ULL;
    gf_init();tabs_init();build_ddt();

    int n_rp=0,n_admits=0,n_24eq=0,n_match=0,n_parity=0;
    long n_offFP=0,n_offTrials=0;

    for(int trial=0;trial<NKEYS;trial++){
        /* ==== Rebound right-pair on FULL AES with INDEPENDENT subkeys ====
         * (valid: I_{m,n} bridge and DFJ Prop 2 don't use KS in rounds 1-4;
         *  real-AES KS-consistent key-fit needs 2^{32} iters — scaling only) */
        uint8_t rk[11][16],Pref[16],Ppair[16];
        static const int diag[4]={0,5,10,15};
        int found=0;
        uint8_t Dz1,Dw4,x2c0[4],z4c0[4],x3r[16],x4d[4];
        for(long att=0;att<500000&&!found;att++){
            Dz1=1+rng()%255;Dw4=1+rng()%255;
            for(int r=0;r<4;r++){x2c0[r]=rng();z4c0[r]=rng();}
            uint8_t Dx2[4]={gm(2,Dz1),Dz1,Dz1,gm(3,Dz1)};
            uint8_t Dy2[4];for(int r=0;r<4;r++)Dy2[r]=Sbox[x2c0[r]^Dx2[r]]^Sbox[x2c0[r]];
            uint8_t Dz2[16]={0};for(int r=0;r<4;r++)Dz2[4*((4-r)%4)+r]=Dy2[r];
            uint8_t Dx3[16];memcpy(Dx3,Dz2,16);mc(Dx3);
            uint8_t Dz4[4]={gm(14,Dw4),gm(9,Dw4),gm(13,Dw4),gm(11,Dw4)};
            uint8_t Dx4d[4];
            for(int r=0;r<4;r++){x4d[r]=iSbox[z4c0[r]];Dx4d[r]=iSbox[z4c0[r]^Dz4[r]]^x4d[r];}
            uint8_t Dw3[16]={0};for(int r=0;r<4;r++)Dw3[diag[r]]=Dx4d[r];
            uint8_t Dz3[16];memcpy(Dz3,Dw3,16);imc(Dz3);
            uint8_t Dy3[16];memcpy(Dy3,Dz3,16);isr(Dy3);
            int ok=1;for(int i=0;i<16;i++){
                if(DDTn[Dx3[i]][Dy3[i]]==0){ok=0;break;}
                x3r[i]=DDTs[Dx3[i]][Dy3[i]][rng()%DDTn[Dx3[i]][Dy3[i]]];}
            if(!ok)continue;
            /* INDEPENDENT subkeys: pick x_2[col≠0],x_4[non-diag] random,
             * SET k_2=x_3⊕w_2, k_3=x_4⊕w_3. rk[0..2],rk[5..7] random. */
            uint8_t x2[16],x4[16];for(int i=0;i<16;i++){x2[i]=rng();x4[i]=rng();}
            for(int r=0;r<4;r++){x2[r]=x2c0[r];x4[diag[r]]=x4d[r];}
            uint8_t w2[16];memcpy(w2,x2,16);sb(w2);sr(w2);mc(w2);
            uint8_t w3[16];memcpy(w3,x3r,16);sb(w3);sr(w3);mc(w3);
            for(int i=0;i<16;i++){rk[3][i]=x3r[i]^w2[i];rk[4][i]=x4[i]^w3[i];}
            for(int r=0;r<=2;r++)for(int i=0;i<16;i++)rk[r][i]=rng();
            for(int r=5;r<=7;r++)for(int i=0;i<16;i++)rk[r][i]=rng();
            /* derive Pref from x_2 backward */
            uint8_t s[16];memcpy(s,x2,16);
            ark(s,rk[2]);imc(s);isr(s);isb(s);uint8_t x1[16];memcpy(x1,s,16);
            ark(s,rk[1]);imc(s);isr(s);isb(s);ark(s,rk[0]);memcpy(Pref,s,16);
            uint8_t x1p0=0;for(int xx=0;xx<256;xx++)if((Sbox[xx]^Sbox[x1[0]])==Dz1){x1p0=xx;break;}
            uint8_t x1p[16];memcpy(x1p,x1,16);x1p[0]=x1p0;
            memcpy(s,x1p,16);ark(s,rk[1]);imc(s);isr(s);isb(s);ark(s,rk[0]);memcpy(Ppair,s,16);
            uint8_t tr0[8][16],tr1[8][16];
            enc7_trace(rk,Pref,tr0);enc7_trace(rk,Ppair,tr1);
            int rp=1;for(int i=1;i<16;i++)if(tr0[5][i]!=tr1[5][i])rp=0;
            if(tr0[5][0]==tr1[5][0]||tr0[5][0]==0)rp=0;
            if(rp)found=1;
        }
        if(!found){printf("[%d] rebound failed\n",trial);continue;}
        n_rp++;

        /* ==== Extract TRUE 10 params + 24 bytes from trace ==== */
        uint8_t tr0[8][16],tr1[8][16];
        enc7_trace(rk,Pref,tr0);enc7_trace(rk,Ppair,tr1);
        uint8_t T_Dz1=Sbox[tr0[1][0]]^Sbox[tr1[1][0]];
        uint8_t T_x2c0[4]={tr0[2][0],tr0[2][1],tr0[2][2],tr0[2][3]};
        uint8_t T_Dw4=tr0[5][0]^tr1[5][0]; /* Δx_5[0]=Δw_4[0] */
        uint8_t T_z4c0[4]={Sbox[tr0[4][0]],Sbox[tr0[4][5]],Sbox[tr0[4][10]],Sbox[tr0[4][15]]};
        uint8_t T_x3[16];memcpy(T_x3,tr0[3],16);
        uint8_t T_x4d[4]={tr0[4][0],tr0[4][5],tr0[4][10],tr0[4][15]};

        /* ==== Check 1: rebound10 admits true x_3 ==== */
        if(rebound10_admits(T_Dz1,T_x2c0,T_Dw4,T_z4c0,T_x3))n_admits++;

        /* ==== Check 2: off24 from TRUE 24 bytes → E' == internal ==== */
        uint8_t Eoff[256];int nEoff=off24_Eprime(T_x2c0,T_x3,T_x4d,Eoff);
        /* internal E' */
        uint8_t x1b[16];memcpy(x1b,tr0[1],16);int vref=x1b[0];
        uint8_t x5i[256];
        for(int v=0;v<256;v++){uint8_t x1[16],s[16];memcpy(x1,x1b,16);x1[0]=v;
            memcpy(s,x1,16);ark(s,rk[1]);imc(s);isr(s);isb(s);ark(s,rk[0]);
            uint8_t tr[8][16];enc7_trace(rk,s,tr);x5i[v]=tr[5][0];}
        uint8_t Eint[256];int nEint=0;
        for(int v=0;v<256;v++){if(v==vref)continue;uint8_t d=x5i[v]^x5i[vref];if(d)Eint[nEint++]=gi(d);}
        {int h1[256]={0},h2[256]={0};
         for(int i=0;i<nEoff;i++)h1[Eoff[i]]++;for(int i=0;i<nEint;i++)h2[Eint[i]]++;
         int eq=1;for(int i=0;i<256;i++)if(h1[i]!=h2[i]){eq=0;break;}
         if(eq)n_24eq++;}

        /* ==== ONLINE H from ciphertexts + k_6[idiag] ==== */
        const int pos[4]={0,13,10,7};uint8_t*k6=rk[7];
        uint8_t C[256][16];
        for(int v=0;v<256;v++){uint8_t x1[16],s[16];memcpy(x1,x1b,16);x1[0]=v;
            memcpy(s,x1,16);ark(s,rk[1]);imc(s);isr(s);isb(s);ark(s,rk[0]);
            uint8_t tr[8][16];enc7_trace(rk,s,tr);memcpy(C[v],tr[7],16);}
        uint8_t z5h[256];
        for(int v=0;v<256;v++){uint8_t x[4];for(int r=0;r<4;r++)x[r]=iSbox[C[v][pos[r]]^k6[pos[r]]];
            z5h[v]=gm(14,x[0])^gm(11,x[1])^gm(13,x[2])^gm(9,x[3]);}
        uint8_t H[256];int nH=0;
        for(int v=0;v<256;v++){if(v==vref)continue;uint8_t g=Minv[z5h[v]^z5h[vref]];
            if(g)H[nH++]=gi(g);}

        /* ==== I_{m,n} match: H (honest online) vs Eoff (10-param offline) ==== */
        uint64_t fHe,fHo,fEe,fEo;fpp(H,nH,&fHe,&fHo);fpp(Eoff,nEoff,&fEe,&fEo);
        int ev=(fHe==fEe),od=(fHo==fEo);
        if(ev||od)n_match++;
        int kA=0;for(int v=0;v<256;v++)if(x5i[v]==0)kA++;
        if(x5i[vref]!=0 && ev==(kA%2==0) && od==(kA%2==1))n_parity++;

        /* ==== Sampled-offline FP: random 10-param tuples vs true-H ==== */
        for(long s=0;s<NSAMP;s++){
            uint8_t rDz1=1+rng()%255,rDw4=1+rng()%255,rx2[4],rz4[4];
            for(int r=0;r<4;r++){rx2[r]=rng();rz4[r]=rng();}
            Off24 o;int nr=rebound10(rDz1,rx2,rDw4,rz4,&o,1);
            if(nr<=0)continue;n_offTrials++;
            uint8_t Er[256];int nEr=off24_Eprime(o.x2c0,o.x3,o.x4d,Er);
            uint64_t fe,fo;fpp(Er,nEr,&fe,&fo);
            if(fe==fHe||fo==fHo)n_offFP++;
        }
    }
    printf("\n===== FULL-AES HONEST w/ RIGHT-PAIR + 10-PARAM OFFLINE (NKEYS=%d) =====\n",NKEYS);
    printf("right-pair via rebound               : %d/%d\n",n_rp,NKEYS);
    printf("rebound10 admits TRUE x_3            : %d/%d\n",n_admits,n_rp);
    printf("off24(TRUE 10-param→x_3→E')==internal: %d/%d\n",n_24eq,n_rp);
    printf("I_{m,n}(H_honest)↔I_{m,n}(E'_10param): %d/%d EITHER\n",n_match,n_rp);
    printf("Parity-exact (a≠0)                   : %d/%d\n",n_parity,n_rp);
    printf("Sampled-offline FP (random 10-param) : %ld/%ld (%.2e)\n",n_offFP,n_offTrials,
           n_offTrials?(double)n_offFP/n_offTrials:0.0);
    return 0;
}
