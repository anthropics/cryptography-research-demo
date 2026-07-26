// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* saes_imn.c — Small-AES GF(2^4) end-to-end key recovery via I_{m,n} bridge.
 *
 * Pipeline (per random 64-bit key):
 *   1. Find right pair (known-key shortcut, ~2^12 trials).
 *   2. Build δ-set on x_1[0], encrypt → 16 ciphertexts.
 *   3. OFFLINE E'={1/Δx_5[0]} from internal (ground truth); fp_E_ev/od.
 *   4. ONLINE: enumerate K_on=(k_{-1}[diag],k_6[idiag]) via pair-DDT.
 *      For each: H={1/M⁻¹(Δz_5[0])} from ciphertexts; fp_H_ev/od.
 *      Match EITHER vs fp_E.
 *   5. Assert true K_on survives. Count survivors.
 *   6. Key completion: brute remaining k_6 bytes, KS⁻¹, trial-encrypt.
 *
 * Also: honest-online test (50 keys, no attack, just match check).
 *
 * Build:  gcc -O3 -march=native -o saes_imn saes_imn.c
 */
#include "saes.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static uint64_t rng_state=0x243f6a8885a308d3ULL;
static uint32_t rng(void){ rng_state=rng_state*6364136223846793005ULL+1; return rng_state>>32; }
static uint8_t rn(void){ return rng()&15; }

/* ---------- I_{m,n} fingerprint over GF(16) ---------- */
#define NPM 6
static const int Pm_exp[NPM]={3,5,7,6,9,11}; /* popcount≥2; includes Frobenius dups for robustness */

static void compute_Pm(const uint8_t*S,int n,uint8_t Pm[NPM]){
    memset(Pm,0,NPM);
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
        uint8_t d=S[i]^S[j]; if(!d)continue;
        for(int k=0;k<NPM;k++) Pm[k]^=sgf_pow(d,Pm_exp[k]);
    }
}
static uint32_t fingerprint(const uint8_t Pm[NPM]){
    int n0=-1; for(int k=0;k<NPM;k++)if(Pm[k]){n0=k;break;}
    if(n0<0) return 0;
    uint32_t fp=1; int n=Pm_exp[n0]; uint8_t Pn=Pm[n0];
    for(int k=0;k<NPM;k++){ if(k==n0)continue;
        uint8_t I=sgf_mul(sgf_pow(Pm[k],n),sgf_pow(sgf_inv(Pn),Pm_exp[k]));
        fp=fp*17+I;
    }
    fp=fp*17+(n0+1);
    return fp;
}
static void fp_pair(const uint8_t*S,int n,uint32_t*fev,uint32_t*fod){
    uint8_t Pm[NPM]; compute_Pm(S,n,Pm); *fev=fingerprint(Pm);
    uint8_t Pm2[NPM]; memcpy(Pm2,Pm,NPM);
    for(int i=0;i<n;i++){uint8_t s=S[i]; if(!s)continue;
        for(int k=0;k<NPM;k++)Pm2[k]^=sgf_pow(s,Pm_exp[k]);}
    *fod=fingerprint(Pm2);
}

/* ---------- DDT for S-box ---------- */
static uint8_t DDT_sol[16][16][4]; /* DDT_sol[din][dout][i]=input values x with SB(x)^SB(x^din)=dout */
static uint8_t DDT_n[16][16];
static void build_ddt(void){
    memset(DDT_n,0,sizeof(DDT_n));
    for(int din=0;din<16;din++)for(int x=0;x<16;x++){
        int dout=sSB[x]^sSB[x^din];
        if(DDT_n[din][dout]<4) DDT_sol[din][dout][DDT_n[din][dout]]=x;
        DDT_n[din][dout]++;
    }
}

/* ---------- partial decrypt: Δz_5[rr] from C + k_6[4 bytes] ---------- */
/* z_5[rr,0] = (iMC row rr)·x_6[·,0] ⊕ u_5[rr,0].  x_6[r,0]=iSB(C[pos_r]^k_6[pos_r]).
 * But for rr≠0, z_5[rr,0]=y_5[rr,rr] after SR... wait, y_5[r,0]=z_5[r,(−r)%4].
 * Actually we want z_5[0,0] (match byte). For cascade rr=0..3 we observe
 * z_5[rr]:=y_5[rr,0]=z_5[rr,(-rr)%4], which needs x_6[·,(-rr)%4] via MC⁻¹...
 * For simplicity in this PoC, cascade will observe Δx_5[rr,0] directly
 * (needs u_5[rr,0] — brute-force 16 each in completion phase).
 *
 * Primary match: rr=0, Δz_5[0,0] via k_6[0,13,10,7]. */
static const int pos0[4]={0,13,10,7};
static const int posR[4][4]={{0,13,10,7},{12,9,6,3},{8,5,2,15},{4,1,14,11}};
static const uint8_t imcR[4][4]={{0xE,0xB,0xD,0x9},{0x9,0xE,0xB,0xD},
                                 {0xD,0x9,0xE,0xB},{0xB,0xD,0x9,0xE}};
static inline uint8_t pdec_z5(int rr,const uint8_t C[16],const uint8_t k6_4[4]){
    uint8_t x0=siSB[C[posR[rr][0]]^k6_4[0]],x1=siSB[C[posR[rr][1]]^k6_4[1]],
            x2=siSB[C[posR[rr][2]]^k6_4[2]],x3=siSB[C[posR[rr][3]]^k6_4[3]];
    return sgf_mul(imcR[rr][0],x0)^sgf_mul(imcR[rr][1],x1)^
           sgf_mul(imcR[rr][2],x2)^sgf_mul(imcR[rr][3],x3);
}
static inline uint8_t pdec_z5_0(const uint8_t C[16], const uint8_t k6_4[4]){
    return pdec_z5(0,C,k6_4);
}
/* multiset equality (as histograms) */
static int mset_eq(const uint8_t*A,int nA,const uint8_t*B,int nB){
    if(nA!=nB)return 0; int h[16]={0};
    for(int i=0;i<nA;i++)h[A[i]]++; for(int i=0;i<nB;i++)h[B[i]]--;
    for(int i=0;i<16;i++)if(h[i])return 0; return 1;
}

/* ---------- self-test ---------- */
static void selftest(void){
    /* S-box bijective */
    int seen[16]={0}; for(int x=0;x<16;x++)seen[sSB[x]]++;
    for(int i=0;i<16;i++)assert(seen[i]==1);
    /* Minv linear */
    assert(sMinv[0]==0);
    for(int a=0;a<16;a++)for(int b=0;b<16;b++) assert(sMinv[a^b]==(sMinv[a]^sMinv[b]));
    /* MC∘iMC=I */
    for(int t=0;t<100;t++){uint8_t c[4]={rn(),rn(),rn(),rn()},c2[4];memcpy(c2,c,4);
        smc_col(c2);simc_col(c2);assert(memcmp(c,c2,4)==0);}
    /* enc/dec roundtrip */
    for(int t=0;t<20;t++){
        uint8_t key[16],rk[11][16],P[16],tr[8][16];
        for(int i=0;i<16;i++){key[i]=rn();P[i]=rn();}
        skey_expand(key,rk); senc7_trace(rk,P,tr);
        /* decrypt */
        uint8_t s[16]; memcpy(s,tr[7],16); sark(s,rk[7]); sisr(s); sisb(s);
        assert(memcmp(s,tr[6],16)==0);
        for(int r=5;r>=0;r--){ sark(s,rk[r+1]); simc(s); sisr(s); sisb(s);
            assert(memcmp(s,tr[r],16)==0); }
        sark(s,rk[0]); assert(memcmp(s,P,16)==0);
        /* KS invert */
        uint8_t rk2[11][16]; skey_invert_from7(rk[7],rk2);
        assert(memcmp(rk2[0],key,16)==0);
    }
    /* bridge h=a⊕a²d */
    int ok=0,tot=0;
    for(int a=1;a<16;a++)for(int x=1;x<16;x++){if(x==a)continue;
        uint8_t dz=sSB[x]^sSB[a], g=sMinv[dz], h=sgf_inv(g);
        uint8_t d=sgf_inv(x^a), pred=a^sgf_mul(sgf_mul(a,a),d);
        tot++; if(h==pred)ok++;
    }
    assert(ok==tot);
    printf("[selftest] S-box bij, Minv linear, MC/iMC, enc/dec, KS-inv, bridge %d/%d OK\n",ok,tot);
}

/* ========================================================================= */
/*                         PART A: honest-online test                        */
/* ========================================================================= */
static void honest_test(int NKEYS){
    int n_pdec=0,n_ev=0,n_od=0,n_either=0,n_a0=0,n_par=0,n_fp=0,n_fpt=0;
    for(int trial=0;trial<NKEYS;trial++){
        uint8_t key[16],rk[11][16]; for(int i=0;i<16;i++)key[i]=rn();
        skey_expand(key,rk);
        const uint8_t *km1=rk[0],*k0=rk[1],*k6=rk[7];
        uint8_t x1b[16]; for(int i=0;i<16;i++)x1b[i]=rn();
        uint8_t P[16][16],C[16][16],x5_0[16];
        for(int v=0;v<16;v++){
            uint8_t x1[16];memcpy(x1,x1b,16);x1[0]=v;
            uint8_t s[16];memcpy(s,x1,16);sark(s,k0);simc(s);sisr(s);sisb(s);sark(s,km1);
            memcpy(P[v],s,16);
            uint8_t tr[8][16];senc7_trace(rk,P[v],tr);
            memcpy(C[v],tr[7],16);x5_0[v]=tr[5][0];
            assert(memcmp(tr[1],x1,16)==0);
        }
        int vref=x1b[0]; uint8_t a=x5_0[vref];
        /* honest pdec */
        uint8_t k6_4[4]={k6[0],k6[13],k6[10],k6[7]};
        uint8_t z5h[16]; int pdec=1;
        for(int v=0;v<16;v++) z5h[v]=pdec_z5_0(C[v],k6_4);
        for(int v=0;v<16;v++)
            if((z5h[v]^z5h[vref])!=(sSB[x5_0[v]]^sSB[x5_0[vref]]))pdec=0;
        if(pdec)n_pdec++;
        /* H, E' */
        uint8_t H[16],E[16];int nH=0,nE=0;
        for(int v=0;v<16;v++){if(v==vref)continue;
            uint8_t g=sMinv[z5h[v]^z5h[vref]]; if(g)H[nH++]=sgf_inv(g);
            uint8_t dx=x5_0[v]^x5_0[vref]; if(dx)E[nE++]=sgf_inv(dx);}
        int kA=0;for(int v=0;v<16;v++)if(x5_0[v]==0)kA++;
        uint32_t fHe,fHo,fEe,fEo; fp_pair(H,nH,&fHe,&fHo); fp_pair(E,nE,&fEe,&fEo);
        int ev=(fHe==fEe),od=(fHo==fEo);
        if(a==0)n_a0++;
        if(ev)n_ev++; if(od)n_od++; if(ev||od)n_either++;
        if(a!=0){ if(ev==(kA%2==0)&&od==(kA%2==1))n_par++; }
        /* wrong-k6 FP */
        for(int w=0;w<20;w++){
            uint8_t wk[4]={rn(),rn(),rn(),rn()}; if(memcmp(wk,k6_4,4)==0)wk[0]^=1;
            uint8_t Hw[16];int nw=0;
            uint8_t zr=pdec_z5_0(C[vref],wk);
            for(int v=0;v<16;v++){if(v==vref)continue;
                uint8_t g=sMinv[pdec_z5_0(C[v],wk)^zr]; if(g)Hw[nw++]=sgf_inv(g);}
            uint32_t fe,fo;fp_pair(Hw,nw,&fe,&fo);n_fpt++;
            if(fe==fEe||fo==fEo)n_fp++;
        }
    }
    printf("\n===== SMALL-AES HONEST-ONLINE (NKEYS=%d) =====\n",NKEYS);
    printf("Δz_5[0] pdec==internal : %d/%d\n",n_pdec,NKEYS);
    printf("even-match             : %d/%d\n",n_ev,NKEYS);
    printf("odd-match (∪{0})       : %d/%d\n",n_od,NKEYS);
    printf("EITHER                 : %d/%d\n",n_either,NKEYS);
    printf("  a=0 edge             : %d\n",n_a0);
    printf("Parity exact (a≠0)     : %d/%d\n",n_par,NKEYS-n_a0);
    printf("wrong-k6 FP            : %d/%d (%.4f)\n",n_fp,n_fpt,(double)n_fp/n_fpt);
}

/* ========================================================================= */
/*                    PART B: end-to-end key recovery                        */
/* ========================================================================= */

/* Enumerate k_{-1}[diag] candidates from pair (P,P'): need Δw_0[1,2,3]=0,
 * i.e. Δz_0[col0]=δ·[E,9,D,B] for some δ. Δz_0[r,0]=Δy_0[r,r]=SB(P[d_r]^k)^SB(P'[d_r]^k).
 * For each δ∈GF(16)*, solve 4 DDT equations independently. */
static int enum_km1(const uint8_t P[16],const uint8_t Pp[16], uint8_t out[][4]){
    static const int diag[4]={0,5,10,15};
    static const uint8_t coef[4]={0xE,0x9,0xD,0xB};
    int n=0;
    for(int delta=1;delta<16;delta++){
        uint8_t cand[4][4]; int nc[4];
        int ok=1;
        for(int r=0;r<4;r++){
            uint8_t din=P[diag[r]]^Pp[diag[r]], dout=sgf_mul(coef[r],delta);
            nc[r]=DDT_n[din][dout]; if(nc[r]==0){ok=0;break;}
            for(int i=0;i<nc[r]&&i<4;i++) cand[r][i]=DDT_sol[din][dout][i]^P[diag[r]];
        }
        if(!ok)continue;
        for(int i0=0;i0<nc[0];i0++)for(int i1=0;i1<nc[1];i1++)
        for(int i2=0;i2<nc[2];i2++)for(int i3=0;i3<nc[3];i3++){
            out[n][0]=cand[0][i0];out[n][1]=cand[1][i1];
            out[n][2]=cand[2][i2];out[n][3]=cand[3][i3];n++;
        }
    }
    return n;
}
/* Enumerate k_6[idiag] candidates from pair (C,C'): need Δz_5[1,2,3]=0 ⟹
 * Δx_6[col0]=MC[:,0]·β=[2β,β,β,3β]. Δx_6[r,0]=iSB(C[pos_r]^k)^iSB(C'[pos_r]^k).
 * Use DDT of iSB. */
static int enum_k6(const uint8_t C[16],const uint8_t Cp[16], uint8_t out[][4]){
    static const uint8_t coef[4]={2,1,1,3};
    int n=0;
    for(int beta=1;beta<16;beta++){
        uint8_t cand[4][4]; int nc[4]; int ok=1;
        for(int r=0;r<4;r++){
            /* iSB DDT: find y s.t. iSB(y)^iSB(y^din)=dout ⟺ SB DDT with roles swapped */
            uint8_t din=C[pos0[r]]^Cp[pos0[r]], dout=sgf_mul(coef[r],beta);
            /* y s.t. iSB(y)^iSB(y^din)=dout ⟺ x=iSB(y), SB(x)^SB(x^dout)=din */
            int m=DDT_n[dout][din]; if(m==0){ok=0;break;}
            nc[r]=m;
            for(int i=0;i<m&&i<4;i++) cand[r][i]=sSB[DDT_sol[dout][din][i]]^C[pos0[r]];
        }
        if(!ok)continue;
        for(int i0=0;i0<nc[0];i0++)for(int i1=0;i1<nc[1];i1++)
        for(int i2=0;i2<nc[2];i2++)for(int i3=0;i3<nc[3];i3++){
            out[n][0]=cand[0][i0];out[n][1]=cand[1][i1];
            out[n][2]=cand[2][i2];out[n][3]=cand[3][i3];n++;
        }
    }
    return n;
}

/* Rebound right-pair construction: returns a (key, P, P') with
 * ΔP[diag only], Δx_1=[δ1,0..], Δx_5=[δ5,0..]. Known-key shortcut. */
static long rb_ddt_ok=0, rb_keyfit_tries=0;
static int rebound_rightpair(uint8_t key_out[16], uint8_t P_out[16], uint8_t Pp_out[16]){
    static const int diag[4]={0,5,10,15};
    for(long attempt=0;attempt<10000000;attempt++){
        /* 10 params */
        uint8_t Dz1=1+rn()%15, Dw4=1+rn()%15;
        uint8_t x2c0[4]={rn(),rn(),rn(),rn()}, z4c0[4]={rn(),rn(),rn(),rn()};
        /* forward → Δx_3 */
        uint8_t Dx2c0[4]={sgf_mul(2,Dz1),sgf_mul(1,Dz1),sgf_mul(1,Dz1),sgf_mul(3,Dz1)};
        uint8_t Dy2c0[4]; for(int r=0;r<4;r++)Dy2c0[r]=sSB[x2c0[r]^Dx2c0[r]]^sSB[x2c0[r]];
        uint8_t Dz2[16]={0}; /* SR: z_2[r,(-r)%4... no: z[r,c]=y[r,(c+r)%4]. y_2[r,0] → z_2[r,(4-r)%4]? */
        /* y_2 has col0 active. z_2[r,c]=y_2[r,(c+r)%4]. Active when (c+r)%4==0 ⟹ c=(4-r)%4. */
        for(int r=0;r<4;r++) Dz2[4*((4-r)%4)+r]=Dy2c0[r];
        uint8_t Dx3[16]; memcpy(Dx3,Dz2,16); smc(Dx3); /* Δw_2=Δx_3 */
        /* backward → Δy_3 */
        uint8_t Dz4c0[4]={sgf_mul(0xE,Dw4),sgf_mul(0x9,Dw4),sgf_mul(0xD,Dw4),sgf_mul(0xB,Dw4)};
        uint8_t x4d[4],Dx4d[4];
        for(int r=0;r<4;r++){
            uint8_t y=z4c0[r],Dy=Dz4c0[r]; /* y_4[r,r]=z_4[r,0] */
            x4d[r]=siSB[y]; Dx4d[r]=siSB[y^Dy]^x4d[r];
        }
        uint8_t Dw3[16]={0}; for(int r=0;r<4;r++)Dw3[diag[r]]=Dx4d[r];
        uint8_t Dz3[16]; memcpy(Dz3,Dw3,16); simc(Dz3);
        uint8_t Dy3[16]; memcpy(Dy3,Dz3,16); sisr(Dy3);
        /* DDT rebound: x_3[i] s.t. SB(x)⊕SB(x⊕Dx3[i])=Dy3[i] */
        uint8_t x3[16]; int ok=1;
        for(int i=0;i<16;i++){
            if(Dx3[i]==0&&Dy3[i]==0){x3[i]=rn();continue;}
            if(DDT_n[Dx3[i]][Dy3[i]]==0){ok=0;break;}
            x3[i]=DDT_sol[Dx3[i]][Dy3[i]][rng()%DDT_n[Dx3[i]][Dy3[i]]];
        }
        if(!ok)continue;
        rb_ddt_ok++;
        /* Now loop x_2[col≠0] until k_2→k_3 gives x_4[diag]==target */
        for(long it=0;it<(1L<<18);it++){
            rb_keyfit_tries++;
            uint8_t x2[16]; for(int i=0;i<16;i++)x2[i]=rn();
            for(int r=0;r<4;r++)x2[r]=x2c0[r];
            uint8_t w2[16]; memcpy(w2,x2,16); ssb(w2);ssr(w2);smc(w2);
            uint8_t k2[16]; for(int i=0;i<16;i++)k2[i]=x3[i]^w2[i];
            /* k_2=rk[3]. Need rk[4]=k_3. Quick KS step: */
            uint8_t k3[16];
            {uint8_t t[4]={k2[12],k2[13],k2[14],k2[15]};
             uint8_t tt=t[0];t[0]=sSB[t[1]];t[1]=sSB[t[2]];t[2]=sSB[t[3]];t[3]=sSB[tt];
             t[0]^=sRcon[4];
             for(int i=0;i<4;i++)k3[i]=k2[i]^t[i];
             for(int c=1;c<4;c++)for(int i=0;i<4;i++)k3[4*c+i]=k2[4*c+i]^k3[4*(c-1)+i];}
            uint8_t w3[16]; memcpy(w3,x3,16); ssb(w3);ssr(w3);smc(w3);
            int m=1; for(int r=0;r<4;r++)if((w3[diag[r]]^k3[diag[r]])!=x4d[r]){m=0;break;}
            if(!m)continue;
            uint8_t rk[11][16]; skey_from(3,k2,rk);
            /* FOUND. Derive P from x_2 backward. */
            uint8_t s[16]; memcpy(s,x2,16);
            sark(s,rk[2]);simc(s);sisr(s);sisb(s); /* x_1 */
            uint8_t x1[16]; memcpy(x1,s,16);
            sark(s,rk[1]);simc(s);sisr(s);sisb(s);sark(s,rk[0]); /* P */
            memcpy(P_out,s,16); memcpy(key_out,rk[0],16);
            /* verify */
            uint8_t tr[8][16]; senc7_trace(rk,P_out,tr);
            if(memcmp(tr[2],x2,16)||memcmp(tr[3],x3,16)) continue; /* shouldn't happen */
            /* P' = P with Δx_1=[δ1,0..] where Δy_1[0]=Dz1 */
            /* Need x_1[0]' s.t. SB(x_1[0]')⊕SB(x_1[0])=Dz1 */
            int nf=DDT_n[0][0]; (void)nf;
            uint8_t x1p0=0; int f=0;
            for(int xx=0;xx<16;xx++)if((sSB[xx]^sSB[x1[0]])==Dz1){x1p0=xx;f=1;break;}
            if(!f)continue;
            uint8_t x1p[16]; memcpy(x1p,x1,16); x1p[0]=x1p0;
            memcpy(s,x1p,16);sark(s,rk[1]);simc(s);sisr(s);sisb(s);sark(s,rk[0]);
            memcpy(Pp_out,s,16);
            /* final check: Δx_5[1..15]=0 */
            uint8_t tr2[8][16]; senc7_trace(rk,Pp_out,tr2);
            int rp=1; for(int i=1;i<16;i++)if(tr[5][i]!=tr2[5][i])rp=0;
            if(tr[5][0]==tr2[5][0])rp=0;
            if(tr[5][0]==0)continue; /* avoid a=0 */
            if(rp) return 1;
        }
    }
    return 0;
}

static int e2e_attack(int NKEYS, int verbose){
    int n_rp=0,n_true_surv=0,n_key_ok=0;
    long tot_surv=0,tot_kon=0;
    for(int trial=0;trial<NKEYS;trial++){
        uint8_t key[16],rk[11][16],Pref[16],Ppair[16],Cref[16],Cpair[16];

        /* ---- 1. rebound right-pair (generates key+pair together) ---- */
        if(!rebound_rightpair(key,Pref,Ppair)){
            if(verbose)printf("[%d] rebound failed\n",trial);continue;}
        skey_expand(key,rk);
        const uint8_t *km1=rk[0],*k0=rk[1],*k6=rk[7];
        uint8_t tr0[8][16],tr1[8][16];
        senc7_trace(rk,Pref,tr0);senc7_trace(rk,Ppair,tr1);
        memcpy(Cref,tr0[7],16);memcpy(Cpair,tr1[7],16);
        n_rp++;

        /* build δ-set on x_1[0] around P_ref, record x_5[0] */
        uint8_t x1b[16]; memcpy(x1b,tr0[1],16);
        int vref=x1b[0]; uint8_t x5_0[16];
        for(int v=0;v<16;v++){
            uint8_t x1[16],s[16];memcpy(x1,x1b,16);x1[0]=v;
            memcpy(s,x1,16);sark(s,k0);simc(s);sisr(s);sisb(s);sark(s,km1);
            uint8_t tr[8][16];senc7_trace(rk,s,tr);x5_0[v]=tr[5][0];
        }
        uint8_t a=x5_0[vref];

        /* ---- 2. OFFLINE E' (ground truth) ---- */
        uint8_t E[16];int nE=0;
        for(int v=0;v<16;v++){if(v==vref)continue;
            uint8_t dx=x5_0[v]^a; if(dx)E[nE++]=sgf_inv(dx);}
        uint32_t fEe,fEo; fp_pair(E,nE,&fEe,&fEo);

        /* ---- 3. ONLINE: enumerate K_on via pair-DDT ---- */
        static uint8_t km1c[4096][4],k6c[4096][4];
        int nkm1=enum_km1(Pref,Ppair,km1c);
        int nk6 =enum_k6(Cref,Cpair,k6c);
        tot_kon+=nkm1*nk6;

        /* find true K_on indices */
        int tkm1=-1,tk6=-1;
        uint8_t km1t[4]={km1[0],km1[5],km1[10],km1[15]};
        uint8_t k6t[4]={k6[0],k6[13],k6[10],k6[7]};
        for(int i=0;i<nkm1;i++)if(memcmp(km1c[i],km1t,4)==0)tkm1=i;
        for(int i=0;i<nk6;i++) if(memcmp(k6c[i],k6t,4)==0)tk6=i;
        assert(tkm1>=0&&tk6>=0);

        /* For each K_on cand: rebuild δ-set, query oracle (=encrypt), compute H, match. */
        static const int diag[4]={0,5,10,15};
        static const uint8_t coef[4]={0xE,0x9,0xD,0xB};
        int surv=0,true_surv=0;
        typedef struct{int ikm1,ik6;}Surv;
        static Surv survs[1<<16];
        for(int ikm1=0;ikm1<nkm1;ikm1++){
            /* build δ-set plaintexts from Pref + km1c[ikm1] */
            uint8_t y0r[4];
            for(int r=0;r<4;r++)y0r[r]=sSB[Pref[diag[r]]^km1c[ikm1][r]];
            uint8_t Cd[16][16];
            for(int t=0;t<16;t++){
                uint8_t Pt[16];memcpy(Pt,Pref,16);
                for(int r=0;r<4;r++)Pt[diag[r]]=siSB[y0r[r]^sgf_mul(coef[r],t)]^km1c[ikm1][r];
                senc7(rk,Pt,Cd[t]); /* oracle query */
            }
            for(int ik6=0;ik6<nk6;ik6++){
                uint8_t z5[16];
                for(int t=0;t<16;t++)z5[t]=pdec_z5_0(Cd[t],k6c[ik6]);
                /* ref = t=0 = Pref */
                uint8_t H[16];int nH=0;
                for(int t=1;t<16;t++){uint8_t g=sMinv[z5[t]^z5[0]];if(g)H[nH++]=sgf_inv(g);}
                uint32_t fHe,fHo;fp_pair(H,nH,&fHe,&fHo);
                if(fHe==fEe||fHo==fEo){
                    if(surv<(1<<16)){survs[surv].ikm1=ikm1;survs[surv].ik6=ik6;}
                    surv++;
                    if(ikm1==tkm1&&ik6==tk6)true_surv=1;
                }
            }
        }
        tot_surv+=surv;
        if(true_surv)n_true_surv++;

        /* ---- 4. SECONDARY filter (multiset + λ-recover) on I_{m,n} survivors.
         * Small-field aid: I_{m,n}≈12 bits on GF(16) vs ~88 on GF(256). ---- */
        /* Ground-truth Δx_5[r,0] over δ-set (= offline match values) */
        uint8_t Dx5[4][16]; /* Dx5[r][v] = x_5[r,0]_v ⊕ x_5[r,0]_ref (from internal) */
        {   uint8_t x5r[4][16];
            for(int v=0;v<16;v++){
                uint8_t x1[16],s[16];memcpy(x1,x1b,16);x1[0]=v;
                memcpy(s,x1,16);sark(s,k0);simc(s);sisr(s);sisb(s);sark(s,km1);
                uint8_t tr[8][16];senc7_trace(rk,s,tr);
                for(int r=0;r<4;r++)x5r[r][v]=tr[5][r];
            }
            for(int r=0;r<4;r++)for(int v=0;v<16;v++)Dx5[r][v]=x5r[r][v]^x5r[r][vref];
        }
        int surv2=0,true_surv2=0;
        typedef struct{uint8_t km1[4],k6p[4][4];}FullCand;
        static FullCand fc[65536]; int nfc=0;
        for(int s=0;s<surv&&s<(1<<16);s++){
            int ikm1=survs[s].ikm1,ik6=survs[s].ik6;
            /* rebuild δ-set for this km1 cand, get ciphertexts */
            uint8_t y0r[4],Cd[16][16];
            for(int r=0;r<4;r++)y0r[r]=sSB[Pref[diag[r]]^km1c[ikm1][r]];
            for(int t=0;t<16;t++){
                uint8_t Pt[16];memcpy(Pt,Pref,16);
                for(int r=0;r<4;r++)Pt[diag[r]]=siSB[y0r[r]^sgf_mul(coef[r],t)]^km1c[ikm1][r];
                senc7(rk,Pt,Cd[t]);
            }
            /* rr=0 secondary: Dz5[0]_t vs predicted from Dx5[0]+a */
            uint8_t Dz[16]; uint8_t zref=pdec_z5(0,Cd[0],k6c[ik6]);
            for(int t=0;t<16;t++)Dz[t]=pdec_z5(0,Cd[t],k6c[ik6])^zref;
            int pass=0;
            for(int aa=1;aa<16;aa++){
                uint8_t Pred[16]; for(int v=0;v<16;v++)Pred[v]=sSB[Dx5[0][v]^aa]^sSB[aa];
                /* compare as multisets over t vs v (reindexed) */
                if(mset_eq(Dz,16,Pred,16)){pass=1;break;}
            }
            if(!pass)continue;
            surv2++;
            if(ikm1==tkm1&&ik6==tk6)true_surv2=1;
            /* ---- 5. Cascade rr=1,2,3: brute k_6[posR[rr]] (2^16), mset+λ match ---- */
            uint8_t k6pR[4][256][4]; int nk6pR[4]={0};
            memcpy(k6pR[0][0],k6c[ik6],4);nk6pR[0]=1;
            int cascade_ok=1;
            for(int rr=1;rr<4;rr++){
                for(int g=0;g<65536;g++){
                    uint8_t k4[4]={g&15,(g>>4)&15,(g>>8)&15,(g>>12)&15};
                    uint8_t Dzr[16],zr=pdec_z5(rr,Cd[0],k4);
                    for(int t=0;t<16;t++)Dzr[t]=pdec_z5(rr,Cd[t],k4)^zr;
                    for(int aa=0;aa<16;aa++){ /* include aa=0 for edge */
                        uint8_t Pr[16];for(int v=0;v<16;v++)Pr[v]=sSB[Dx5[rr][v]^aa]^sSB[aa];
                        if(mset_eq(Dzr,16,Pr,16)){
                            if(nk6pR[rr]<256)memcpy(k6pR[rr][nk6pR[rr]],k4,4);
                            nk6pR[rr]++;break;
                        }
                    }
                }
                if(nk6pR[rr]==0||nk6pR[rr]>=256)cascade_ok=0;
                /* sanity: true k_6[posR[rr]] must be among survivors when km1,k6[rr=0] true */
                if(ikm1==tkm1&&ik6==tk6){
                    uint8_t tk[4]={k6[posR[rr][0]],k6[posR[rr][1]],k6[posR[rr][2]],k6[posR[rr][3]]};
                    int f=0;for(int i=0;i<nk6pR[rr]&&i<256;i++)if(memcmp(k6pR[rr][i],tk,4)==0)f=1;
                    if(!f&&verbose)printf("    [rr=%d] TRUE k_6[pos] NOT in cascade survivors! nk6pR=%d\n",rr,nk6pR[rr]);
                }
            }
            if(!cascade_ok&&verbose)printf("    cascade overflow/empty: %d %d %d\n",nk6pR[1],nk6pR[2],nk6pR[3]);
            /* ---- 6. Assemble full k_6, KS⁻¹, trial-enc ---- */
            for(int i1=0;i1<nk6pR[1]&&i1<256;i1++)
            for(int i2=0;i2<nk6pR[2]&&i2<256;i2++)
            for(int i3=0;i3<nk6pR[3]&&i3<256;i3++){
                uint8_t k6full[16];
                for(int r=0;r<4;r++){
                    k6full[posR[0][r]]=k6pR[0][0][r];
                    k6full[posR[1][r]]=k6pR[1][i1][r];
                    k6full[posR[2][r]]=k6pR[2][i2][r];
                    k6full[posR[3][r]]=k6pR[3][i3][r];
                }
                uint8_t rk2[11][16]; skey_from(7,k6full,rk2);
                /* check k_{-1}[diag] consistency */
                if(rk2[0][0]!=km1c[ikm1][0]||rk2[0][5]!=km1c[ikm1][1]||
                   rk2[0][10]!=km1c[ikm1][2]||rk2[0][15]!=km1c[ikm1][3])continue;
                /* trial-enc */
                uint8_t Ct[16];senc7(rk2,Pref,Ct);
                if(memcmp(Ct,Cref,16)==0){
                    if(nfc<65536){memcpy(fc[nfc].km1,km1c[ikm1],4);
                        for(int rr=0;rr<4;rr++)memcpy(fc[nfc].k6p[rr],(rr==0?k6pR[0][0]:rr==1?k6pR[1][i1]:rr==2?k6pR[2][i2]:k6pR[3][i3]),4);}
                    nfc++;
                }
            }
        }
        int key_ok=(nfc>=1);
        /* verify recovered key(s) include true */
        int true_in_fc=0;
        for(int i=0;i<nfc;i++){
            uint8_t k6f[16];for(int rr=0;rr<4;rr++)for(int r=0;r<4;r++)k6f[posR[rr][r]]=fc[i].k6p[rr][r];
            if(memcmp(k6f,k6,16)==0)true_in_fc=1;
        }
        if(key_ok&&true_in_fc)n_key_ok++;

        if(verbose)printf("[%d] Kon=%d surv(Imn)=%d surv(mset)=%d true∈Imn=%d true∈mset=%d keys=%d true_key∈=%d\n",
                          trial,nkm1*nk6,surv,surv2,true_surv,true_surv2,nfc,true_in_fc);
    }
    printf("\n===== SMALL-AES END-TO-END (NKEYS=%d) =====\n",NKEYS);
    printf("right-pair found                      : %d/%d\n",n_rp,NKEYS);
    printf("true K_on survives I_{m,n} filter     : %d/%d\n",n_true_surv,n_rp);
    printf("mean |K_on cands| (pair-DDT)          : %.1f\n",(double)tot_kon/n_rp);
    printf("mean |I_{m,n} survivors|              : %.2f\n",(double)tot_surv/n_rp);
    printf("I_{m,n} filter ratio                  : %.4f (small-field artifact)\n",(double)tot_surv/tot_kon);
    printf("FULL 64-bit KEY RECOVERED (trial-enc) : %d/%d\n",n_key_ok,n_rp);
    return n_key_ok;
}

int main(int argc,char**argv){
    int NKEYS=50;
    if(argc>1)NKEYS=atoi(argv[1]);
    if(argc>2)rng_state^=strtoull(argv[2],0,0)*0x9e3779b97f4a7c15ULL;
    saes_init(); build_ddt(); selftest();
    honest_test(NKEYS);
    e2e_attack(NKEYS,NKEYS<=20);
    return 0;
}
