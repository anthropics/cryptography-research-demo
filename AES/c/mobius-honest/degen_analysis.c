// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* degen_analysis.c — characterize when all P_m=0 on small-AES δ-sets.
 * Hypothesis: degenerate ⟺ odd-support(H) ∈ {∅, GF(q)*, GF(q)}. */
#include "saes.h"
#include <stdio.h>
#include <stdlib.h>
static uint64_t rng_s=0xdeadbeef12345678ULL;
static uint32_t rng(void){rng_s=rng_s*6364136223846793005ULL+1;return rng_s>>32;}
static uint8_t rn(void){return rng()&15;}

#define NPM 6
static const int Pe[NPM]={3,5,7,6,9,11};
static void cPm(const uint8_t*S,int n,uint8_t Pm[NPM]){
    for(int k=0;k<NPM;k++)Pm[k]=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){uint8_t d=S[i]^S[j];if(!d)continue;
        for(int k=0;k<NPM;k++)Pm[k]^=sgf_pow(d,Pe[k]);}
}
int main(void){
    saes_init();
    int NK=5000;
    int degen=0, oddsup_full=0, oddsup_empty=0, oddsup_other=0;
    for(int t=0;t<NK;t++){
        uint8_t key[16],rk[11][16];for(int i=0;i<16;i++)key[i]=rn();
        skey_expand(key,rk);
        uint8_t x1b[16];for(int i=0;i<16;i++)x1b[i]=rn();
        uint8_t x5[16];
        for(int v=0;v<16;v++){uint8_t x1[16],s[16];memcpy(x1,x1b,16);x1[0]=v;
            memcpy(s,x1,16);sark(s,rk[1]);simc(s);sisr(s);sisb(s);sark(s,rk[0]);
            uint8_t tr[8][16];senc7_trace(rk,s,tr);x5[v]=tr[5][0];}
        int vref=x1b[0];uint8_t a=x5[vref];if(a==0)continue;
        /* E'={1/Δx5} */
        uint8_t E[16];int nE=0;
        for(int v=0;v<16;v++){if(v==vref)continue;uint8_t dx=x5[v]^a;if(dx)E[nE++]=sgf_inv(dx);}
        uint8_t Pm[NPM];cPm(E,nE,Pm);
        int dg=1;for(int k=0;k<NPM;k++)if(Pm[k])dg=0;
        if(!dg)continue;
        degen++;
        /* odd-support of E */
        int h[16]={0};for(int i=0;i<nE;i++)h[E[i]]^=1;
        int os=0;for(int i=0;i<16;i++)os+=h[i];
        if(os==0)oddsup_empty++;
        else if(os==15&&h[0]==0)oddsup_full++;
        else oddsup_other++;
        if(degen<=10){
            printf("  degen#%d: |os|=%d os={",degen,os);
            for(int i=0;i<16;i++)if(h[i])printf("%x ",i);printf("} ");
            /* check if affine subspace (coset of GF(2)-linear) */
            int is_aff=0;
            if(os>=1){
                uint8_t base=0;for(int i=0;i<16;i++)if(h[i]){base=i;break;}
                /* shift to contain 0 */
                int h2[16]={0};for(int i=0;i<16;i++)if(h[i])h2[i^base]=1;
                /* check closure under XOR */
                is_aff=1;
                for(int a=0;a<16;a++)if(h2[a])for(int b=0;b<16;b++)if(h2[b])
                    if(!h2[a^b]){is_aff=0;break;}
            }
            printf("affine=%d nE=%d\n",is_aff,nE);
        }
    }
    /* Also: for RANDOM 15-element multisets, how often degenerate? */
    int rdg=0;
    for(int t=0;t<100000;t++){
        uint8_t R[15];for(int i=0;i<15;i++)R[i]=1+rn()%15;
        uint8_t Pm[NPM];cPm(R,15,Pm);int dg=1;for(int k=0;k<NPM;k++)if(Pm[k])dg=0;
        if(dg)rdg++;
    }
    printf("random 15-mset degen: %d/100000 (%.4f%%)\n",rdg,100.0*rdg/100000);
    printf("degen/NK = %d/%d (%.3f%%)\n",degen,NK,100.0*degen/NK);
    printf("  odd-support=∅: %d, =GF(16)*: %d, other: %d\n",
           oddsup_empty,oddsup_full,oddsup_other);
    return 0;
}
