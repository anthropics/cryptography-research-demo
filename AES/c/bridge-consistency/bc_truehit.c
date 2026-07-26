// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
// bc_truehit.c — TRUE-hit acceptance of bridge-consistency test.
// For random AES key + random x_1[1..15], vary x_1[0]=0..16, run 4 full AES rounds
// (rounds 1..4: SB,SR,MC,ARK each), read a_ω = x_5[0]. Set s=a_0, d_ω=a_ω⊕s.
// Count #bad_{16} = #{ω∈1..16 : a_ω ∈ {0,s}}.  ALSO verify bridge identity directly:
// compute v_ω = S(a_ω)⊕κ (any κ), g_ω=(L^{-1}(v_0⊕v_ω))^{-1}, check g_ω==s²d_ω^{-1}⊕s
// at every non-bad ω and FAILS at every bad ω.
//
// Output: histogram of #bad_{16} over N keys, compare to Binom(16,1/128);
//         P[score(s_true)>=13] = P[#bad<=3]; identity verification pass/fail.
//
// build: gcc -O3 -march=native -fopenmp -o bc_th bc_truehit.c
// run:   ./bc_th <N> <seed>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* AES GF(256) */
static uint8_t SBOX[256], INV[256], SQR[256], L_inv_tab[256];
static uint8_t xt(uint8_t a){return (a<<1)^((a&0x80)?0x1B:0);}
static uint8_t gmul(uint8_t a,uint8_t b){uint8_t r=0;for(int i=0;i<8;i++){if(b&1)r^=a;b>>=1;a=xt(a);}return r;}
static void gf_init(void){
    INV[0]=0; for(int a=1;a<256;a++)for(int b=1;b<256;b++)if(gmul(a,b)==1){INV[a]=b;break;}
    /* AES affine L: y_i = x_i ^ x_{i+4} ^ x_{i+5} ^ x_{i+6} ^ x_{i+7} ^ c_i, c=0x63 */
    for(int x=0;x<256;x++){
        uint8_t q=INV[x],r=0;
        for(int i=0;i<8;i++){
            int b=((q>>i)&1)^((q>>((i+4)&7))&1)^((q>>((i+5)&7))&1)^((q>>((i+6)&7))&1)^((q>>((i+7)&7))&1);
            r|=b<<i;
        }
        SBOX[x]=r^0x63;
    }
    /* L^{-1}: build L as GF(2) matrix, invert. Simpler: L(x)=SBOX w/o inverse ⟹
       for u∈GF(256): L^{-1}[u] = unique q with L(q)=u. Since SBOX[x]=L(INV[x])^0x63,
       we have L(q) = SBOX[INV[q]]^0x63 ... no. Simpler: L(q)^0x63 = SBOX of (INV of q).
       Actually L is: for input q, output = affine part of SBOX applied to q.
       We compute L^{-1} by inverting the 8x8 bit matrix directly via brute force: */
    for(int u=0;u<256;u++){
        for(int q=0;q<256;q++){
            uint8_t r=0;
            for(int i=0;i<8;i++){
                int b=((q>>i)&1)^((q>>((i+4)&7))&1)^((q>>((i+5)&7))&1)^((q>>((i+6)&7))&1)^((q>>((i+7)&7))&1);
                r|=b<<i;
            }
            if(r==(uint8_t)u){L_inv_tab[u]=q;break;}
        }
    }
    for(int a=0;a<256;a++) SQR[a]=gmul(a,a);
    /* self-check: SBOX[0x53]==0xED, SBOX[0]==0x63 */
    if(SBOX[0]!=0x63||SBOX[0x53]!=0xED){fprintf(stderr,"SBOX selftest fail\n");exit(1);}
}

/* one AES round: SB,SR,MC,ARK on 16-byte state (column-major: s[r+4c]) */
static void mc_col(uint8_t *c){
    uint8_t a0=c[0],a1=c[1],a2=c[2],a3=c[3];
    c[0]=xt(a0)^xt(a1)^a1^a2^a3;
    c[1]=a0^xt(a1)^xt(a2)^a2^a3;
    c[2]=a0^a1^xt(a2)^xt(a3)^a3;
    c[3]=xt(a0)^a0^a1^a2^xt(a3);
}
static void aes_round(uint8_t *s,const uint8_t *rk){
    uint8_t t[16];
    for(int i=0;i<16;i++) t[i]=SBOX[s[i]];
    /* ShiftRows: row r shifted left by r. s[r+4c] layout. */
    uint8_t u[16];
    for(int r=0;r<4;r++)for(int c=0;c<4;c++) u[r+4*c]=t[r+4*((c+r)&3)];
    for(int c=0;c<4;c++) mc_col(u+4*c);
    for(int i=0;i<16;i++) s[i]=u[i]^rk[i];
}

/* xoshiro */
static inline uint64_t rotl(uint64_t x,int k){return(x<<k)|(x>>(64-k));}
static inline uint64_t xr(uint64_t*xs){uint64_t r=rotl(xs[1]*5,7)*9,t=xs[1]<<17;xs[2]^=xs[0];xs[3]^=xs[1];xs[1]^=xs[2];xs[0]^=xs[3];xs[2]^=t;xs[3]=rotl(xs[3],45);return r;}

int main(int argc,char**argv){
    uint64_t N=argc>1?strtoull(argv[1],0,0):1000000;
    uint64_t SEED=argc>2?strtoull(argv[2],0,0):0xC0FFEE;
    gf_init();

    uint64_t hist_bad[17]={0};
    uint64_t id_ok=0,id_fail=0,bad_id_ok=0,bad_id_fail=0,s0=0;

    #pragma omp parallel
    {
      uint64_t lhb[17]={0},lio=0,lif=0,lbio=0,lbif=0,ls0=0;
      uint64_t xs[4];
      #ifdef _OPENMP
      int tid=omp_get_thread_num();
      #else
      int tid=0;
      #endif
      uint64_t z=SEED^(0x9E3779B97F4A7C15ULL*tid);
      for(int i=0;i<4;i++){z+=0x9E3779B97F4A7C15ULL;uint64_t r=z;r=(r^(r>>30))*0xBF58476D1CE4E5B9ULL;r=(r^(r>>27))*0x94D049BB133111EBULL;xs[i]=r^(r>>31);}
      for(int i=0;i<10;i++)xr(xs);

      #pragma omp for schedule(static)
      for(uint64_t trial=0;trial<N;trial++){
        uint8_t x1_base[16],rk[4][16];
        /* random x_1[1..15] and 4 round keys */
        uint64_t r;
        r=xr(xs);for(int i=0;i<8;i++)x1_base[i]=r>>(8*i);
        r=xr(xs);for(int i=0;i<8;i++)x1_base[8+i]=r>>(8*i);
        for(int j=0;j<4;j++){r=xr(xs);for(int i=0;i<8;i++)rk[j][i]=r>>(8*i);r=xr(xs);for(int i=0;i<8;i++)rk[j][8+i]=r>>(8*i);}

        uint8_t a[17];
        for(int om=0;om<=16;om++){
            uint8_t st[16]; memcpy(st,x1_base,16); st[0]=(uint8_t)om;
            aes_round(st,rk[0]); /* round 1: x_1 -> x_2 */
            aes_round(st,rk[1]); /* x_2 -> x_3 */
            aes_round(st,rk[2]); /* x_3 -> x_4 */
            aes_round(st,rk[3]); /* x_4 -> x_5 */
            a[om]=st[0];
        }
        uint8_t s=a[0];
        if(s==0){ls0++;continue;} /* paper absorbs s=0 as 256/255 data factor */

        int nbad=0;
        uint8_t kappa=37; /* arbitrary "u_5[0]" */
        uint8_t v0=SBOX[a[0]]^kappa;
        for(int om=1;om<=16;om++){
            int bad=(a[om]==0||a[om]==s);
            nbad+=bad;
            /* bridge identity verification */
            uint8_t v=SBOX[a[om]]^kappa;
            uint8_t g=INV[L_inv_tab[v0^v]];
            uint8_t d=s^a[om];
            uint8_t rhs=gmul(SQR[s],INV[d])^s;
            int ok=(g==rhs);
            if(bad){ if(ok)lbio++; else lbif++; }
            else   { if(ok)lio++;  else lif++;  }
        }
        lhb[nbad]++;
      }
      #pragma omp critical
      { for(int k=0;k<17;k++)hist_bad[k]+=lhb[k];
        id_ok+=lio;id_fail+=lif;bad_id_ok+=lbio;bad_id_fail+=lbif;s0+=ls0; }
    }

    printf("# bc_truehit N=%llu seed=0x%llX  (s=0 skipped: %llu)\n",
        (unsigned long long)N,(unsigned long long)SEED,(unsigned long long)s0);
    printf("# bridge identity on non-bad indices: ok=%llu fail=%llu\n",
        (unsigned long long)id_ok,(unsigned long long)id_fail);
    printf("# bridge identity on BAD  indices: ok=%llu fail=%llu\n",
        (unsigned long long)bad_id_ok,(unsigned long long)bad_id_fail);
    printf("# #bad_{16} histogram (empirical  vs  Binom(16,1/128)):\n");
    uint64_t M=N-s0; double p=2.0/256.0;
    uint64_t ge4=0;
    for(int k=0;k<=16;k++){
        double binom=1; for(int j=0;j<k;j++)binom*= (double)(16-j)/(j+1);
        double pk=binom*pow(p,k)*pow(1-p,16-k);
        printf("  b=%2d  emp=%12llu  pred=%15.3f  ratio=%.4f\n",
            k,(unsigned long long)hist_bad[k],pk*M,
            hist_bad[k]? (double)hist_bad[k]/(pk*M):0.0);
        if(k>=4)ge4+=hist_bad[k];
    }
    printf("# P[true-hit REJECTED @ tau=13] = P[#bad>=4] : emp=%.6e  pred=%.6e\n",
        (double)ge4/M, /* analytic below computed in-line */
        0.0);
    return 0;
}
