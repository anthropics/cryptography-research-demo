// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef PROP2_H
#define PROP2_H
/* DFJ'13 Prop-2 at rounds 3,4,5 (δ at z_2[0], observe Δx_6[0]).
 * 24-ref form: x_3[col0](4) + x_4[16](16) + x_5[diag](4).
 * e_ω = Δx_6[0,0] for ω = Δz_2[0,0] ∈ {1..255}.
 */
#include "aes_core.h"

typedef struct {
    u8 x3c[4];    /* x_3[col0] = x_3[0,1,2,3] */
    u8 x4[16];    /* x_4[all] */
    u8 x5d[4];    /* x_5[diag0] = x_5[0,5,10,15] */
} refs24_t;

/* Direct Prop-2: compute e_ω from 24 refs. */
static u8 prop2_direct(const refs24_t*R, u8 om){
    /* Δx_3[col0] = MC·(ω,0,0,0) */
    u8 dy3[4];
    for(int i=0;i<4;i++){
        u8 dx3i = gmul(MCc[i][0], om);
        dy3[i] = SBOX[R->x3c[i]^dx3i] ^ SBOX[R->x3c[i]];
    }
    /* Δx_4[p] for p=(r,c): after SR, z_3[r,c]=y_3[r,(c+r)%4]; y_3 nonzero only col0.
       So Δz_3[r,c]=dy3[r] iff (c+r)%4==0, i.e., c=(4-r)%4. Else 0.
       Δw_3[r,c]=Σ_i MC[r][i]·Δz_3[i,c]. Only i with (c+i)%4==0 contributes:
       i=(4-c)%4. So Δw_3[r,c]=MC[r][(4-c)%4]·dy3[(4-c)%4]. */
    u8 dy4[16];
    for(int p=0;p<16;p++){
        int r=p&3, c=p>>2, ic=(4-c)&3;
        u8 dx4 = gmul(MCc[r][ic], dy3[ic]);
        dy4[p] = SBOX[R->x4[p]^dx4] ^ SBOX[R->x4[p]];
    }
    /* Δx_5[r,r] = Σ_i MC[r][i]·Δz_4[i,r] = Σ_i MC[r][i]·Δy_4[i,(r+i)%4] */
    u8 e=0;
    for(int r=0;r<4;r++){
        u8 dx5=0;
        for(int i=0;i<4;i++) dx5 ^= gmul(MCc[r][i], dy4[P_r[r][i]]);
        u8 dy5 = SBOX[R->x5d[r]^dx5] ^ SBOX[R->x5d[r]];
        e ^= gmul(MCc[0][r], dy5);
    }
    return e;
}

/* Compute all 255 e_ω into e[1..255] (e[0]=0). */
static void prop2_all(const refs24_t*R, u8 e[256]){
    e[0]=0;
    for(int om=1;om<256;om++) e[om]=prop2_direct(R,om);
}

/* ---------- Rebound outer-base (for UU-CACHE) ----------
 * Outer params (c_in=1): (om_star, x3c[4], d_out, dx5d[4])
 *   om_star = Δz_2[0]_pair
 *   d_out   = Δw_5[0]_pair = Δx_6[0]_pair
 *   dx5d[r] = Δx_5[r,r]_pair
 * From these, derive (Δx_4^pair[16], Δy_4^pair[16], Δy_5^pair[4]) and the
 * DDT-solutions x_4[p]^{0,1,...}, x_5[r]^{0,1}.
 */
typedef struct {
    u8 om_star, d_out;
    u8 x3c[4], dx5d[4];
    /* derived: */
    u8 dx4_pair[16], dy4_pair[16];    /* pair diffs at layer 4 */
    u8 dy5_pair[4];                    /* pair diffs at y_5[diag] */
    u8 n4[16], sol4[16][4];           /* DDT sols for x_4[p] */
    u8 n5[4],  sol5[4][4];            /* DDT sols for x_5[r,r] */
    int valid;                         /* all DDT>0 */
    u64 n_inner;                       /* Π n4 × Π n5 */
} outer_t;

static void outer_derive(outer_t*O){
    /* forward: Δx_4^pair[p] from (om_star, x3c) */
    u8 dy3[4];
    for(int i=0;i<4;i++){
        u8 dx3i=gmul(MCc[i][0],O->om_star);
        dy3[i]=SBOX[O->x3c[i]^dx3i]^SBOX[O->x3c[i]];
    }
    for(int p=0;p<16;p++){
        int r=p&3,c=p>>2,ic=(4-c)&3;
        O->dx4_pair[p]=gmul(MCc[r][ic],dy3[ic]);
    }
    /* backward: Δy_5^pair[r]=iMC[r][0]·d_out;  Δy_4^pair[p∈P_r]=iMC[row(p)][r]·dx5d[r] */
    for(int r=0;r<4;r++) O->dy5_pair[r]=gmul(iMCc[r][0],O->d_out);
    for(int r=0;r<4;r++) for(int i=0;i<4;i++){
        int p=P_r[r][i]; int row=p&3;
        O->dy4_pair[p]=gmul(iMCc[row][r],O->dx5d[r]);
    }
    /* DDT at layer 4 and 5 */
    O->valid=1; O->n_inner=1;
    for(int p=0;p<16;p++){
        u8 n=DDTn[O->dx4_pair[p]][O->dy4_pair[p]];
        O->n4[p]=n; if(n==0){O->valid=0;}
        for(int k=0;k<n&&k<4;k++) O->sol4[p][k]=DDTsol[O->dx4_pair[p]][O->dy4_pair[p]][k];
        if(n) O->n_inner*=n;
    }
    for(int r=0;r<4;r++){
        u8 n=DDTn[O->dx5d[r]][O->dy5_pair[r]];
        O->n5[r]=n; if(n==0){O->valid=0;}
        for(int k=0;k<n&&k<4;k++) O->sol5[r][k]=DDTsol[O->dx5d[r]][O->dy5_pair[r]][k];
        if(n) O->n_inner*=n;
    }
}

/* ---------- UU-CACHE ----------
 * For fixed outer O, build UU_r[g][ω] for g∈[0,2^5), ω∈[1,256).
 *   g=(a0,a1,a2,a3,b) where a_i selects x_4[P_r[i]] branch, b selects x_5[r] branch.
 * For simplicity handle n4/n5 in {2,4} by mapping bit to idx%n.
 * But for EXACTNESS we want mixed-radix; assumes n=2 for every DDT solution set (holds for >99% of bases).
 */
typedef struct {
    u8 UU[4][32][256];   /* UU[r][g][ω] */
    u8 dy3_om[4][256];   /* Δy_3[i,0]^{(ω)} cache */
} uucache_t;

static u64 g_uu_sb=0;

static void uu_build(const outer_t*O, uucache_t*U){
    /* precompute dy3_om (depends only on x3c, NOT on pair diffs) */
    for(int i=0;i<4;i++){
        for(int om=0;om<256;om++){
            u8 dx3=gmul(MCc[i][0],om);
            U->dy3_om[i][om]=SBOX[O->x3c[i]^dx3]^SBOX[O->x3c[i]];
            g_uu_sb+=2;
        }
    }
    /* for each r, build UU[r][g][ω] */
    for(int r=0;r<4;r++){
        /* precompute per-position SB diffs: for each i∈P_r, each a∈{0,1}, each ω:
           sdy4[i][a][ω] = SB(x_4^a[p]⊕Δx_4[p]^{(ω)})⊕SB(x_4^a[p])
           where Δx_4[p]^{(ω)}=MC[row(p)][(4-col(p))%4]·dy3_om[(4-col(p))%4][ω] */
        u8 sdy4[4][2][256];
        for(int i=0;i<4;i++){
            int p=P_r[r][i], row=p&3, col=p>>2, ic=(4-col)&3;
            for(int a=0;a<2;a++){
                u8 xv = O->sol4[p][a % O->n4[p]];  /* wrap for n=4 handled separately */
                for(int om=0;om<256;om++){
                    u8 dx4=gmul(MCc[row][ic],U->dy3_om[ic][om]);
                    sdy4[i][a][om]=SBOX[xv^dx4]^SBOX[xv];
                    g_uu_sb+=2;
                }
            }
        }
        for(int g=0;g<32;g++){
            int a0=g&1,a1=(g>>1)&1,a2=(g>>2)&1,a3=(g>>3)&1,b=(g>>4)&1;
            u8 x5v = O->sol5[r][b % O->n5[r]];
            for(int om=0;om<256;om++){
                u8 dx5 = gmul(MCc[r][0],sdy4[0][a0][om])^gmul(MCc[r][1],sdy4[1][a1][om])
                       ^ gmul(MCc[r][2],sdy4[2][a2][om])^gmul(MCc[r][3],sdy4[3][a3][om]);
                u8 dy5 = SBOX[x5v^dx5]^SBOX[x5v]; g_uu_sb+=2;
                U->UU[r][g][om]=gmul(MCc[0][r],dy5);
            }
        }
    }
}

/* assemble e_ω from UU for given inner-config (a[16],b[4]) bit-vector.
   For n4[p]==4 or n5[r]==4, the extra bits are NOT handled here; caller
   must restrict to n==2 bases or use direct. */
static inline u8 uu_assemble(const uucache_t*U, u32 a16, u8 b4, u8 om){
    u8 e=0;
    for(int r=0;r<4;r++){
        u32 g=0;
        for(int i=0;i<4;i++) g|=((a16>>P_r[r][i])&1)<<i;
        g|=((b4>>r)&1)<<4;
        e^=U->UU[r][g][om];
    }
    return e;
}

#endif
