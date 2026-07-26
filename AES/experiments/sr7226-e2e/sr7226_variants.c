// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* sr7226.c — small-scale AES SR(7,2,2,6) + I_{m,n} Mobius-invariant checks.
 *
 * Cipher: 2x2 state of 6-bit cells over GF(2^6) = GF(2)[x]/(x^6+x+1), 7 rounds.
 *   State is column-major: idx = row + 2*col  ->  cells {0,1}=col0, {2,3}=col1.
 *   SB : S(t) = L(t^-1) XOR 0x2A,  L a fixed GF(2)-linear bijection on 6 bits
 *        (circulant row (1,1,0,1,0,0): L(b)_i = b_i ^ b_{i+1} ^ b_{i+3}).
 *   SR : row 0 fixed, row 1 rotl 1  ==>  swap cells 1 <-> 3.
 *   MC : circ(2,3) over GF(64); det = 2^2 ^ 3^2 = 1  ==>  MC^{-1} = MC.
 *   Last round omits MC.  Whitening key is rk[0] (= paper's k_{-1}).
 *
 * What this file verifies (report.tex §3.1, on the 6-bit toy):
 *   (a) bridge identity  g_w = s^2 * d_w^{-1} XOR s     (on w with a_w not in {0,s})
 *   (b) core I_{m,n} invariant on good w
 *   (c) full delta-set + add-0 parity patch: exactly one parity matches, = (#bad mod 2)
 *   (d) first-nonzero-P_m fingerprint: online/offline zero-patterns agree, values match
 *   (e) wrong-k_6 sweep: only the true k_6[0,3] matches
 *   (f) fingerprint distribution on random inputs ~ uniform on (GF(64)^x)^4
 *
 * Build:  cc -O2 -o sr7226 sr7226.c
 * Run:    ./sr7226                    # 200 correctness trials
 *         ./sr7226 wrong [N]          # N wrong-k6 sweeps (default 5)
 *         ./sr7226 dist  [N]          # N-sample distribution test (default 200000)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#endif

/* ------------------------------------------------------------------ */
/* GF(2^6) with modulus x^6 + x + 1  (0x43)                            */
/* ------------------------------------------------------------------ */
#define Q   64
#define QM1 63
typedef uint8_t gf;

static gf  EXP[2*QM1];
static int LOG[Q];
static gf  INV[Q];

static gf gf_mul_slow(gf a, gf b) {
    unsigned r = 0;
    for (int i = 0; i < 6; i++) {
        if (b & 1) r ^= a;
        b >>= 1;
        a <<= 1;
        if (a & 0x40) a ^= 0x43;
    }
    return (gf)r;
}
static void gf_init(void) {
    /* find a primitive element */
    gf g = 0;
    for (gf cand = 2; cand < Q; cand++) {
        gf x = 1; int ok = 1;
        for (int i = 1; i < QM1; i++) { x = gf_mul_slow(x, cand); if (x == 1) { ok = 0; break; } }
        if (ok) { g = cand; break; }
    }
    assert(g);
    gf x = 1;
    for (int i = 0; i < QM1; i++) { EXP[i] = x; LOG[x] = i; x = gf_mul_slow(x, g); }
    for (int i = 0; i < QM1; i++) EXP[i+QM1] = EXP[i];
    LOG[0] = -1;
    INV[0] = 0;
    for (int a = 1; a < Q; a++) INV[a] = EXP[QM1 - LOG[a]];
}
static inline gf gf_mul(gf a, gf b) { return (a && b) ? EXP[LOG[a] + LOG[b]] : 0; }
static inline gf gf_inv(gf a)       { return INV[a]; }
static inline gf gf_pow(gf a, int e){
    if (!a) return e ? 0 : 1;
    return EXP[((long)LOG[a] * (e % QM1) + (long)QM1*QM1) % QM1];
}
static gf rnd6(void){ return (gf)(rand()&(Q-1)); }

/* ------------------------------------------------------------------ */
/* S-box:  S(t) = L(t^{-1}) ^ AFFC                                     */
/* ------------------------------------------------------------------ */
#define AFFC 0x2A
static gf LTBL[Q], LINV[Q], SBOX[Q], SINV[Q];

static gf rol6(gf b, int r) { return (gf)(((b << r) | (b >> (6 - r))) & (Q-1)); }
static void sbox_init(void) {
    for (int b = 0; b < Q; b++)
        LTBL[b] = (gf)((b ^ rol6(b,1) ^ rol6(b,3)) & (Q-1));
    /* verify L bijective, build L^{-1} */
    int seen[Q] = {0};
    for (int b = 0; b < Q; b++) { assert(!seen[LTBL[b]]); seen[LTBL[b]] = 1; LINV[LTBL[b]] = (gf)b; }
    for (int b = 0; b < Q; b++) {
        SBOX[b] = (gf)(LTBL[gf_inv((gf)b)] ^ AFFC);
    }
    for (int b = 0; b < Q; b++) SINV[SBOX[b]] = (gf)b;
    assert(SBOX[0] == AFFC);
}

/* ------------------------------------------------------------------ */
/* Round functions (2x2 state, column-major indices 0..3)              */
/* ------------------------------------------------------------------ */
#define NCELLS  4
#define NROUNDS 7
#define MCa 2
#define MCb 3

static inline void sub_cells (gf *s) { for (int i=0;i<4;i++) s[i]=SBOX[s[i]]; }
static inline void isub_cells(gf *s) { for (int i=0;i<4;i++) s[i]=SINV[s[i]]; }
static inline void shift_rows(gf *s) { gf t=s[1]; s[1]=s[3]; s[3]=t; }   /* SR = SR^{-1} */
static inline void mix_cols  (gf *s) {                                   /* MC = MC^{-1} */
    for (int c=0;c<2;c++){
        gf a=s[2*c], b=s[2*c+1];
        s[2*c  ]=gf_mul(MCa,a)^gf_mul(MCb,b);
        s[2*c+1]=gf_mul(MCb,a)^gf_mul(MCa,b);
    }
}
static inline void add_key(gf *s, const gf *k){ for(int i=0;i<4;i++) s[i]^=k[i]; }

/* key schedule: rk[0..NROUNDS], each 4 cells */
static gf RCON[NROUNDS];
static void key_schedule(const gf *mk, gf rk[NROUNDS+1][4]) {
    memcpy(rk[0], mk, 4);
    for (int i=0;i<NROUNDS;i++){
        gf t0 = (gf)(SBOX[rk[i][3]] ^ RCON[i]);
        gf t1 = SBOX[rk[i][2]];
        rk[i+1][0]=rk[i][0]^t0;
        rk[i+1][1]=rk[i][1]^t1;
        rk[i+1][2]=rk[i][2]^rk[i+1][0];
        rk[i+1][3]=rk[i][3]^rk[i+1][1];
    }
}

/* encrypt, capturing x_i (round inputs) into xs[0..NROUNDS-1] */
static void encrypt_trace(const gf *pt, gf rk[NROUNDS+1][4],
                          gf *ct, gf xs[NROUNDS][4]) {
    gf s[4]; memcpy(s,pt,4);
    add_key(s, rk[0]);               memcpy(xs[0],s,4);          /* x_0 */
    for (int r=0;r<NROUNDS-1;r++){
        sub_cells(s); shift_rows(s); mix_cols(s); add_key(s,rk[r+1]);
        memcpy(xs[r+1],s,4);                                     /* x_{r+1} */
    }
    sub_cells(s); shift_rows(s); add_key(s,rk[NROUNDS]);         /* last: no MC */
    memcpy(ct,s,4);
}

/* ------------------------------------------------------------------ */
/* delta-set: 64 plaintexts whose x_1 varies only in cell 0.           */
/* Built with knowledge of the true key (invert round 0).              */
/* ------------------------------------------------------------------ */
static void build_delta_set(gf rk[NROUNDS+1][4], const gf base_x1[4],
                            gf pts[Q][4]) {
    for (int i=0;i<Q;i++){
        gf s[4]={ (gf)i, base_x1[1], base_x1[2], base_x1[3] };
        add_key(s,rk[1]); mix_cols(s); shift_rows(s); isub_cells(s); /* -> x_0 */
        add_key(s,rk[0]);                                            /* -> P   */
        memcpy(pts[i],s,4);
    }
}

/* ------------------------------------------------------------------ */
/* Online v = MC^{-1}(x_6)[0] from ciphertext + k_6 guess at cells 0,3 */
/* y_6[col0] = SR^{-1}(C xor k_6)[col0] uses (C^k6)[0] and (C^k6)[3]   */
/* ------------------------------------------------------------------ */
static inline gf online_v(const gf *ct, const gf *k6){
    gf y0 = (gf)(ct[0]^k6[0]);
    gf y1 = (gf)(ct[3]^k6[3]);
    gf x60 = SINV[y0], x61 = SINV[y1];               /* x_6[0], x_6[1] */
    return (gf)(gf_mul(MCa,x60) ^ gf_mul(MCb,x61));  /* MC^{-1}=MC */
}

/* ------------------------------------------------------------------ */
/* I_{m,n} machinery.  E = (5,11,13,23,31): Frobenius-coset reps,      */
/* coprime to 63, wt>=2.  Max exponent 31.                             */
/* ------------------------------------------------------------------ */
static const int E_LIST[] = {5,11,13,23,31};
#define NE ((int)(sizeof E_LIST/sizeof *E_LIST))
#define RMAX 32

#ifdef __AVX2__
static gf POW[Q][RMAX] __attribute__((aligned(32)));   /* POW[g][r] = g^r, r=0..31 */
#else
static gf POW[Q][RMAX];
#endif
static void pow_init(void){
    for(int g=0; g<Q; g++) for(int r=0;r<RMAX;r++) POW[g][r]=gf_pow((gf)g,r);
}

/* p_r = XOR_w vals[w]^r  (w from 1 to n-1) */
static void single_psums(const gf *vals, int n, gf p[RMAX]){
#ifdef __AVX2__
    __m256i acc=_mm256_setzero_si256();
    for(int w=1; w<n; w++)
        acc=_mm256_xor_si256(acc,_mm256_load_si256((const __m256i*)POW[vals[w]]));
    _mm256_storeu_si256((__m256i*)p,acc);
#else
    memset(p,0,RMAX);
    for(int w=1; w<n; w++){
        const gf *row = POW[vals[w]];
        for(int r=0;r<RMAX;r++) p[r]^=row[r];
    }
#endif
}
/* Lucas reduction, m odd */
static gf Pm_from_p(const gf p[RMAX], int m){
    gf pm = (p[0]&1) ? p[m] : 0;
    for(int k=(m-1)&m; k; k=(k-1)&m){
        int j=m-k;
        if(k<j) pm ^= gf_mul(p[k],p[j]);
    }
    return pm;
}
/* naive O(n^2) for cross-check */
static gf Pm_naive(const gf *vals, int n, int m){
    gf acc=0;
    for(int i=1;i<n;i++) for(int j=i+1;j<n;j++) acc ^= gf_pow((gf)(vals[i]^vals[j]),m);
    return acc;
}

/* first-nonzero fingerprint: <=4 ratios in GF(64)^x, or ncells=0 on fail */
typedef struct { int ncells; gf cell[NE-1]; } fp_t;

static fp_t fingerprint_nz(const gf *vals, int n){
    gf p[RMAX]; single_psums(vals,n,p);
    gf P[NE]; int nz[NE], k=0;
    for(int i=0;i<NE;i++){ P[i]=Pm_from_p(p,E_LIST[i]); if(P[i]) nz[k++]=i; }
    fp_t f={0};
    if(k<2) return f;
    int m0=E_LIST[nz[0]]; gf Pm0=P[nz[0]];
    for(int j=1;j<k && f.ncells<NE-1;j++){
        int ni=E_LIST[nz[j]]; gf Pn=P[nz[j]];
        f.cell[f.ncells++]=gf_mul(gf_pow(Pm0,ni), gf_inv(gf_pow(Pn,m0)));
    }
    return f;
}
static int fp_eq(fp_t a, fp_t b){
    if(a.ncells!=b.ncells) return 0;
    for(int i=0;i<a.ncells;i++) if(a.cell[i]!=b.cell[i]) return 0;
    return 1;
}

/* ================================================================== */
/* χ-canonicalization fingerprint (report.tex appendix)                */
/*   χ(V) ∈ GF(2)^64, bit d = parity of |{ω : v_ω = d}|.               */
/*   Canonical form = AGL(1,64)-orbit representative.                  */
/* ================================================================== */
static gf  MULTBL[Q][Q];     /* a*b */
static gf  TRTBL[Q];         /* absolute trace GF(64)->GF(2) */
static gf  HTTBL[Q];         /* half-trace: HT(c)^2 ^ HT(c) = c when Tr(c)=0 */

static gf TAU1=0;  /* any fixed element with Tr=1 */
static void chi_init(void){
    for(int a=0;a<Q;a++) for(int b=0;b<Q;b++) MULTBL[a][b]=gf_mul((gf)a,(gf)b);
    for(int c=0;c<Q;c++){
        gf t=0,x=(gf)c; for(int i=0;i<6;i++){ t^=x; x=gf_mul(x,x); } TRTBL[c]=(gf)(t&1);
    }
    /* solve u^2+u=c by table (half-trace formula is for odd-degree fields) */
    for(int c=0;c<Q;c++) HTTBL[c]=0xFF;
    for(int u=0;u<Q;u++){ gf c=(gf)(gf_mul((gf)u,(gf)u)^u); if(HTTBL[c]==0xFF) HTTBL[c]=(gf)u; }
    for(int c=0;c<Q;c++) assert((TRTBL[c]==0)==(HTTBL[c]!=0xFF));
    for(int c=1;c<Q;c++) if(TRTBL[c]){ TAU1=(gf)c; break; }
    assert(TAU1);
}

static inline uint64_t permute_chi(uint64_t chi, gf alpha, gf beta){
    uint64_t out=0;
    while(chi){
        int d=__builtin_ctzll(chi);
        chi &= chi-1;
        out |= 1ULL << (MULTBL[alpha][d]^beta);
    }
    return out;
}

/* Brute-force canonical: min over all (α,β) ∈ AGL(1,64). */
static uint64_t chi_canon_brute(uint64_t chi){
    uint64_t best=~0ULL;
    for(int a=1;a<Q;a++) for(int b=0;b<Q;b++){
        uint64_t v=permute_chi(chi,(gf)a,(gf)b);
        if(v<best) best=v;
    }
    return best;
}

static long CHI_BRUTE_FALLBACKS=0;
/* Moving-frame canonical.  Returns canonical χ (or falls back to brute). */
/* n_odd = |V| mod 2.  p[] = single-index power sums (already computed).  */
static uint64_t chi_canon_fast(uint64_t chi, const gf p[RMAX], int n_odd){
    /* pick first P_m != 0 with gcd(m,63)=1 */
    static const int M[5]   ={ 5,11,13,23,31};
    static const int MINV[5]={25,40,29,52, 2};   /* (63 - inv(m,63)) mod 63 */
    gf alpha=0;
    for(int i=0;i<5;i++){
        gf Pm=Pm_from_p(p,M[i]);
        if(Pm){ alpha=gf_pow(Pm,MINV[i]); break; }
    }
    if(!alpha){ __atomic_add_fetch(&CHI_BRUTE_FALLBACKS,1,__ATOMIC_RELAXED); return chi_canon_brute(chi); }  /* rare */

    gf S1=p[1];
    if(n_odd){
        gf beta=MULTBL[alpha][S1];
        return permute_chi(chi,alpha,beta);
    }
    /* n even: Artin–Schreier solve  u^2 ⊕ u = c  with c = S3/S1^3.
     * If Tr(c)=1 the equation is unsolvable; XOR in a fixed τ with Tr(τ)=1
     * (both sides take the same branch since c is AGL-invariant).
     * β* = α*·S1·u; two roots {u, u⊕1} → two candidates → min.        */
    if(!S1){
        /* ~1/Q: S1,S3,S5 all translation-invariant when S1=0; take min over β. */
        uint64_t best=~0ULL;
        for(int b=0;b<Q;b++){ uint64_t v=permute_chi(chi,alpha,(gf)b); if(v<best)best=v; }
        return best;
    }
    gf S3=p[3];
    gf c = gf_mul(S3, gf_inv(gf_mul(S1,gf_mul(S1,S1))));
    if(TRTBL[c]) c^=TAU1;
    gf u = HTTBL[c];
    gf aS1 = MULTBL[alpha][S1];
    gf b0 = MULTBL[aS1][u], b1 = (gf)(b0^aS1);
    uint64_t r0=permute_chi(chi,alpha,b0), r1=permute_chi(chi,alpha,b1);
    return r0<r1 ? r0 : r1;
}

/* Fingerprint of a value-list (index 0 ignored as usual). */
static uint64_t chi_fp(const gf *vals, int n){
    uint64_t chi=0; gf p[RMAX];
    for(int w=1;w<n;w++) chi ^= 1ULL<<vals[w];
    single_psums(vals,n,p);
    return chi_canon_fast(chi,p,(n-1)&1);
}
static uint64_t chi_fp_brute(const gf *vals, int n){
    uint64_t chi=0;
    for(int w=1;w<n;w++) chi ^= 1ULL<<vals[w];
    { __atomic_add_fetch(&CHI_BRUTE_FALLBACKS,1,__ATOMIC_RELAXED); return chi_canon_brute(chi); }
}

/* ------------------------------------------------------------------ */
/* χ verification mode                                                 */
/* ------------------------------------------------------------------ */
static void chi_test(int N){
    /* (i) AGL-invariance of brute (reference) and of fast (used in attack),   */
    /*     plus completeness: different orbits -> different fp.                */
    int brute_ok=0, fast_ok=0, fast_even_bad=0;
    for(int t=0;t<N;t++){
        gf V[Q+1]; V[0]=0; for(int w=1;w<=Q;w++) V[w]=rnd6();
        gf a=(gf)(1+rand()%QM1), b=rnd6();
        gf W[Q+1]; W[0]=0; for(int w=1;w<=Q;w++) W[w]=(gf)(MULTBL[a][V[w]]^b);
        /* brute (lex-min) is AGL-invariant by construction */
        assert(chi_fp_brute(V,Q)==chi_fp_brute(W,Q));
        assert(chi_fp_brute(V,Q+1)==chi_fp_brute(W,Q+1));
        brute_ok++;
        /* fast (moving-frame) must also be AGL-invariant */
        uint64_t fV0=chi_fp(V,Q),   fW0=chi_fp(W,Q);
        uint64_t fV1=chi_fp(V,Q+1), fW1=chi_fp(W,Q+1);
        if(fV0==fW0) fast_ok++; else assert(!"fast odd-n not AGL-invariant");
        if(fV1!=fW1) fast_even_bad++;
    }
    printf("chi brute AGL-inv: %d/%d OK;  fast odd-n AGL-inv: %d/%d OK;  "
           "fast even-n mismatches: %d\n", brute_ok,N,fast_ok,N,fast_even_bad);

    /* (ii) bridge match: online g vs offline dinv under real cipher */
    int bridge_ok=0, skip=0;
    for(int t=0;t<N;t++){
        gf mk[4]={rnd6(),rnd6(),rnd6(),rnd6()}, rk[NROUNDS+1][4];
        key_schedule(mk,rk);
        gf base_x1[4]={rnd6(),rnd6(),rnd6(),rnd6()}, pts[Q][4];
        build_delta_set(rk,base_x1,pts);
        gf ct[Q][4],x5_0[Q],xs[NROUNDS][4];
        for(int i=0;i<Q;i++){encrypt_trace(pts[i],rk,ct[i],xs);x5_0[i]=xs[5][0];}
        gf s=x5_0[0]; if(!s){skip++;continue;}
        gf dinv[Q+1],g[Q+1],v[Q];
        for(int w=0;w<Q;w++) dinv[w]=gf_inv(s^x5_0[w]);
        for(int w=0;w<Q;w++) v[w]=online_v(ct[w],rk[NROUNDS]);
        for(int w=0;w<Q;w++) g[w]=gf_inv(LINV[v[0]^v[w]]);
        dinv[Q]=0; g[Q]=0;
        uint64_t on0=chi_fp(g,Q),   on1=chi_fp(g,Q+1);
        uint64_t of0=chi_fp(dinv,Q),of1=chi_fp(dinv,Q+1);
        int bad=0; for(int w=1;w<Q;w++) if(x5_0[w]==0||x5_0[w]==s) bad++;
        assert((bad&1)? (on1==of1) : (on0==of0));
        bridge_ok++;
    }
    printf("chi bridge match (online g == offline dinv, parity-patched): %d/%d OK (%d s=0 skips)\n",
           bridge_ok,N,skip);
}

/* ================================================================== */
/* DDT + ConstructTable (DFJ Prop.1 analogue for SR(7,2,2,6))          */
/* ================================================================== */
/* DDT[Δin][Δout] : list of x with S(x)^S(x^Δin)=Δout.                  */
/* For inversion-type S-box: 0, 2, or 4 solutions, paired {x, x^Δin}.   */
static gf  DDTsol[Q][Q][4];
static uint8_t DDTn[Q][Q];
static void ddt_init(void){
    memset(DDTn,0,sizeof DDTn);
    for(int din=1;din<Q;din++) for(int x=0;x<Q;x++){
        gf dout=(gf)(SBOX[x]^SBOX[x^din]);
        uint8_t *n=&DDTn[din][dout];
        if(*n<4) DDTsol[din][dout][*n]=(gf)x;
        (*n)++;
    }
}

/* 6-parameter tuple + 4 branch bits. */
typedef struct {
    gf dy1;          /* Δy_1[0]               */
    gf x2_0, x2_1;   /* x_2[0,1] of reference */
    gf dw4;          /* Δw_4[0] = Δx_5[0]     */
    gf z4_0, z4_1;   /* z_4[0,1] = S(x_4[0]), S(x_4[3]) of reference */
} params_t;

/* Forward: Δx_3[0..3] from (dy1, x2). */
static void ct_forward_dx3(const params_t *P, gf dx3[4]){
    gf dx2_0=MULTBL[MCa][P->dy1], dx2_1=MULTBL[MCb][P->dy1];
    gf dy2_0=(gf)(SBOX[P->x2_0]^SBOX[P->x2_0^dx2_0]);
    gf dy2_1=(gf)(SBOX[P->x2_1]^SBOX[P->x2_1^dx2_1]);
    /* z_2 active at cells 0 (=y2[0]) and 3 (=y2[1]); MC each column */
    dx3[0]=MULTBL[MCa][dy2_0]; dx3[1]=MULTBL[MCb][dy2_0];
    dx3[2]=MULTBL[MCb][dy2_1]; dx3[3]=MULTBL[MCa][dy2_1];
}
/* Backward: Δy_3[0..3] from (dw4, z4). */
static void ct_backward_dy3(const params_t *P, gf dy3[4], gf *x4_0, gf *x4_3){
    gf dz4_0=MULTBL[MCa][P->dw4], dz4_1=MULTBL[MCb][P->dw4];
    *x4_0=SINV[P->z4_0]; *x4_3=SINV[P->z4_1];
    gf dx4_0=(gf)(*x4_0 ^ SINV[P->z4_0^dz4_0]);
    gf dx4_3=(gf)(*x4_3 ^ SINV[P->z4_1^dz4_1]);
    /* Δx_4[1]=Δx_4[2]=0 (differential).  Δw_3=Δx_4.  Δz_3=MC^{-1}Δw_3.  */
    gf dz3_0=MULTBL[MCa][dx4_0], dz3_1=MULTBL[MCb][dx4_0];
    gf dz3_2=MULTBL[MCb][dx4_3], dz3_3=MULTBL[MCa][dx4_3];
    /* Δy_3 = SR^{-1}(Δz_3)  (swap 1<->3) */
    dy3[0]=dz3_0; dy3[1]=dz3_3; dy3[2]=dz3_2; dy3[3]=dz3_1;
}

/* Given params + 4 branch bits, compute d_ω = Δx_5[0]_ω for all Δu∈GF(64).
 * Returns 0 if DDT infeasible, else 1 (and fills d[64], d[0]=0).
 * Reference is Δu=0 (the pair's first element).                         */
static int construct_d(const params_t *P, int branch, gf d[Q]){
    gf dx3[4], dy3[4], x4_0r, x4_3r;
    ct_forward_dx3(P,dx3);
    ct_backward_dy3(P,dy3,&x4_0r,&x4_3r);
    /* DDT solve x_3[j]_ref.  2 branch bits per cell (covers DDT entries up to 4). */
    gf x3r[4];
    for(int j=0;j<4;j++){
        int n=DDTn[dx3[j]][dy3[j]];
        int b=(branch>>(2*j))&3;
        if(b>=n) return 0;
        x3r[j]=DDTsol[dx3[j]][dy3[j]][b];
    }
    /* derive constants C_j and k_3[0,3] */
    gf Sx20=SBOX[P->x2_0], Sx21=SBOX[P->x2_1];
    gf C0=(gf)(x3r[0]^MULTBL[MCa][Sx20]);
    gf C1=(gf)(x3r[1]^MULTBL[MCb][Sx20]);
    gf C2=(gf)(x3r[2]^MULTBL[MCb][Sx21]);
    gf C3=(gf)(x3r[3]^MULTBL[MCa][Sx21]);
    gf k3_0=(gf)(x4_0r ^ (MULTBL[MCa][SBOX[x3r[0]]]^MULTBL[MCb][SBOX[x3r[3]]]));
    gf k3_3=(gf)(x4_3r ^ (MULTBL[MCb][SBOX[x3r[2]]]^MULTBL[MCa][SBOX[x3r[1]]]));
    gf Sx40r=P->z4_0, Sx43r=P->z4_1;
    /* loop Δu */
    for(int du=0;du<Q;du++){
        gf x20=(gf)(P->x2_0^MULTBL[MCa][du]);
        gf x21=(gf)(P->x2_1^MULTBL[MCb][du]);
        gf s0=SBOX[x20], s1=SBOX[x21];
        gf x30=(gf)(MULTBL[MCa][s0]^C0), x31=(gf)(MULTBL[MCb][s0]^C1);
        gf x32=(gf)(MULTBL[MCb][s1]^C2), x33=(gf)(MULTBL[MCa][s1]^C3);
        gf x40=(gf)(MULTBL[MCa][SBOX[x30]]^MULTBL[MCb][SBOX[x33]]^k3_0);
        gf x43=(gf)(MULTBL[MCb][SBOX[x32]]^MULTBL[MCa][SBOX[x31]]^k3_3);
        d[du]=(gf)(MULTBL[MCa][SBOX[x40]^Sx40r]^MULTBL[MCb][SBOX[x43]^Sx43r]);
    }
    return 1;
}

/* ctable mode: for random keys, find a right pair, extract true params,
 * run construct_d, assert it matches the true d-sequence from internals. */
static void ctable_test(int N){
    int ok=0, nopair=0;
    for(int t=0;t<N;t++){
        gf mk[4]={rnd6(),rnd6(),rnd6(),rnd6()}, rk[NROUNDS+1][4];
        key_schedule(mk,rk);
        /* pick random x_1 with cell 0 = some ref value; we need a PAIR that
           satisfies the 1->2->4->2->1 truncated differential.  Easiest: fix
           x_1[1..3], vary x_1[0] over all 64, encrypt, find ω,ω' with
           Δx_4[1]=Δx_4[2]=0 AND Δw_4[1]=0 (i.e., Δx_5 active only cell 0). */
        gf pts[Q][4]; gf xs_all[Q][NROUNDS][4], ct[Q][4];
        int ref=-1, oth=-1;
        for(int tries=0; tries<512 && ref<0; tries++){
            gf base_x1[4]={0,rnd6(),rnd6(),rnd6()};
            build_delta_set(rk,base_x1,pts);
            for(int i=0;i<Q;i++) encrypt_trace(pts[i],rk,ct[i],xs_all[i]);
            for(int i=0;i<Q && ref<0;i++) for(int j=i+1;j<Q;j++){
                gf *xa=xs_all[i][4], *xb=xs_all[j][4];
                if(xa[1]!=xb[1]||xa[2]!=xb[2]) continue;            /* Δx_4[1,2]=0 */
                gf *x5a=xs_all[i][5], *x5b=xs_all[j][5];
                if(x5a[1]!=x5b[1]||x5a[2]!=x5b[2]||x5a[3]!=x5b[3]) continue; /* Δx_5 only [0] */
                ref=i; oth=j; break;
            }
        }
        if(ref<0){ nopair++; continue; }
        /* extract true params from internals at ref */
        params_t P;
        P.dy1  = (gf)(SBOX[xs_all[ref][1][0]]^SBOX[xs_all[oth][1][0]]);
        P.x2_0 = xs_all[ref][2][0];
        P.x2_1 = xs_all[ref][2][1];
        P.dw4  = (gf)(xs_all[ref][5][0]^xs_all[oth][5][0]);
        P.z4_0 = SBOX[xs_all[ref][4][0]];
        P.z4_1 = SBOX[xs_all[ref][4][3]];
        /* true d[ω] indexed by Δu = S(x_1[0]_ω)^S(x_1[0]_ref) */
        gf s=xs_all[ref][5][0];
        gf dtrue[Q]={0};
        gf Sref=SBOX[xs_all[ref][1][0]];
        for(int w=0;w<Q;w++){
            gf du=(gf)(SBOX[xs_all[w][1][0]]^Sref);
            dtrue[du]=(gf)(s^xs_all[w][5][0]);
        }
        /* true branch bits: which DDT partner is x_3[j]_ref */
        gf dx3[4],dy3[4],x40r,x43r;
        ct_forward_dx3(&P,dx3); ct_backward_dy3(&P,dy3,&x40r,&x43r);
        int branch=0;
        for(int j=0;j<4;j++){
            gf xt=xs_all[ref][3][j]; int n=DDTn[dx3[j]][dy3[j]], found=-1;
            for(int k=0;k<n && k<4;k++) if(DDTsol[dx3[j]][dy3[j]][k]==xt) found=k;
            assert(found>=0);
            branch |= found<<(2*j);
        }
        gf d[Q];
        assert(construct_d(&P,branch,d));
        assert(!memcmp(d,dtrue,Q));
        ok++; (void)oth;
    }
    printf("ConstructTable: true params -> true d-sequence  %d/%d OK  (%d no-right-pair)\n",
           ok,N,nopair);
}

/* ================================================================== */
/* Offline table build: enumerate all (params, branch), insert χ fp    */
/* (both parities) into a Bloom filter backed by an mmap'd file.       */
/* ================================================================== */
/* ---- DFJ-variant multiset hash (u_5 guessed online, plain difference multiset) ---- */
static int DFJ_MODE=0;                 /* 1 = build/attack with multiset hash instead of chi */
static int COLD_MODE=0;                /* 1 = no cross-branch amortization in the table build (naive) */
static uint64_t R64MS[Q];
static void r64_init(void){             /* fixed seed so build and attack runs agree */
    uint64_t x=0x9E3779B97F4A7C15ULL;
    for(int i=0;i<Q;i++){ x^=x<<13; x^=x>>7; x^=x<<17; R64MS[i]=x*0xff51afd7ed558ccdULL; }
}
static inline uint64_t mset_hash(const gf *d, int n){ uint64_t h=0; for(int i=0;i<n;i++) h+=R64MS[d[i]]; return h; }
static uint64_t *BLOOM=NULL;
static uint64_t  BLOOM_MASK=0;    /* BLOOM_BITS-1, BLOOM_BITS a power of 2 */
#ifndef BLOOM_K
#define BLOOM_K 5
#endif

static inline uint64_t mix64(uint64_t x){
    x^=x>>33; x*=0xff51afd7ed558ccdULL;
    x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL;
    x^=x>>33; return x;
}
static inline void bloom_set(uint64_t h){
    uint64_t a=mix64(h), b=mix64(h^0x9e3779b97f4a7c15ULL);
    for(int k=0;k<BLOOM_K;k++){
        uint64_t bit=(a+(uint64_t)k*b)&BLOOM_MASK;
        __atomic_or_fetch(&BLOOM[bit>>6],1ULL<<(bit&63),__ATOMIC_RELAXED);
    }
}
static inline int bloom_test(uint64_t h){
    uint64_t a=mix64(h), b=mix64(h^0x9e3779b97f4a7c15ULL);
    for(int k=0;k<BLOOM_K;k++){
        uint64_t bit=(a+(uint64_t)k*b)&BLOOM_MASK;
        if(!(BLOOM[bit>>6]&(1ULL<<(bit&63)))) return 0;
    }
    return 1;
}
static uint64_t BLOOM_BYTES=0;
static void bloom_map(const char *path, int log2bits, int create){
    uint64_t bits=1ULL<<log2bits; BLOOM_BYTES=bits>>3; BLOOM_MASK=bits-1;
    if(create){
        BLOOM=mmap(NULL,BLOOM_BYTES,PROT_READ|PROT_WRITE,
                   MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE,-1,0);
        if(BLOOM==MAP_FAILED){perror("mmap anon");exit(1);}
        fprintf(stderr,"bloom: anon  %.1f GB  (will write to %s at end)\n",
                BLOOM_BYTES/1073741824.0,path);
    }else{
        int fd=open(path,O_RDONLY); if(fd<0){perror("open");exit(1);}
        BLOOM=mmap(NULL,BLOOM_BYTES,PROT_READ,MAP_PRIVATE,fd,0);
        if(BLOOM==MAP_FAILED){perror("mmap");exit(1);}
        fprintf(stderr,"bloom: %s  %.1f GB  RO\n",path,BLOOM_BYTES/1073741824.0);
    }
}
static void bloom_save(const char *path){
    int fd=open(path,O_WRONLY|O_CREAT|O_TRUNC,0644);
    if(fd<0){perror("open w");exit(1);}
    uint64_t off=0;
    while(off<BLOOM_BYTES){
        ssize_t w=write(fd,(char*)BLOOM+off,BLOOM_BYTES-off>1<<30?1<<30:BLOOM_BYTES-off);
        if(w<=0){perror("write");exit(1);} off+=w;
    }
    close(fd);
    fprintf(stderr,"bloom: saved %.1f GB to %s\n",BLOOM_BYTES/1073741824.0,path);
}

static double wall(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts);
    return ts.tv_sec+1e-9*ts.tv_nsec;
}

/* Enumerate param work-units.  One wu = one (dy1, x2_0) pair.
 * wu_lo..wu_hi in [0, 63*64).  Returns #d-sequences inserted.          */
static long build_table(int wu_lo, int wu_hi, int nthreads){
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#else
    (void)nthreads;
#endif
    long total=0; volatile int done=0;
    double t0=wall();
    #pragma omp parallel for schedule(dynamic,1) reduction(+:total)
    for(int wu=wu_lo; wu<wu_hi; wu++){
        int dy1=1+wu/Q, x2_0=wu%Q;
        long loc=0;
        /* precompute s0[Δu]=S(x2_0^2Δu), s1_row later per x2_1 */
        gf s0[Q]; for(int du=0;du<Q;du++) s0[du]=SBOX[x2_0^MULTBL[MCa][du]];
        gf dx2_0=MULTBL[MCa][dy1], dx2_1=MULTBL[MCb][dy1];
        gf dy2_0=(gf)(SBOX[x2_0]^SBOX[x2_0^dx2_0]);
        gf dx3_0=MULTBL[MCa][dy2_0], dx3_1=MULTBL[MCb][dy2_0];

        for(int x2_1=0;x2_1<Q;x2_1++){
            gf s1[Q]; for(int du=0;du<Q;du++) s1[du]=SBOX[x2_1^MULTBL[MCb][du]];
            gf dy2_1=(gf)(SBOX[x2_1]^SBOX[x2_1^dx2_1]);
            gf dx3_2=MULTBL[MCb][dy2_1], dx3_3=MULTBL[MCa][dy2_1];
            /* precompute y3-terms per branch of each cell */
            /* x3[0]=2*s0+C0, x3[1]=3*s0+C1, x3[2]=3*s1+C2, x3[3]=2*s1+C3 */
            for(int dw4=1;dw4<Q;dw4++){
                gf dz4_0=MULTBL[MCa][dw4], dz4_1=MULTBL[MCb][dw4];
            for(int z4_0=0;z4_0<Q;z4_0++){
                gf x4_0r=SINV[z4_0];
                gf dx4_0=(gf)(x4_0r^SINV[z4_0^dz4_0]);
                gf dy3_0=MULTBL[MCa][dx4_0], dy3_3=MULTBL[MCb][dx4_0];
                int n0=DDTn[dx3_0][dy3_0], n3=DDTn[dx3_3][dy3_3];
                if(!n0||!n3) continue;
            for(int z4_1=0;z4_1<Q;z4_1++){
                gf x4_3r=SINV[z4_1];
                gf dx4_3=(gf)(x4_3r^SINV[z4_1^dz4_1]);
                gf dy3_1=MULTBL[MCa][dx4_3], dy3_2=MULTBL[MCb][dx4_3];
                int n1=DDTn[dx3_1][dy3_1], n2=DDTn[dx3_2][dy3_2];
                if(!n1||!n2) continue;
                gf Sx40r=(gf)z4_0, Sx43r=(gf)z4_1;
                for(int b0=0;b0<n0;b0++){
                  gf x3r0=DDTsol[dx3_0][dy3_0][b0];
                  gf C0=(gf)(x3r0^MULTBL[MCa][s0[0]]);
                for(int b3=0;b3<n3;b3++){
                  gf x3r3=DDTsol[dx3_3][dy3_3][b3];
                  gf C3=(gf)(x3r3^MULTBL[MCa][s1[0]]);
                  gf k3_0=(gf)(x4_0r^(MULTBL[MCa][SBOX[x3r0]]^MULTBL[MCb][SBOX[x3r3]]));
                  /* x4[0](Δu)=2*S(2*s0+C0)^3*S(2*s1+C3)^k3_0 */
                  gf y4a[Q];
                  for(int du=0;du<Q;du++)
                      y4a[du]=SBOX[(MULTBL[MCa][SBOX[MULTBL[MCa][s0[du]]^C0]]
                                   ^MULTBL[MCb][SBOX[MULTBL[MCa][s1[du]]^C3]]^k3_0)];
                for(int b1=0;b1<n1;b1++){
                  gf x3r1=DDTsol[dx3_1][dy3_1][b1];
                  gf C1=(gf)(x3r1^MULTBL[MCb][s0[0]]);
                for(int b2=0;b2<n2;b2++){
                  gf x3r2=DDTsol[dx3_2][dy3_2][b2];
                  gf C2=(gf)(x3r2^MULTBL[MCb][s1[0]]);
                  gf k3_3=(gf)(x4_3r^(MULTBL[MCb][SBOX[x3r2]]^MULTBL[MCa][SBOX[x3r1]]));
                  if(COLD_MODE){ /* naive: rebuild the a-side sequence per branch tuple (no amortization) */
                      for(int du=0;du<Q;du++)
                          y4a[du]=SBOX[(MULTBL[MCa][SBOX[MULTBL[MCa][s0[du]]^C0]]
                                       ^MULTBL[MCb][SBOX[MULTBL[MCa][s1[du]]^C3]]^k3_0)];
                  }
                  gf dinv[Q+1], dseq[Q];
                  for(int du=0;du<Q;du++){
                      gf y4b=SBOX[(MULTBL[MCb][SBOX[MULTBL[MCb][s1[du]]^C2]]
                                  ^MULTBL[MCa][SBOX[MULTBL[MCb][s0[du]]^C1]]^k3_3)];
                      gf d=(gf)(MULTBL[MCa][y4a[du]^Sx40r]^MULTBL[MCb][y4b^Sx43r]);
                      dseq[du]=d; dinv[du]=INV[d];
                  }
                  dinv[Q]=0;
                  if(DFJ_MODE){ bloom_set(mset_hash(dseq,Q)); }
                  else { bloom_set(chi_fp(dinv,Q)); bloom_set(chi_fp(dinv,Q+1)); }
                  loc++;
                }}}}
            }}}
        }
        total+=loc;
        #pragma omp atomic
        done++;
        fprintf(stderr,"[build] wu %4d/%d  %.1fs wall  %ld dseq  (brute-fallbacks: %ld)\n",
                done,wu_hi-wu_lo,wall()-t0,loc,CHI_BRUTE_FALLBACKS);
    }
    return total;
}

static void build_mode(const char *path,int log2bits,int wu_lo,int wu_hi,int nthreads){
    bloom_map(path,log2bits,1);
    double t0=wall();
    long n=build_table(wu_lo,wu_hi,nthreads);
    double tb=wall()-t0;
    uint64_t words=BLOOM_BYTES>>3, pop=0;
    for(uint64_t i=0;i<words;i++) pop+=__builtin_popcountll(BLOOM[i]);
    double load=(double)pop/(BLOOM_MASK+1);
    double fpr=1; for(int k=0;k<BLOOM_K;k++) fpr*=load;
    printf("build: %ld d-seq in %.1fs (%.0f/s), bloom load=%.4f, est FP rate=%.2e\n",
           n,tb,n/tb,load,fpr);
    bloom_save(path);
}

/* ================================================================== */
/* Full online attack (χ fp, DDT-constrained guesses, NO u_5[0]).      */
/*                                                                     */
/* 1. Encrypt 2^12 diagonal structure (P[0,3] vary, P[1,2] fixed).     */
/* 2. Find candidate pairs: ΔC[1]=ΔC[2]=0.                             */
/* 3. Per pair: DDT-constrain k_{-1}[0,3] (Δx_1[1]=0) and k_6[0,3]     */
/*    (Δz_5[1]=0); for each combo build δ-set (ref = pair[0]),         */
/*    compute g, χ fp (both parities), Bloom lookup.                   */
/* 4. On hit: brute mk[1,2], verify via key-schedule + trial encrypt.  */
/* ================================================================== */
static uint8_t DDTinvN[Q][Q];          /* for S^{-1} */
static gf      DDTinvSol[Q][Q][4];
static void ddtinv_init(void){
    memset(DDTinvN,0,sizeof DDTinvN);
    for(int din=1;din<Q;din++) for(int y=0;y<Q;y++){
        gf dout=(gf)(SINV[y]^SINV[y^din]);
        uint8_t *n=&DDTinvN[din][dout];
        if(*n<4) DDTinvSol[din][dout][*n]=(gf)y;
        (*n)++;
    }
}

static void attack_full(const char *bloom_path, int log2bits, unsigned seed_mk,
                        int self_only){
    srand(seed_mk);
    gf MK[4]={rnd6(),rnd6(),rnd6(),rnd6()}, RK[NROUNDS+1][4];
    key_schedule(MK,RK);
    if(!self_only) bloom_map(bloom_path,log2bits,0);
    long Tpairs=0,Tlk=0,Thit=0,Tks=0; int Tstruct=0;
    int recovered=0; gf REC[4]={0};
    double t0=wall();

  for(int str=0; str<16 && !recovered; str++){
    Tstruct++;
    gf Pfix1=rnd6(), Pfix2=rnd6();
    static gf CT[Q][Q][4];            /* CT[P0][P3] */
    gf xs[NROUNDS][4];
    for(int a=0;a<Q;a++) for(int b=0;b<Q;b++){
        gf pt[4]={(gf)a,Pfix1,Pfix2,(gf)b};
        encrypt_trace(pt,RK,CT[a][b],xs);
    }
    gf KP[4]={7,Pfix1,Pfix2,13}, KC[4]; encrypt_trace(KP,RK,KC,xs);
    if(str==0)
      printf("secret key = %02x %02x %02x %02x  (k_{-1}=mk; k_6[0,3]=%d,%d)\n",
             MK[0],MK[1],MK[2],MK[3],RK[NROUNDS][0],RK[NROUNDS][3]);

    /* ---- Bloom source ---- */
    uint64_t SELF_FP[2]={0,0};
    if(self_only){
        /* single-entry "table": compute true offline fp from internals of
           a right pair located via known internals (pipeline validation only). */
        /* we need a right pair inside THIS structure. Search with internals. */
        int refP0=-1,refP3=-1,othP0=-1,othP3=-1;
        for(int a=0;a<Q && refP0<0;a++)for(int b=0;b<Q && refP0<0;b++)
        for(int c=0;c<Q && refP0<0;c++)for(int dd=0;dd<Q;dd++){
            if(a==c&&b==dd) continue;
            gf *C1=CT[a][b],*C2=CT[c][dd];
            if(C1[1]!=C2[1]||C1[2]!=C2[2]) continue;
            /* full-internal check */
            gf p1[4]={(gf)a,Pfix1,Pfix2,(gf)b},p2[4]={(gf)c,Pfix1,Pfix2,(gf)dd};
            gf x1a[NROUNDS][4],x1b[NROUNDS][4],t1[4],t2[4];
            encrypt_trace(p1,RK,t1,x1a);encrypt_trace(p2,RK,t2,x1b);
            if(x1a[1][1]!=x1b[1][1]) continue;                    /* Δx_1[1]=0 */
            if(x1a[5][1]!=x1b[5][1]||x1a[5][2]!=x1b[5][2]||x1a[5][3]!=x1b[5][3]) continue;
            refP0=a;refP3=b;othP0=c;othP3=dd;break;
        }
        if(refP0<0){fprintf(stderr,"  [struct %d: no right pair]\n",str);continue;}
        /* compute true d-seq (δ-set around ref with true k_{-1}) */
        gf zref0=SBOX[refP0^MK[0]], zref3=SBOX[refP3^MK[3]];
        gf x5[Q],s;
        for(int t=0;t<Q;t++){
            gf P0=(gf)(SINV[zref0^MULTBL[MCa][t]]^MK[0]);
            gf P3=(gf)(SINV[zref3^MULTBL[MCb][t]]^MK[3]);
            gf pt[4]={P0,Pfix1,Pfix2,P3},ct[4];
            encrypt_trace(pt,RK,ct,xs); x5[t]=xs[5][0];
        }
        s=x5[0];
        if(!s){fprintf(stderr,"  [struct %d: s=0]\n",str);continue;}
        gf dinv[Q+1]; for(int w=0;w<Q;w++)dinv[w]=gf_inv(s^x5[w]); dinv[Q]=0;
        SELF_FP[0]=chi_fp(dinv,Q); SELF_FP[1]=chi_fp(dinv,Q+1);
        printf("self_only: right pair P=(%d,%d)/(%d,%d)  s=%02x  fp=%016llx/%016llx\n",
               refP0,refP3,othP0,othP3,s,
               (unsigned long long)SELF_FP[0],(unsigned long long)SELF_FP[1]);
    }
    #define LOOKUP(h) (self_only ? ((h)==SELF_FP[0]||(h)==SELF_FP[1]) : bloom_test(h))

    /* ---- online ---- */
    gf T0[Q],T1[Q],GLUT[Q];
    for(int x=0;x<Q;x++){T0[x]=gf_mul(MCa,SINV[x]);T1[x]=gf_mul(MCb,SINV[x]);}
    for(int x=0;x<Q;x++) GLUT[x]=gf_inv(LINV[x]);
    gf inv2=gf_inv(2), r32=gf_mul(3,inv2);

    long npairs=0, nlk=0, nhit=0, nks=0;

    /* candidate-pair search: hash by (C[1],C[2]) */
    for(int a=0;a<Q;a++)for(int b=0;b<Q;b++)
    for(int c=a;c<Q;c++)for(int dd=(c==a?b+1:0);dd<Q;dd++){
        gf *Cr=CT[a][b],*Co=CT[c][dd];
        if(Cr[1]!=Co[1]||Cr[2]!=Co[2]) continue;
        npairs++;
        gf dP0=(gf)(a^c), dP3=(gf)(b^dd);
        if(!dP0||!dP3) continue;   /* need both active for DDT */
        gf dC0=(gf)(Cr[0]^Co[0]), dC3=(gf)(Cr[3]^Co[3]);
        if(!dC0||!dC3) continue;

        /* ---- enumerate k_{-1}[0,3] with Δx_1[1]=0 ---- */
        for(int km0=0;km0<Q;km0++){
            gf dy00=(gf)(SBOX[a^km0]^SBOX[c^km0]);
            gf need_dy03=MULTBL[r32][dy00];
            int nk3=DDTn[dP3][need_dy03];
            for(int ik3=0;ik3<nk3;ik3++){
                gf km3=(gf)(DDTsol[dP3][need_dy03][ik3]^b);  /* x=P3^km3 */
                /* build δ-set ciphertexts (indexed by t) */
                gf zref0=SBOX[a^km0], zref3=SBOX[b^km3];
                gf ct0[Q],ct3[Q];
                for(int t=0;t<Q;t++){
                    gf P0=(gf)(SINV[zref0^MULTBL[MCa][t]]^km0);
                    gf P3=(gf)(SINV[zref3^MULTBL[MCb][t]]^km3);
                    ct0[t]=CT[P0][P3][0]; ct3[t]=CT[P0][P3][3];
                }
                /* ---- enumerate k_6[0,3] with Δz_5[1]=0 ---- */
                for(int k60=0;k60<Q;k60++){
                    gf dx60=(gf)(SINV[Cr[0]^k60]^SINV[Co[0]^k60]);
                    gf need_dx61=MULTBL[r32][dx60];
                    int nk63=DDTinvN[dC3][need_dx61];
                    for(int ik63=0;ik63<nk63;ik63++){
                        gf k63=(gf)(DDTinvSol[dC3][need_dx61][ik63]^Cr[3]);
                        int hit_here=0;
                        if(DFJ_MODE){
                            /* DFJ variant: guess u_5[0] (64x more lookups), plain difference multiset */
                            gf v[Q];
                            for(int w=0;w<Q;w++) v[w]=(gf)(T0[ct0[w]^k60]^T1[ct3[w]^k63]);
                            for(int u5=0;u5<Q;u5++){
                                nlk++;
                                gf a0=SINV[v[0]^u5];
                                gf dd[Q]; for(int w=0;w<Q;w++) dd[w]=(gf)(a0^SINV[v[w]^u5]);
                                if(LOOKUP(mset_hash(dd,Q))){ hit_here=1; break; }
                            }
                            if(!hit_here) continue;
                            nhit++;
                        } else {
                        nlk++;
                        /* compute g, χ fp, lookup */
                        gf v0=(gf)(T0[ct0[0]^k60]^T1[ct3[0]^k63]);
                        gf g[Q+1]; g[0]=0;
                        for(int w=1;w<Q;w++){
                            gf vw=(gf)(T0[ct0[w]^k60]^T1[ct3[w]^k63]);
                            g[w]=GLUT[v0^vw];
                        }
                        g[Q]=0;
                        uint64_t f0=chi_fp(g,Q), f1=chi_fp(g,Q+1);
                        if(!LOOKUP(f0) && !LOOKUP(f1)) continue;
                        nhit++;
                        }
                        /* key solve: mk[0]=km0,mk[3]=km3, brute mk[1,2] */
                        gf cand[4]={(gf)km0,0,0,km3}, rkc[NROUNDS+1][4];
                        if(!getenv("NO_DISPOSE")) for(int m1=0;m1<Q;m1++)for(int m2=0;m2<Q;m2++){
                            cand[1]=(gf)m1;cand[2]=(gf)m2; nks++;
                            key_schedule(cand,rkc);
                            if(rkc[NROUNDS][0]!=k60||rkc[NROUNDS][3]!=k63)continue;
                            gf tc[4]; encrypt_trace(KP,rkc,tc,xs);
                            if(memcmp(tc,KC,4)) continue;
                            gf pt2[4]={41,Pfix1,Pfix2,5},tc2[4],kc2[4];
                            encrypt_trace(pt2,RK,kc2,xs);
                            encrypt_trace(pt2,rkc,tc2,xs);
                            if(!memcmp(tc2,kc2,4)){memcpy(REC,cand,4);recovered++;}
                        }
                    }
                }
            }
        }
    }
    Tpairs+=npairs; Tlk+=nlk; Thit+=nhit; Tks+=nks;
    fprintf(stderr,"  [struct %d] pairs=%ld lk=%ld hits=%ld rec=%d\n",
            str,npairs,nlk,nhit,recovered);
  } /* structures */

    double secs=wall()-t0;
    printf("\n-- attack_full (χ, DDT-constrained, no u_5[0]) --\n");
    printf("  structures      : %d\n",Tstruct);
    printf("  candidate pairs : %ld\n",Tpairs);
    printf("  χ lookups       : %ld\n",Tlk);
    printf("  bloom/self hits : %ld\n",Thit);
    printf("  key-sched tries : %ld\n",Tks);
    printf("  key recovered   : %s  -> %02x %02x %02x %02x  (%d hits)\n",
           recovered?"YES":"NO",REC[0],REC[1],REC[2],REC[3],recovered);
    printf("  time            : %.2f s\n",secs);
    if(!getenv("NO_DISPOSE")) assert(recovered && !memcmp(REC,MK,4));
    #undef LOOKUP
}

/* ------------------------------------------------------------------ */
/* One correctness trial                                               */
/* ------------------------------------------------------------------ */
static int trial(int verbose){
    gf mk[4]={rnd6(),rnd6(),rnd6(),rnd6()}, rk[NROUNDS+1][4];
    key_schedule(mk,rk);
    gf base_x1[4]={rnd6(),rnd6(),rnd6(),rnd6()}, pts[Q][4];
    build_delta_set(rk,base_x1,pts);

    gf ct[Q][4], x5_0[Q]; gf xs[NROUNDS][4];
    for(int i=0;i<Q;i++){ encrypt_trace(pts[i],rk,ct[i],xs); x5_0[i]=xs[5][0]; }

    gf s=x5_0[0];
    if(s==0){ if(verbose) printf("s=0 -> skip (absorbed as %d/%d data factor)\n",Q,Q-1); return -1; }

    gf d[Q],dinv[Q],v[Q],g[Q];
    for(int w=0;w<Q;w++){ d[w]=s^x5_0[w]; dinv[w]=gf_inv(d[w]); }
    for(int w=0;w<Q;w++) v[w]=online_v(ct[w],rk[NROUNDS]);
    for(int w=0;w<Q;w++) g[w]=gf_inv(LINV[v[0]^v[w]]);

    /* (a) bridge identity on good w */
    int bad=0;
    gf s2=gf_mul(s,s);
    for(int w=1;w<Q;w++){
        if(x5_0[w]==0 || x5_0[w]==s){ bad++; continue; }
        assert(g[w]==(gf)(gf_mul(s2,dinv[w])^s));
    }

    /* (c) full delta-set, add-0 parity patch */
    gf g1[Q+1],di1[Q+1];
    memcpy(g1,g,Q);  g1[Q]=0;
    memcpy(di1,dinv,Q); di1[Q]=0;
    fp_t on0=fingerprint_nz(g,Q),    on1=fingerprint_nz(g1,Q+1);
    fp_t of0=fingerprint_nz(dinv,Q), of1=fingerprint_nz(di1,Q+1);
    int m00=fp_eq(on0,of0), m11=fp_eq(on1,of1);
    /* On the predicted parity, P_m(online)=s^{2m}P_m(offline) exactly
       => same zero-pattern => same ncells => fp_eq must hold. */
    fp_t *onp=(bad&1)?&on1:&on0, *ofp=(bad&1)?&of1:&of0;
    if(!((bad&1)?m11:m00) || onp->ncells!=ofp->ncells){
        printf("FAIL key=%02x%02x%02x%02x s=%02x bad=%d m00=%d m11=%d "
               "on.nc=(%d,%d) of.nc=(%d,%d)\n",
               mk[0],mk[1],mk[2],mk[3],s,bad,m00,m11,
               on0.ncells,on1.ncells,of0.ncells,of1.ncells);
        for(int i=0;i<NE;i++){
            gf pgi,pdi; gf pp[RMAX];
            single_psums((bad&1)?g1:g,(bad&1)?Q+1:Q,pp);  pgi=Pm_from_p(pp,E_LIST[i]);
            single_psums((bad&1)?di1:dinv,(bad&1)?Q+1:Q,pp); pdi=Pm_from_p(pp,E_LIST[i]);
            printf("  e=%2d  P(g)=%02x  P(dinv)=%02x  s^{2e}*P(dinv)=%02x\n",
                   E_LIST[i],pgi,pdi,gf_mul(gf_pow(s,2*E_LIST[i]),pdi));
        }
        assert(0);
    }
    assert(m00||m11);
    /* (d) zero-pattern agreement + all cells in GF^x */
    for(int i=0;i<onp->ncells;i++) assert(onp->cell[i]!=0);
    if(onp->ncells==0){ if(verbose) printf("  (degenerate: all P_e=0, skip)\n"); return -2; }

    if(verbose)
        printf("key=%02x%02x%02x%02x  s=%02x  #bad=%d  raw=%d add0=%d  ncells=%d  OK\n",
               mk[0],mk[1],mk[2],mk[3],s,bad,m00,m11,onp->ncells);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Wrong-k_6 sweep                                                     */
/* ------------------------------------------------------------------ */
static void wrong_k6_test(void){
    gf mk[4]={rnd6(),rnd6(),rnd6(),rnd6()}, rk[NROUNDS+1][4];
    key_schedule(mk,rk);
    gf base_x1[4]={rnd6(),rnd6(),rnd6(),rnd6()}, pts[Q][4];
    build_delta_set(rk,base_x1,pts);
    gf ct[Q][4], x5_0[Q]; gf xs[NROUNDS][4];
    for(int i=0;i<Q;i++){ encrypt_trace(pts[i],rk,ct[i],xs); x5_0[i]=xs[5][0]; }
    gf s=x5_0[0]; if(s==0){ wrong_k6_test(); return; }
    gf dinv[Q]; for(int w=0;w<Q;w++) dinv[w]=gf_inv(s^x5_0[w]);
    gf di1[Q+1]; memcpy(di1,dinv,Q); di1[Q]=0;
    fp_t of0=fingerprint_nz(dinv,Q), of1=fingerprint_nz(di1,Q+1);
    if(of0.ncells<NE-1 || of1.ncells<NE-1){ wrong_k6_test(); return; }  /* want full-width fp for the demo */
    const gf *k6=rk[NROUNDS];

    int pos[2]={0,3};
    printf("--- wrong-k6 sweep  key=%02x%02x%02x%02x  true k6[0,3]=(%d,%d) ---\n",
           mk[0],mk[1],mk[2],mk[3],k6[0],k6[3]);

    /* helper */
    #define HIT(kg) ({ \
        gf vv[Q],gg[Q+1]; \
        for(int w=0;w<Q;w++){ vv[w]=online_v(ct[w],kg); } \
        for(int w=0;w<Q;w++){ gg[w]=gf_inv(LINV[vv[0]^vv[w]]); } \
        gg[Q]=0; \
        fp_eq(fingerprint_nz(gg,Q),of0)||fp_eq(fingerprint_nz(gg,Q+1),of1); })

    gf kg[4]; memcpy(kg,k6,4);
    assert(HIT(kg));
    for(int pi=0;pi<2;pi++){
        int p=pos[pi], hits=0, extras=0;
        for(int b=0;b<Q;b++){
            kg[p]=(gf)b;
            if(HIT(kg)){ hits++; if(b!=k6[p]) extras++; }
        }
        kg[p]=k6[p];
        printf("  sweep k6[%d]: true byte hits; %d spurious (24-bit fp => ~%.1e/sweep)\n",
               p,extras, 2.0*Q/((double)QM1*QM1*QM1*QM1));
        assert(hits>=1 && extras<=1);
    }
    /* exhaustive over all 64^2 wrong (k6[0],k6[3]) */
    int wrong_hits=0, Ntot=Q*Q-1;
    for(int a=0;a<Q;a++) for(int b=0;b<Q;b++){
        if(a==k6[0]&&b==k6[3]) continue;
        kg[0]=(gf)a; kg[3]=(gf)b;
        wrong_hits+=HIT(kg);
    }
    printf("  %d wrong (k6[0],k6[3]) pairs: %d spurious hits (expect ~%.3f)\n",
           Ntot,wrong_hits, 2.0*Ntot/((double)QM1*QM1*QM1*QM1));
    assert(wrong_hits<=3);
    #undef HIT
}

/* ------------------------------------------------------------------ */
/* Full key-recovery demo (mode "attack").                             */
/*                                                                     */
/* Online loop:                                                        */
/*   for k_{-1}[0,3] in 64^2:  (= mk[0], mk[3])                         */
/*     build delta-set from the 2^12 ciphertext structure               */
/*     for k_6[0,3] in 64^2:                                            */
/*       compute v,g,fp  (NO u_5[0] guess)  -> lookup                   */
/*   -> 2^24 lookups.  DFJ baseline adds a u_5[0] loop: 2^30 lookups.   */
/*                                                                     */
/* The "offline table" here is the single true-param entry (which the   */
/* real 2^40 DFJ table would contain); see header comment re: sizing.   */
/* On hit, solve mk[1,2] from the key-schedule constraint rk[7][0,3].   */
/* ------------------------------------------------------------------ */
static void attack(unsigned seed_mk){
    /* -------- secret side: pick key, encrypt the 2^12 structure ----- */
    srand(seed_mk);
    gf MK[4]={rnd6(),rnd6(),rnd6(),rnd6()}, RK[NROUNDS+1][4];
    key_schedule(MK,RK);
    gf Pfix1=rnd6(), Pfix2=rnd6();
    static gf CT[Q][Q][4];             /* CT[P0][P3] */
    gf xs[NROUNDS][4];
    for(int a=0;a<Q;a++) for(int b=0;b<Q;b++){
        gf pt[4]={(gf)a,Pfix1,Pfix2,(gf)b};
        encrypt_trace(pt,RK,CT[a][b],xs);
    }
    /* one known pt/ct for final verification */
    gf KP[4]={7,Pfix1,Pfix2,13}, KC[4]; encrypt_trace(KP,RK,KC,xs);

    /* -------- "offline": the true d-sequence fingerprint ------------ */
    /* Build true delta-set (same recipe as build_delta_set_from_guess   */
    /* with the TRUE k_{-1}), read x_5[0] from the trace.                */
    gf x5_0[Q];
    for(int t=0;t<Q;t++){
        gf P0=(gf)(SINV[gf_mul(MCa,(gf)t)]^MK[0]);
        gf P3=(gf)(SINV[gf_mul(MCb,(gf)t)]^MK[3]);
        gf pt[4]={P0,Pfix1,Pfix2,P3};
        encrypt_trace(pt,RK,CT[P0][P3],xs);   /* (already in CT, but need xs) */
        x5_0[t]=xs[5][0];
    }
    gf s=x5_0[0];
    if(s==0){ printf("s=0 for this key, rerolling\n"); attack(seed_mk+1); return; }
    gf dinv[Q+1]; for(int w=0;w<Q;w++) dinv[w]=gf_inv(s^x5_0[w]); dinv[Q]=0;
    fp_t OFF0=fingerprint_nz(dinv,Q), OFF1=fingerprint_nz(dinv,Q+1);
    if(OFF0.ncells<NE-1 && OFF1.ncells<NE-1){
        printf("degenerate offline fp, rerolling\n"); attack(seed_mk+1); return;
    }
    printf("secret key = %02x %02x %02x %02x   true k6[0,3]=(%d,%d)  s=%02x\n",
           MK[0],MK[1],MK[2],MK[3],RK[NROUNDS][0],RK[NROUNDS][3],s);

    /* small LUTs */
    static gf T0[Q],T1[Q],GLUT[Q];
    for(int x=0;x<Q;x++){ T0[x]=gf_mul(MCa,SINV[x]); T1[x]=gf_mul(MCb,SINV[x]); }
    for(int x=0;x<Q;x++) GLUT[x]=gf_inv(LINV[x]);

    /* -------- online: exhaustive over k_{-1}[0,3] x k_6[0,3] --------- */
    long lookups=0, spurious=0;
    int recovered=0; gf REC[4]={0};
    clock_t t0=clock();

    for(int km1_0=0;km1_0<Q;km1_0++){
      fprintf(stderr,"\r[mobius] %2d/64  %.1fs",km1_0,(double)(clock()-t0)/CLOCKS_PER_SEC);
      for(int km1_3=0;km1_3<Q;km1_3++){
        /* delta-set indices into CT */
        gf ct0[Q], ct3[Q];
        for(int t=0;t<Q;t++){
            gf P0=(gf)(SINV[gf_mul(MCa,(gf)t)]^km1_0);
            gf P3=(gf)(SINV[gf_mul(MCb,(gf)t)]^km1_3);
            ct0[t]=CT[P0][P3][0]; ct3[t]=CT[P0][P3][3];
        }
        for(int k6_0=0;k6_0<Q;k6_0++) for(int k6_3=0;k6_3<Q;k6_3++){
            lookups++;
            gf v0=(gf)(T0[ct0[0]^k6_0]^T1[ct3[0]^k6_3]);
            gf g[Q+1]; g[0]=0;
            for(int w=1;w<Q;w++){
                gf vw=(gf)(T0[ct0[w]^k6_0]^T1[ct3[w]^k6_3]);
                g[w]=GLUT[v0^vw];
            }
            g[Q]=0;
            fp_t f0=fingerprint_nz(g,Q);
            fp_t f1=fingerprint_nz(g,Q+1);
            if(!(fp_eq(f0,OFF0)||fp_eq(f1,OFF1))) continue;

            /* hit: candidate (mk[0],mk[3],k6[0],k6[3]). Solve mk[1,2]. */
            int true_hit = (km1_0==MK[0]&&km1_3==MK[3]
                          &&k6_0==RK[NROUNDS][0]&&k6_3==RK[NROUNDS][3]);
            gf cand[4]={(gf)km1_0,0,0,(gf)km1_3}, rkc[NROUNDS+1][4];
            int solved=0;
            for(int m1=0;m1<Q;m1++) for(int m2=0;m2<Q;m2++){
                cand[1]=(gf)m1; cand[2]=(gf)m2;
                key_schedule(cand,rkc);
                if(rkc[NROUNDS][0]!=k6_0||rkc[NROUNDS][3]!=k6_3) continue;
                /* 12-bit filter passed; verify */
                gf tc[4]; encrypt_trace(KP,rkc,tc,xs);
                if(!memcmp(tc,KC,4)){
                    memcpy(REC,cand,4); recovered++; solved=1;
                }
            }
            if(!true_hit && !solved) spurious++;
        }
    }}
    fprintf(stderr,"\r[mobius] done                \n");
    double secs=(double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("\n-- Mobius (no u_5[0]) --\n");
    printf("  lookups       : %ld  (= 64^4 = 2^24)\n",lookups);
    printf("  spurious hits : %ld  (24-bit fp vs 1 offline entry => expect ~%0.1f)\n",
           spurious, 2.0*lookups/((double)QM1*QM1*QM1*QM1));
    printf("  key recovered : %s  -> %02x %02x %02x %02x\n",
           recovered?"YES":"NO",REC[0],REC[1],REC[2],REC[3]);
    printf("  time          : %.1f s\n",secs);
    fflush(stdout);
    assert(recovered && !memcmp(REC,MK,4));

    /* -------- DFJ baseline: same outer loops + inner u_5[0] guess ---- */
    /* Multiset hash: commutative Σ_w R64[d_w]  (no sort).               */
    printf("\n-- DFJ baseline (WITH u_5[0] guess) --\n"); fflush(stdout);
    static uint64_t R64[Q];
    for(int i=0;i<Q;i++) R64[i]=((uint64_t)rand()<<33)^((uint64_t)rand()<<11)^rand();
    uint64_t OFFMS=0; for(int w=0;w<Q;w++) OFFMS+=R64[s^x5_0[w]];

    long dlk=0, dsp=0; int drec=0; gf DREC[4]={0};
    t0=clock();
    for(int km1_0=0;km1_0<Q;km1_0++){
      fprintf(stderr,"\r[dfj]    %2d/64  %.1fs",km1_0,(double)(clock()-t0)/CLOCKS_PER_SEC);
      for(int km1_3=0;km1_3<Q;km1_3++){
        gf ct0[Q], ct3[Q];
        for(int t=0;t<Q;t++){
            gf P0=(gf)(SINV[gf_mul(MCa,(gf)t)]^km1_0);
            gf P3=(gf)(SINV[gf_mul(MCb,(gf)t)]^km1_3);
            ct0[t]=CT[P0][P3][0]; ct3[t]=CT[P0][P3][3];
        }
        for(int k6_0=0;k6_0<Q;k6_0++) for(int k6_3=0;k6_3<Q;k6_3++){
            gf v[Q];
            for(int w=0;w<Q;w++) v[w]=(gf)(T0[ct0[w]^k6_0]^T1[ct3[w]^k6_3]);
            for(int u5=0;u5<Q;u5++){
                dlk++;
                gf a0=SINV[v[0]^u5];
                uint64_t h=0;
                for(int w=0;w<Q;w++) h+=R64[a0^SINV[v[w]^u5]];
                if(h!=OFFMS) continue;
                int th=(km1_0==MK[0]&&km1_3==MK[3]
                      &&k6_0==RK[NROUNDS][0]&&k6_3==RK[NROUNDS][3]);
                gf cand[4]={(gf)km1_0,0,0,(gf)km1_3}, rkc[NROUNDS+1][4];
                int sv=0;
                for(int m1=0;m1<Q;m1++) for(int m2=0;m2<Q;m2++){
                    cand[1]=(gf)m1; cand[2]=(gf)m2; key_schedule(cand,rkc);
                    if(rkc[NROUNDS][0]!=k6_0||rkc[NROUNDS][3]!=k6_3) continue;
                    gf tc[4]; encrypt_trace(KP,rkc,tc,xs);
                    if(!memcmp(tc,KC,4)){ memcpy(DREC,cand,4); drec++; sv=1; }
                }
                if(!th&&!sv) dsp++;
            }
        }
    }}
    fprintf(stderr,"\r[dfj]    done                \n");
    secs=(double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("  lookups       : %ld  (= 64^5 = 2^30)\n",dlk);
    printf("  spurious hits : %ld\n",dsp);
    printf("  key recovered : %s  -> %02x %02x %02x %02x\n",
           drec?"YES":"NO",DREC[0],DREC[1],DREC[2],DREC[3]);
    printf("  time          : %.1f s\n",secs);
    printf("\n==> Mobius: %ld lookups vs DFJ: %ld lookups  (%.0fx fewer)\n",
           lookups,dlk,(double)dlk/lookups);
    assert(drec && !memcmp(DREC,MK,4));
}

/* ------------------------------------------------------------------ */
/* Distribution test                                                   */
/* ------------------------------------------------------------------ */
static void dist_test(int N){
    int hist[NE-1][Q]; memset(hist,0,sizeof hist);
    int skipped=0, cnt=0;
    for(int t=0;t<N;t++){
        gf vals[Q]; vals[0]=0;
        for(int w=1;w<Q;w++) vals[w]=rnd6();
        fp_t f=fingerprint_nz(vals,Q);
        if(f.ncells<NE-1){ skipped++; continue; }
        for(int j=0;j<NE-1;j++) hist[j][f.cell[j]]++;
        cnt++;
    }
    printf("--- distribution: %d samples, %d skipped (<%d nonzero P_e) ---\n",cnt,skipped,NE);
    double exp=(double)cnt/QM1;
    for(int j=0;j<NE-1;j++){
        assert(hist[j][0]==0);
        double chi2=0; int lo=1<<30,hi=0;
        for(int v=1;v<Q;v++){
            double d=hist[j][v]-exp; chi2+=d*d/exp;
            if(hist[j][v]<lo)lo=hist[j][v];
            if(hist[j][v]>hi)hi=hist[j][v];
        }
        printf("  cell[%d]  chi2(%d dof)=%7.1f   min/exp/max = %d/%.1f/%d\n",
               j,QM1-1,chi2,lo,exp,hi);
    }
    /* per-cell collision prob */
    for(int j=0;j<NE-1;j++){
        long coll=0; for(int v=1;v<Q;v++) coll+=(long)hist[j][v]*(hist[j][v]-1)/2;
        double p=(double)coll/((double)cnt*(cnt-1)/2);
        printf("  cell[%d] collision prob = %.5f  (ideal 1/63 = %.5f)\n",j,p,1.0/QM1);
    }
}

/* ------------------------------------------------------------------ */
int main(int argc, char**argv){
    unsigned seed=(unsigned)time(NULL);
    for(int i=1;i<argc;i++) if(!strncmp(argv[i],"seed=",5)) seed=(unsigned)atoi(argv[i]+5);
    srand(seed);
    gf_init(); sbox_init(); pow_init(); chi_init(); ddt_init(); ddtinv_init(); r64_init();
    /* RCON = successive powers of primitive element */
    for(int i=0;i<NROUNDS;i++) RCON[i]=EXP[i];

    /* sanity: MC self-inverse */
    { gf s[4]={rnd6(),rnd6(),rnd6(),rnd6()},t[4]; memcpy(t,s,4);
      mix_cols(t); mix_cols(t); assert(!memcmp(s,t,4)); }
    /* sanity: encrypt/decrypt-by-reencrypt roundtrip via trace */
    /* cross-check fast P_m vs naive */
    { gf v[Q]; v[0]=0; for(int w=1;w<Q;w++) v[w]=rnd6();
      gf p[RMAX]; single_psums(v,Q,p);
      for(int i=0;i<NE;i++) assert(Pm_from_p(p,E_LIST[i])==Pm_naive(v,Q,E_LIST[i]));
      gf v1[Q+1]; memcpy(v1,v,Q); v1[Q]=0;
      single_psums(v1,Q+1,p);
      for(int i=0;i<NE;i++) assert(Pm_from_p(p,E_LIST[i])==Pm_naive(v1,Q+1,E_LIST[i]));
      printf("fast P_m matches naive (|V|=63 and 64)  OK\n"); }

    const char *mode = (argc>1)?argv[1]:"trial";
    if(!strcmp(mode,"attack")){
        attack(seed);
    } else if(!strcmp(mode,"wrong")){
        int N=(argc>2)?atoi(argv[2]):5;
        for(int i=0;i<N;i++) wrong_k6_test();
    } else if(!strcmp(mode,"buildcold")){
        DFJ_MODE=1; COLD_MODE=1;
        const char *path=(argc>2)?argv[2]:"bloom_cold.bin";
        int log2bits=(argc>3)?atoi(argv[3]):24;
        int wu_lo=(argc>4)?atoi(argv[4]):0, wu_hi=(argc>5)?atoi(argv[5]):1, thr=(argc>6)?atoi(argv[6]):1;
        build_mode(path,log2bits,wu_lo,wu_hi,thr);
    } else if(!strcmp(mode,"builddfj")){
        DFJ_MODE=1;
        const char *path=(argc>2)?argv[2]:"bloom_dfj.bin";
        int log2bits=(argc>3)?atoi(argv[3]):24;
        int wu_lo   =(argc>4)?atoi(argv[4]):0;
        int wu_hi   =(argc>5)?atoi(argv[5]):1;
        int thr     =(argc>6)?atoi(argv[6]):1;
        build_mode(path,log2bits,wu_lo,wu_hi,thr);
    } else if(!strcmp(mode,"attack_dfj")){
        DFJ_MODE=1;
        const char *path=(argc>2)?argv[2]:"bloom_dfj.bin";
        int log2bits=(argc>3)?atoi(argv[3]):40;
        attack_full(path,log2bits,seed,0);
    } else if(!strcmp(mode,"build")){
        /* ./sr7226f build <bloomfile> <log2bits> <wu_lo> <wu_hi> [threads]
         * full space: wu_lo=0 wu_hi=4032 (=63*64).                      */
        const char *path=(argc>2)?argv[2]:"bloom.bin";
        int log2bits=(argc>3)?atoi(argv[3]):24;
        int wu_lo   =(argc>4)?atoi(argv[4]):0;
        int wu_hi   =(argc>5)?atoi(argv[5]):1;
        int thr     =(argc>6)?atoi(argv[6]):1;
        build_mode(path,log2bits,wu_lo,wu_hi,thr);
    } else if(!strcmp(mode,"bloomfp")){
        /* ./sr7226f bloomfp <bloomfile> <log2bits> [N]
         * Measure the Bloom's empirical FP rate with (a) uniform-random
         * uint64 queries and (b) chi_fp of uniform-random 63-element lists,
         * so we can separate "Bloom artifact" from "real-table correlation". */
        bloom_map(argv[2],atoi(argv[3]),0);
        long N=(argc>4)?atol(argv[4]):100000000L;
        long hA=0,hB=0; uint64_t x=0x243f6a8885a308d3ULL;
        for(long i=0;i<N;i++){ x=mix64(x+i); hA+=bloom_test(x); }
        for(long i=0;i<N;i++){
            gf v[Q+1]; v[0]=0; for(int w=1;w<Q;w++)v[w]=rnd6(); v[Q]=0;
            hB += bloom_test(chi_fp(v,Q)) | bloom_test(chi_fp(v,Q+1));
        }
        printf("bloomfp: N=%ld  uniform-uint64 p=%.6e  random-chi (2-parity) p=%.6e\n",
               N,(double)hA/N,(double)hB/N);
    } else if(!strcmp(mode,"attack_full")){
        /* ./sr7226f attack_full <bloomfile> <log2bits> [seed] [self] */
        const char *path=(argc>2)?argv[2]:"bloom.bin";
        int log2bits=(argc>3)?atoi(argv[3]):40;
        int self=(argc>5 && !strcmp(argv[5],"self"));
        attack_full(path,log2bits,seed,self);
    } else if(!strcmp(mode,"chi")){
        chi_test((argc>2)?atoi(argv[2]):1000);
    } else if(!strcmp(mode,"ctable")){
        ctable_test((argc>2)?atoi(argv[2]):1000);
    } else if(!strcmp(mode,"dist")){
        dist_test((argc>2)?atoi(argv[2]):200000);
    } else {
        int N=(argc>1)?atoi(argv[1]):200, ok=0,s0=0,deg=0;
        for(int i=0;i<N;i++){ int r=trial(i<10); if(r==-1)s0++; else if(r==-2)deg++; else ok++; }
        printf("\n%d/%d trials OK  (s=0 skips: %d, expect ~%.1f;  all-P_e=0 degen: %d)\n",
               ok,N,s0,N/64.0,deg);
    }
    return 0;
}
