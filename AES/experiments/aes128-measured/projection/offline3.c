// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* offline3.c -- REAL-SCALE three-way OFFLINE per-entry cost benchmark for the
 * 7-round AES-128 attack table build, on genuine Prop-2 parameter bases and
 * genuine 255-element difference sequences  d_w = Delta x_5[0]  (w=1..255).
 *
 *   (A) OURS           : DDT-Gray incremental walk (GFNI/AVX-512, ~1.06
 *                        scalar-equivalent S-box-class lookups per element)
 *                        over the 2^20 branch states of a 10-byte base, then
 *                        the brute-free chi-star fingerprint pair (fp0,fp1)
 *                        via the vector kernel (kernel_v2_vec4).
 *   (B) DFJ-OPTIMIZED  : the SAME DDT-Gray incremental walk, but the entry
 *                        is DFJ's object: a multiset hash of the raw
 *                        difference sequence (order-free sum of 64-bit
 *                        random-table values R64[d_w]; no canonicalization,
 *                        no u_5 absorption).
 *   (C) DFJ-AS-PUBLISHED: per-entry COLD recomputation: the 24-ref state is
 *                        rebuilt from the 10-byte base + the 20 DDT branch
 *                        bits and the full round-2..4 difference trace is
 *                        run for each of the 255 elements (no Gray
 *                        amortization, no per-omega caches), then the same
 *                        multiset hash.  Scalar plain C.
 *
 * All three are streamed and hash-accumulated (nothing stored; nothing
 * optimized away).  Counters count scalar-EQUIVALENT S-box-class lookups
 * (one GFNI 64-lane affine = 64 S-box-class lookups) and other 256-entry
 * table reads (GF-mul / R64 hash) per (entry, element), exactly in the
 * ddt_gray_e2e.c counting style.
 *
 * Modes:
 *   gate           : correctness gate on 2 bases x 8192 states (=16384 >1e4
 *                    entries): A-walk sequence == C-cold sequence (byte and
 *                    hash equal), B's multiset hash (walk) == multiset hash of
 *                    C's cold sequence, A fp (vector) == kernel_full_v2
 *                    (scalar) bit-exact.
 *   bench V k reps : time builder V in {A,B,C} over k bases x 2^20 Gray
 *                    states (entries = k<<20), reps repetitions; rdtscp ticks
 *                    + imul core-clock calibration; median [min,max].
 *
 * Build: gcc -O3 -march=native -I../offline -o offline3 offline3.c -lm
 * Run  : taskset -c <idle core> ./offline3 gate
 *        taskset -c <idle core> ./offline3 bench A 16 5
 */
#include "common.h"   /* real Prop-2 bases, DDT-Gray vector walk, cold_E, kernels */

/* ---------------- counters (scalar-equivalent lookup counts) ---------------- */
static u64 CNT_SB = 0;      /* S-box-class lookups (SBOX/iSBOX evaluations)      */
static u64 CNT_TAB = 0;     /* other 256-entry table reads (GF-mul MUL[], R64) */
static volatile u64 SINK;

/* ---------------- core-clock calibration (no PMU in the VM) ---------------- */
static double calib_ratio(void){   /* core cycles per TSC tick via imul latency 3 */
    u64 best=~0ULL;
    for(int rep=0;rep<200;rep++){
        u64 y=3,t0=rdtsc_now();
        for(long i=0;i<20000;i++){
            __asm__ volatile("imulq %1,%0\n\timulq %1,%0\n\timulq %1,%0\n\timulq %1,%0\n\timulq %1,%0":"+r"(y):"r"(5UL));
        }
        u64 t1=rdtsc_now();
        __asm__ volatile(""::"r"(y));
        if(t1-t0<best)best=t1-t0;
    }
    return (3.0*100000.0)/(double)best;
}
static int cmpd(const void*a,const void*b){double x=*(const double*)a,y=*(const double*)b;return x<y?-1:x>y;}

/* ---------------- DFJ multiset hash of the raw sequence ---------------- */
static u64 R64[256];
static inline u64 mset_hash64(const u8 *E){
    /* order-free: sum over the 255 elements of a random 64-bit table entry */
    u64 h = 0;
    for(int w = 1; w < 256; w++) h += R64[E[w]];
    return h;
}

/* ---------------- (C) COUNTED cold trace, DFJ as published ---------------- *
 * identical arithmetic to cold_E() (rounds 2..4 difference trace from the
 * 24-ref state) but with explicit lookup counters.  Per element: 4 SBOX at
 * y_2[col0] + 16 SBOX at y_3 + 4 SBOX at y_4[diag] = 24 S-box-class
 * lookups, plus 4+16+16+4 = 40 GF-mul (MUL[][]) table reads.  Per entry
 * (amortized over 255 elements): 24 reference S-boxes. */
static inline u8 SBc(u8 x){ CNT_SB++; return SBOX[x]; }
static inline u8 GMc(u8 a, u8 b){ CNT_TAB++; return MUL[a][b]; }
static void cold_E_counted(const ref24_t *R, u8 E[256]){
    u8 sb_x2[4], sb_x3[4][4], sb_x4[4];
    for(int r = 0; r < 4; r++) sb_x2[r] = SBc(R->x2[r]);
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) sb_x3[r][c] = SBc(R->x3[r][c]);
    for(int r = 0; r < 4; r++) sb_x4[r] = SBc(R->x4[r]);
    E[0] = 0;
    for(int w = 1; w < 256; w++){
        u8 dy2[4];
        for(int r = 0; r < 4; r++) dy2[r] = (u8)(SBc(R->x2[r] ^ GMc(MCc[r][0], (u8)w)) ^ sb_x2[r]);
        u8 dy3[4][4];
        for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){
            int rc = (4 - c) & 3;
            u8 dx3 = GMc(MCc[r][rc], dy2[rc]);
            dy3[r][c] = (u8)(SBc(R->x3[r][c] ^ dx3) ^ sb_x3[r][c]);
        }
        u8 dy4[4];
        for(int r = 0; r < 4; r++){
            u8 dx4 = 0;
            for(int rp = 0; rp < 4; rp++) dx4 ^= GMc(MCc[r][rp], dy3[rp][(r + rp) & 3]);
            dy4[r] = (u8)(SBc(R->x4[r] ^ dx4) ^ sb_x4[r]);
        }
        E[w] = (u8)(GMc(2, dy4[0]) ^ GMc(3, dy4[1]) ^ dy4[2] ^ dy4[3]);
    }
}
/* rebuild the entry's 24-ref state from the 10-byte base + Gray branch bits
 * (bits 0..3: x_4[diag r] ^= Dx4[r]; bits 4..19: x_3[rp][cp] ^= Dx3, p=rp+4cp).
 * This is the "cold" per-entry setup (20 byte ops): DFJ's build reads the
 * branch solutions per entry from the DDT tables instead of walking. */
static inline void cold_state(const base_t *B, const ref24_t *R0, u32 g, ref24_t *R){
    *R = *R0;
    for(int r = 0; r < 4; r++) if((g >> r) & 1) R->x4[r] ^= B->Dx4[r];
    for(int p = 0; p < 16; p++) if((g >> (4 + p)) & 1) R->x3[p & 3][p >> 2] ^= B->Dx3[p & 3][p >> 2];
}

/* ---------------- builders over one base (NW Gray states) ----------------
 * No per-entry timing calls (rdtscp is ~25-40 cyc and serializing; it would
 * inflate the measurement).  The whole base is timed from outside; the
 * generate/reduce split is obtained by running the SAME loop with the
 * reduce stage switched off (stage=1) and subtracting medians. */
/* A: walk (counted) + vector chi* fingerprints, 4 at a time */
static u64 build_A(walk_t *W, const base_t *B, u32 NW, int stage){
    u8 ebuf[4][256] __attribute__((aligned(64)));
    const u8 *rows[4][4];
    for(int k = 0; k < 4; k++){ rows[k][0] = ebuf[k]; rows[k][1] = ZROWS; rows[k][2] = ZROWS; rows[k][3] = ZROWS; }
    u64 acc = 0;
    walk_init(W, B);                         CNT_SB += 24 + 255ull * 24;   /* cold start, once per 2^20 */
    int nbuf = 0;
    for(u32 i = 0; i < NW; i++){
        if(i){
            int j = __builtin_ctz(i);
            walk_step(W, j);
            CNT_SB += (j < 4) ? (1 + 255ull) : (1 + 2 * 255ull);   /* incremental Gray step */
        }
        memcpy(ebuf[nbuf], W->E, 256);      /* stage the entry for the 4-wide kernel */
        nbuf++;
        if(nbuf == 4){
            if(stage >= 2){
                u64 f0v[4], f1v[4];
                kernel_v2_vec4(rows, f0v, f1v);          /* brute-free chi* pair x4 */
                for(int k = 0; k < 4; k++) acc ^= f0v[k] + (f1v[k] << 1);
            } else acc ^= ebuf[0][255] + ebuf[1][254] + ebuf[2][253] + ebuf[3][252];
            nbuf = 0;
        }
    }
    return acc;
}
/* B: SAME walk (counted) + multiset hash of the raw sequence (255 R64 reads) */
static u64 build_B(walk_t *W, const base_t *B, u32 NW, int stage){
    u64 acc = 0;
    walk_init(W, B);                         CNT_SB += 24 + 255ull * 24;
    for(u32 i = 0; i < NW; i++){
        if(i){
            int j = __builtin_ctz(i);
            walk_step(W, j);
            CNT_SB += (j < 4) ? (1 + 255ull) : (1 + 2 * 255ull);
        }
        if(stage >= 2){ acc ^= mset_hash64(W->E);   CNT_TAB += 255; }   /* 1 read / element */
        else acc ^= W->E[255];
    }
    return acc;
}
/* C: per-entry cold recomputation from the parameters (no amortization) + multiset hash */
static u64 build_C(const base_t *B, const ref24_t *R0, u32 NW, int stage){
    u8 E[256];
    ref24_t R;
    u64 acc = 0;
    for(u32 i = 0; i < NW; i++){
        u32 g = i ^ (i >> 1);                 /* same entry enumeration order as the Gray walk */
        cold_state(B, R0, g, &R);
        cold_E_counted(&R, E);                 /* 24 SB per element + 24 per entry, cold */
        if(stage >= 2){ acc ^= mset_hash64(E);   CNT_TAB += 255; }
        else acc ^= E[255];
    }
    return acc;
}

/* ===================================================================== */
int main(int argc, char **argv){
    const char *mode = argc > 1 ? argv[1] : "gate";
    aesbench_init();
    rng_t rng; rng_seed2(&rng, 0xA5E50FF1ULL, 0x3333ULL);
    {   /* random 64-bit multiset-hash table (DFJ-style) */
        rng_t r2; rng_seed2(&r2, 0xD51, 0xB6);
        for(int i = 0; i < 256; i++) R64[i] = rng_next(&r2);
    }
    walk_t *W = (walk_t*)aligned_alloc(64, sizeof(walk_t));

    if(!strcmp(mode, "gate")){
        /* ---- correctness gate: 2 bases x 2^13 states = 16384 entries ---- */
        const int NB = 2; const u32 NS = 1u << 13;
        u64 n = 0, ok_seq = 0, ok_hash3 = 0, ok_fp = 0, ok_mset = 0;
        u64 seqh_xor = 0;
        for(int b = 0; b < NB; b++){
            base_t B; base_sample(&B, &rng);
            if(!base_check(&B)){ fprintf(stderr, "FATAL base_check\n"); return 2; }
            /* C's reference base state */
            ref24_t R0;
            for(int r = 0; r < 4; r++) R0.x2[r] = B.bx2[r];
            for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) R0.x3[r][c] = B.bx3[r][c];
            for(int r = 0; r < 4; r++) R0.x4[r] = B.bx4[r];
            walk_init(W, &B);
            for(u32 i = 0; i < NS; i++){
                if(i) walk_step(W, __builtin_ctz(i));
                /* A: the walk sequence */
                const u8 *EA = W->E;
                /* B: the SAME walk, reduced to a multiset hash */
                u64 hB = mset_hash64(EA);
                /* C: cold recomputation from base + gray bits */
                ref24_t R; u8 EC[256];
                cold_state(&B, &R0, i ^ (i >> 1), &R);
                cold_E_counted(&R, EC);
                u64 hC = mset_hash64(EC);
                u64 sA = seq_hash(EA), sC = seq_hash(EC);
                n++;
                if(!memcmp(EA, EC, 256)) ok_seq++;
                if(sA == sC) ok_hash3++;        /* A==B trivially (same array); compare A/B vs C */
                if(hB == hC) ok_mset++;
                /* A fingerprint: vector kernel == scalar reference */
                u64 v0, v1, s0, s1;
                fp_vec(EA, &v0, &v1); fp_ref(EA, &s0, &s1);
                if(v0 == s0 && v1 == s1) ok_fp++;
                seqh_xor ^= sA;
            }
        }
        printf("=== OFFLINE correctness gate (%llu entries = %d bases x 2^13 Gray states) ===\n", (unsigned long long)n, NB);
        printf("[gate o-1] delta-sequence A(walk)==C(cold) bytewise        : %llu/%llu\n", (unsigned long long)ok_seq, (unsigned long long)n);
        printf("[gate o-2] raw-sequence hash A==B==C (FNV64 of d_1..d_255): %llu/%llu  (A,B share the walk array; xor=%016llx)\n",
               (unsigned long long)ok_hash3, (unsigned long long)n, (unsigned long long)seqh_xor);
        printf("[gate o-3] B multiset-hash(walk seq) == hash of C cold seq : %llu/%llu\n", (unsigned long long)ok_mset, (unsigned long long)n);
        printf("[gate o-4] A chi* fp: kernel_v2_vec == kernel_full_v2 (both words, bit-exact): %llu/%llu\n", (unsigned long long)ok_fp, (unsigned long long)n);
        int pass = (ok_seq == n && ok_hash3 == n && ok_mset == n && ok_fp == n);
        printf("=== gate: %s ===\n", pass ? "ALL PASS" : "FAIL");
        return pass ? 0 : 1;
    }

    if(!strcmp(mode, "bench")){
        char V = argc > 2 ? argv[2][0] : 'A';
        int nb = argc > 3 ? atoi(argv[3]) : 16;       /* bases, each 2^20 entries */
        int reps = argc > 4 ? atoi(argv[4]) : 5;
        int wb = argc > 5 ? atoi(argv[5]) : 20;       /* walk bits per base */
        int stage = argc > 6 ? atoi(argv[6]) : 2;    /* 2=full (gen+reduce), 1=generation only */
        const u32 NW = 1u << wb;
        double tsc = tsc_hz();
        /* sample the bases ONCE (untimed; rejection sampling is ~2^16 tries/base) */
        base_t *Bs = (base_t*)malloc(sizeof(base_t) * nb);
        ref24_t *R0s = (ref24_t*)malloc(sizeof(ref24_t) * nb);
        for(int b = 0; b < nb; b++){
            base_sample(&Bs[b], &rng);
            if(!base_check(&Bs[b])){ fprintf(stderr, "FATAL base_check\n"); return 2; }
            for(int r = 0; r < 4; r++) R0s[b].x2[r] = Bs[b].bx2[r];
            for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) R0s[b].x3[r][c] = Bs[b].bx3[r][c];
            for(int r = 0; r < 4; r++) R0s[b].x4[r] = Bs[b].bx4[r];
        }
        u64 entries = (u64)nb * NW;
        printf("=== OFFLINE variant %c stage %d: %d bases x 2^%d states = %llu entries per rep (2^%.2f), %d reps; TSC %.4f GHz ===\n",
               V, stage, nb, wb, (unsigned long long)entries, log2((double)entries), reps, tsc / 1e9);
        double tt[16], rr0[16], rr1[16];
        u64 sb_per = 0, tab_per = 0;
        for(int rep = 0; rep < reps; rep++){
            CNT_SB = 0; CNT_TAB = 0;
            u64 acc = 0;
            double r0 = calib_ratio();
            u64 T0 = rdtsc_now();
            for(int b = 0; b < nb; b++){
                if(V == 'A') acc ^= build_A(W, &Bs[b], NW, stage);
                else if(V == 'B') acc ^= build_B(W, &Bs[b], NW, stage);
                else acc ^= build_C(&Bs[b], &R0s[b], NW, stage);
            }
            u64 T1 = rdtsc_now();
            double r1 = calib_ratio();
            SINK ^= acc;
            tt[rep] = (double)(T1 - T0) / entries;
            rr0[rep] = r0; rr1[rep] = r1;
            sb_per = CNT_SB; tab_per = CNT_TAB;
            printf("  rep %d: ticks/entry %.2f  core/TSC %.3f->%.3f  acc=%016llx\n",
                   rep, tt[rep], r0, r1, (unsigned long long)acc);
            fflush(stdout);
        }
        double ratio = 0; for(int i = 0; i < reps; i++) ratio += (rr0[i] + rr1[i]) * 0.5; ratio /= reps;
        qsort(tt, reps, sizeof(double), cmpd);
        double medt = tt[reps/2];
        double sbpe = (double)sb_per / entries, tabpe = (double)tab_per / entries;
        printf("OFF-%c stage%d  ticks/entry median %.2f [min %.2f max %.2f]\n", V, stage, medt, tt[0], tt[reps-1]);
        printf("OFF-%c stage%d  core-cyc/entry median %.1f  (ratio %.3f)   = %.4f cyc/element\n",
               V, stage, medt * ratio, ratio, medt * ratio / 255.0);
        printf("OFF-%c stage%d  COUNTED: S-box-class lookups/entry %.2f = %.4f per element | other 256-tables/entry %.2f = %.4f per element\n",
               V, stage, sbpe, sbpe / 255.0, tabpe, tabpe / 255.0);
        printf("OFF-%c stage%d  per-core entries/s %.3f M (at median)\n", V, stage, tsc / medt / 1e6);
        free(Bs); free(R0s);
        return 0;
    }
    fprintf(stderr, "modes: gate | bench {A|B|C} nbases reps [walkbits] [stage]\n");
    return 1;
}
