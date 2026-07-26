// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* phi.c -- false-positive / uniformity study of the xi-ordered difference-ratio
 * fingerprint Phi (the paper's appendix "The xi-ordered fingerprint").
 *
 * DEFINITION USED (documented in REPORT.md):
 *   Input: 14 ordered GF(2^8) elements D_1..D_14 (AES field).
 *   Let m := the smallest index in {2,...,14} with D_m != D_1   (normally m=2;
 *   if D_2 == D_1 the relabelling rule swaps in the next distinct element).
 *   Phi(D) := ( (D_j ^ D_1) / (D_m ^ D_1) )  for j in {2,...,14} \ {m}, ascending
 *           = 12 bytes.   If all 14 elements are equal, Phi := 0^12 (degenerate).
 *   Phi is invariant under D -> alpha*D ^ beta (alpha != 0): equalities, hence m,
 *   are preserved and alpha cancels in every ratio.
 *
 * Modes:
 *   phi gate1  n=... seed=...                 : invariance gate
 *   phi gate2  n=... seed=... threads=...     : bridge composition on real data
 *   phi wrong  logn=.. dir=.. threads=.. seed=.. : genuine wrong-candidate Phi's
 *   phi table  logbases=.. walkbits=.. dir=.. threads=.. seed=.. : table-entry Phi's
 *   phi anaw   dir=.. threads=..              : analyze wrong buckets
 *   phi anat   dir=.. threads=..              : analyze table buckets
 *   phi cross  wdir=.. tdir=.. threads=..     : cross-match
 *
 * Build: gcc -O3 -march=native -pthread -o phi phi.c -lm
 */
#include "common.h"
#include <wmmintrin.h>

/* ------------------------------------------------------------------ */
/* Phi                                                                   */
/* ------------------------------------------------------------------ */
/* returns m (0-based index of denominator element, 1 in the normal case),
 * -1 if degenerate (all equal).  fp gets 12 bytes. */
static inline int phi14(const u8 *D, u8 *fp){
    int mi = -1;
    for(int k = 1; k < 14; k++) if(D[k] != D[0]){ mi = k; break; }
    if(mi < 0){ memset(fp, 0, 12); return -1; }
    u8 di = GF_inv[(u8)(D[mi] ^ D[0])];
    int o = 0;
    for(int k = 1; k < 14; k++){
        if(k == mi) continue;
        fp[o++] = gmul((u8)(D[k] ^ D[0]), di);
    }
    return mi;
}

/* ------------------------------------------------------------------ */
/* arg helpers                                                           */
/* ------------------------------------------------------------------ */
static const char *argval(int argc, char **argv, const char *key){
    size_t L = strlen(key);
    for(int i = 2; i < argc; i++) if(!strncmp(argv[i], key, L) && argv[i][L] == '=') return argv[i] + L + 1;
    return NULL;
}
static long long argll(int argc, char **argv, const char *key, long long def){
    const char *v = argval(argc, argv, key); return v ? strtoll(v, 0, 0) : def;
}
static const char *argstr(int argc, char **argv, const char *key, const char *def){
    const char *v = argval(argc, argv, key); return v ? v : def;
}

/* ------------------------------------------------------------------ */
/* AES-NI 7-round encryption (rounds 1..7, k_0 whitening, last round no MC) */
/* ------------------------------------------------------------------ */
typedef struct { __m128i rk[8]; } nik_t;
static inline void nik_set(nik_t *N, const aesctx *ctx){
    for(int r = 0; r < 8; r++) N->rk[r] = _mm_loadu_si128((const __m128i*)ctx->rk[r]);
}
static inline __m128i nik_enc(const nik_t *N, __m128i p){
    __m128i s = _mm_xor_si128(p, N->rk[0]);
    s = _mm_aesenc_si128(s, N->rk[1]);
    s = _mm_aesenc_si128(s, N->rk[2]);
    s = _mm_aesenc_si128(s, N->rk[3]);
    s = _mm_aesenc_si128(s, N->rk[4]);
    s = _mm_aesenc_si128(s, N->rk[5]);
    s = _mm_aesenc_si128(s, N->rk[6]);
    s = _mm_aesenclast_si128(s, N->rk[7]);
    return s;
}
static int nik_selftest(void){
    rng_t R; rng_seed2(&R, 0x1234567, 99);
    for(int t = 0; t < 2000; t++){
        u8 key[16], P[16], C1[16], C2[16];
        for(int i = 0; i < 16; i++){ key[i] = rng8(&R); P[i] = rng8(&R); }
        aesctx ctx; aes_set_key(&ctx, key, 7);
        aes_enc(&ctx, P, C1);
        nik_t N; nik_set(&N, &ctx);
        __m128i c = nik_enc(&N, _mm_loadu_si128((const __m128i*)P));
        _mm_storeu_si128((__m128i*)C2, c);
        if(memcmp(C1, C2, 16)) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* generic bucket writer (byte records)                                 */
/* ------------------------------------------------------------------ */
typedef struct {
    int fd[256];
    pthread_mutex_t mtx[256];
    u64 count[256];
    int recsz;
} gbuck_t;
static int gbuck_open(gbuck_t *S, const char *dir, const char *pfx, int recsz){
    char path[700];
    S->recsz = recsz;
    for(int b = 0; b < 256; b++){
        snprintf(path, sizeof path, "%s/%s%02x.bin", dir, pfx, b);
        S->fd[b] = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if(S->fd[b] < 0){ perror(path); return -1; }
        pthread_mutex_init(&S->mtx[b], NULL);
        S->count[b] = 0;
    }
    return 0;
}
static void gbuck_write(gbuck_t *S, int b, const void *buf, size_t len){
    pthread_mutex_lock(&S->mtx[b]);
    const char *p = (const char*)buf; size_t left = len;
    while(left){
        ssize_t w = write(S->fd[b], p, left);
        if(w <= 0){ perror("write"); abort(); }
        p += w; left -= (size_t)w;
    }
    S->count[b] += len / (size_t)S->recsz;
    pthread_mutex_unlock(&S->mtx[b]);
}
static void gbuck_close(gbuck_t *S){
    for(int b = 0; b < 256; b++){ fsync(S->fd[b]); close(S->fd[b]); }
}
typedef struct {
    u8 *buf[256];
    int n[256], cap;
    gbuck_t *S;
} gtl_t;
static void gtl_init(gtl_t *T, gbuck_t *S, int cap){
    T->S = S; T->cap = cap;
    for(int b = 0; b < 256; b++){ T->buf[b] = (u8*)malloc((size_t)cap * S->recsz); T->n[b] = 0; }
}
static inline void gtl_put(gtl_t *T, int b, const void *rec){
    memcpy(T->buf[b] + (size_t)T->n[b] * T->S->recsz, rec, T->S->recsz);
    if(++T->n[b] == T->cap){ gbuck_write(T->S, b, T->buf[b], (size_t)T->n[b] * T->S->recsz); T->n[b] = 0; }
}
static void gtl_flush(gtl_t *T){
    for(int b = 0; b < 256; b++) if(T->n[b]){ gbuck_write(T->S, b, T->buf[b], (size_t)T->n[b] * T->S->recsz); T->n[b] = 0; }
    for(int b = 0; b < 256; b++) free(T->buf[b]);
}

/* ------------------------------------------------------------------ */
/* statistics accumulator                                                */
/* ------------------------------------------------------------------ */
#define NJOINT 4
static const int JOINT_A[NJOINT] = {0, 0, 5, 3};
static const int JOINT_B[NJOINT] = {1, 11, 6, 9};
typedef struct {
    u64 hist[12][256];
    u32 *joint;          /* NJOINT * 65536 */
    u64 n, relabel, degenerate, zero_in_window, a_zero;
} stats_t;
static void stats_init(stats_t *S){
    memset(S, 0, sizeof *S);
    S->joint = (u32*)calloc((size_t)NJOINT * 65536, sizeof(u32));
}
static inline void stats_add(stats_t *S, const u8 *fp){
    for(int k = 0; k < 12; k++) S->hist[k][fp[k]]++;
    for(int q = 0; q < NJOINT; q++) S->joint[q * 65536 + ((int)fp[JOINT_A[q]] << 8 | fp[JOINT_B[q]])]++;
    S->n++;
}
static void stats_merge(stats_t *dst, const stats_t *src){
    for(int k = 0; k < 12; k++) for(int v = 0; v < 256; v++) dst->hist[k][v] += src->hist[k][v];
    for(size_t i = 0; i < (size_t)NJOINT * 65536; i++) dst->joint[i] += src->joint[i];
    dst->n += src->n; dst->relabel += src->relabel; dst->degenerate += src->degenerate;
    dst->zero_in_window += src->zero_in_window; dst->a_zero += src->a_zero;
}
static void stats_write(const stats_t *S, const char *path){
    FILE *f = fopen(path, "w"); if(!f){ perror(path); return; }
    fprintf(f, "# n=%llu relabel=%llu degenerate=%llu zero_in_window=%llu a_zero=%llu\n",
            (unsigned long long)S->n, (unsigned long long)S->relabel, (unsigned long long)S->degenerate,
            (unsigned long long)S->zero_in_window, (unsigned long long)S->a_zero);
    for(int k = 0; k < 12; k++){
        fprintf(f, "hist %d", k);
        for(int v = 0; v < 256; v++) fprintf(f, " %llu", (unsigned long long)S->hist[k][v]);
        fprintf(f, "\n");
    }
    for(int q = 0; q < NJOINT; q++){
        fprintf(f, "joint %d %d", JOINT_A[q], JOINT_B[q]);
        for(int i = 0; i < 65536; i++) fprintf(f, " %u", S->joint[q * 65536 + i]);
        fprintf(f, "\n");
    }
    fclose(f);
}

/* ------------------------------------------------------------------ */
/* GATE 1: invariance                                                    */
/* ------------------------------------------------------------------ */
static int mode_gate1(int argc, char **argv){
    long long n = argll(argc, argv, "n", 1000000);
    u64 seed = (u64)argll(argc, argv, "seed", 0x51ed1);
    rng_t R; rng_seed2(&R, seed, 1);
    u64 fail = 0, degen = 0, relabel = 0, forced = 0;
    for(long long t = 0; t < n; t++){
        u8 D[14], E[14], f1[12], f2[12];
        for(int i = 0; i < 14; i++) D[i] = rng8(&R);
        /* force repeats in 25% of vectors to exercise the relabelling rule */
        u8 ctl = rng8(&R);
        if(ctl < 64){
            forced++;
            int r = 1 + (ctl & 7);          /* copy D_1 into the next r slots */
            for(int i = 1; i <= r && i < 14; i++) D[i] = D[0];
            if(ctl & 8){ D[5] = D[3]; D[9] = D[3]; }   /* some extra internal ties */
        }
        u8 alpha = rng8nz(&R), beta = rng8(&R);
        for(int i = 0; i < 14; i++) E[i] = (u8)(gmul(alpha, D[i]) ^ beta);
        int m1 = phi14(D, f1), m2 = phi14(E, f2);
        if(m1 < 0) degen++;
        if(m1 != 1) relabel++;
        if(m1 != m2 || memcmp(f1, f2, 12)) fail++;
    }
    printf("[gate1] n=%lld failures=%llu forced_tie_vectors=%llu relabel(m!=2)=%llu degenerate(all-equal)=%llu\n",
           n, (unsigned long long)fail, (unsigned long long)forced, (unsigned long long)relabel, (unsigned long long)degen);
    printf("[gate1] %s\n", fail == 0 ? "PASS" : "FAIL");
    return fail == 0 ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* GATE 2: bridge composition on real right-pair data                   */
/* ------------------------------------------------------------------ */
typedef struct {
    int tid, nthr; long long n; u64 seed;
    u64 inst, clean, clean_match, bad_any, bad_d0, bad_da, bad_both, bad_match, a_zero, a_zero_match;
    u64 label_d0, label_da, relabel_on, relabel_off;
} g2w_t;

static void *gate2_worker(void *arg){
    g2w_t *w = (g2w_t*)arg;
    rng_t R; rng_seed2(&R, w->seed, 777ULL * (w->tid + 1));
    for(long long t = w->tid; t < w->n; t += w->nthr){
        u8 key[16], P[16], C[16];
        for(int i = 0; i < 16; i++){ key[i] = rng8(&R); P[i] = rng8(&R); }
        aesctx ctx; aes_set_key(&ctx, key, 7);
        trace_t T; aes_enc_trace(&ctx, P, C, &T);
        /* offline: true parameter refs -> e_true[om] = Delta x_6[0], om = Delta z_2[0] (eta label) */
        refs24_t Rr;
        for(int i = 0; i < 4; i++) Rr.x3c[i] = T.x[3][i];
        for(int i = 0; i < 16; i++) Rr.x4[i] = T.x[4][i];
        for(int i = 0; i < 4; i++) Rr.x5d[i] = T.x[5][DIAG0[i]];
        u8 e_true[256]; prop2_all(&Rr, e_true);
        u8 a = T.x[6][0];
        u8 Doff[14];
        int bad_d0 = 0, bad_da = 0;
        for(int j = 1; j <= 14; j++){
            u8 d = e_true[j];
            Doff[j-1] = GF_inv[d];
            if(d == 0) bad_d0++;
            else if(d == a) bad_da++;      /* if a==0 this is the same event */
        }
        /* online with TRUE key bytes: delta-set via true km1, peel via true k7a */
        u8 km1[4], k7a[4];
        for(int i = 0; i < 4; i++) km1[i] = ctx.rk[0][DIAG0[i]];
        for(int i = 0; i < 4; i++) k7a[i] = ctx.rk[7][ADIAG[i]];
        u8 z1r[4], w1r[4];
        for(int i = 0; i < 4; i++) z1r[i] = SBOX[(u8)(P[DIAG0[i]] ^ km1[i])];
        for(int r = 0; r < 4; r++) w1r[r] = (u8)(gmul(MCc[r][0], z1r[0]) ^ gmul(MCc[r][1], z1r[1]) ^ gmul(MCc[r][2], z1r[2]) ^ gmul(MCc[r][3], z1r[3]));
        u8 xi = (u8)(w1r[0] ^ ctx.rk[1][0]);            /* = x_2^{(0)}[0] */
        if(xi != T.x[2][0]){ fprintf(stderr, "FATAL: xi mismatch\n"); exit(9); }
        nik_t N; nik_set(&N, &ctx);
        u8 x7ref[4];
        {
            for(int i = 0; i < 4; i++) x7ref[i] = iSBOX[(u8)(C[ADIAG[i]] ^ k7a[i])];
        }
        u8 Don[14];
        for(int j = 1; j <= 14; j++){
            u8 v = (u8)(iSBOX[(u8)(SBOX[xi] ^ j)] ^ xi);     /* eta(v)=j */
            u8 w1[4] = { (u8)(w1r[0] ^ v), w1r[1], w1r[2], w1r[3] };
            u8 z1[4];
            for(int r = 0; r < 4; r++) z1[r] = (u8)(gmul(iMCc[r][0], w1[0]) ^ gmul(iMCc[r][1], w1[1]) ^ gmul(iMCc[r][2], w1[2]) ^ gmul(iMCc[r][3], w1[3]));
            u8 PP[16]; memcpy(PP, P, 16);
            for(int i = 0; i < 4; i++) PP[DIAG0[i]] = (u8)(iSBOX[z1[i]] ^ km1[i]);
            u8 Cv[16];
            _mm_storeu_si128((__m128i*)Cv, nik_enc(&N, _mm_loadu_si128((const __m128i*)PP)));
            u8 dy6 = 0;
            for(int i = 0; i < 4; i++) dy6 ^= gmul(iMCc[0][i], (u8)(iSBOX[(u8)(Cv[ADIAG[i]] ^ k7a[i])] ^ x7ref[i]));
            Don[j-1] = GF_inv[Linv[dy6]];                     /* g_v */
        }
        u8 fon[12], foff[12];
        int mon = phi14(Don, fon), moff = phi14(Doff, foff);
        int match = (mon == moff) && !memcmp(fon, foff, 12);
        if(moff != 1) w->relabel_off++;
        if(mon != 1) w->relabel_on++;
        w->inst++;
        w->label_d0 += (u64)bad_d0; w->label_da += (u64)bad_da;
        if(a == 0){ w->a_zero++; if(match) w->a_zero_match++; }
        if(bad_d0 == 0 && bad_da == 0 && a != 0){
            w->clean++; if(match) w->clean_match++;
        } else {
            w->bad_any++;
            if(bad_d0 && bad_da) w->bad_both++; else if(bad_d0) w->bad_d0++; else if(bad_da) w->bad_da++;
            if(match) w->bad_match++;
        }
    }
    return NULL;
}

static int mode_gate2(int argc, char **argv){
    long long n = argll(argc, argv, "n", 20000);
    int nthr = (int)argll(argc, argv, "threads", 8);
    u64 seed = (u64)argll(argc, argv, "seed", 0xb21d6e);
    g2w_t *ws = (g2w_t*)calloc((size_t)nthr, sizeof(g2w_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    for(int i = 0; i < nthr; i++){ ws[i].tid = i; ws[i].nthr = nthr; ws[i].n = n; ws[i].seed = seed; pthread_create(&th[i], NULL, gate2_worker, &ws[i]); }
    g2w_t S; memset(&S, 0, sizeof S);
    for(int i = 0; i < nthr; i++){
        pthread_join(th[i], NULL);
        S.inst += ws[i].inst; S.clean += ws[i].clean; S.clean_match += ws[i].clean_match;
        S.bad_any += ws[i].bad_any; S.bad_d0 += ws[i].bad_d0; S.bad_da += ws[i].bad_da;
        S.bad_both += ws[i].bad_both; S.bad_match += ws[i].bad_match;
        S.a_zero += ws[i].a_zero; S.a_zero_match += ws[i].a_zero_match;
        S.label_d0 += ws[i].label_d0; S.label_da += ws[i].label_da;
        S.relabel_on += ws[i].relabel_on; S.relabel_off += ws[i].relabel_off;
    }
    double pn = (double)S.inst;
    printf("[gate2] instances=%llu (real AES-128 7r, true key, delta-set at x_2[0], match byte x_6[0], labels eta=1..14)\n",
           (unsigned long long)S.inst);
    printf("[gate2] clean instances (no bad label, a!=0): %llu (%.4f); online==offline Phi on clean: %llu/%llu  %s\n",
           (unsigned long long)S.clean, S.clean / pn, (unsigned long long)S.clean_match, (unsigned long long)S.clean,
           S.clean_match == S.clean ? "PASS" : "FAIL");
    printf("[gate2] exception instances: %llu (%.4f)  [d=0 only: %llu, d=a only: %llu, both: %llu; a(=s)=0 instances: %llu]\n",
           (unsigned long long)S.bad_any, S.bad_any / pn, (unsigned long long)S.bad_d0, (unsigned long long)S.bad_da,
           (unsigned long long)S.bad_both, (unsigned long long)S.a_zero);
    printf("[gate2] per-label rates: d=0: %.5f  d=a: %.5f (each ideal 1/256=%.5f)\n",
           S.label_d0 / pn / 14.0, S.label_da / pn / 14.0, 1.0 / 256.0);
    printf("[gate2] Phi still matches on exception instances: %llu/%llu ; on a=0 instances: %llu/%llu\n",
           (unsigned long long)S.bad_match, (unsigned long long)S.bad_any,
           (unsigned long long)S.a_zero_match, (unsigned long long)S.a_zero);
    printf("[gate2] relabel events (m!=2): online %llu, offline %llu\n",
           (unsigned long long)S.relabel_on, (unsigned long long)S.relabel_off);
    printf("[gate2] reference: (1-2/256)^14=%.4f coverage -> exc %.4f ; (1-1/256)^14=%.4f -> exc %.4f ; appendix eta_cov=0.943 -> exc 0.057\n",
           pow(254.0/256, 14), 1 - pow(254.0/256, 14), pow(255.0/256, 14), 1 - pow(255.0/256, 14));
    int ok = (S.clean_match == S.clean);
    printf("[gate2] %s\n", ok ? "PASS (0 mismatches on clean instances)" : "FAIL");
    free(ws); free(th);
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* COVERAGE (plain + tau-variant) on real right-pair data, 16 labels     */
/* ------------------------------------------------------------------ */
/* Per real right-pair instance (random key + plaintext, true key known to the
 * harness): s := x_6[0]; label j (=1..16) is BAD iff d_j in {0, s} (the bridge
 * g = s^2 d^-1 xor s has no pre-image / collides there).
 * Rules measured (all also require s != 0):
 *   plain : window {1..14} bad-free.
 *   tauC  : adjacent-pair-deletion family, 15 probes: W_t = {1..16}\{t,t+1},
 *           probe order t=15 (clean window = labels 1..14) then t=1..14;
 *           covered iff some W_t is bad-free.  ("clean key + 14 substituted keys")
 *   tauA  : spare-substitution family, 15 probes: W_0 = {1..14},
 *           W_j = ({1..14}\{j}) U {15}, j = 1..14; covered iff some W bad-free.
 *   ideal : any 14-of-16, i.e. <= 2 bad labels among 1..16.
 * Gate (bridge composition): on every covered instance the online Phi over the
 * chosen window (computed through the true key bytes) must equal the offline Phi
 * over the same window.                                                    */
typedef struct {
    int tid, nthr; long long n; u64 seed;
    u64 inst, s_zero;
    u64 lab_d0[16], lab_ds[16];          /* per-label bad categories (d=0 ; d=s with s!=0, d!=0) */
    u64 hb14[15], hb16[17];              /* #bad among first 14 / 16 (all instances) */
    u64 hb14_s[15], hb16_s[17];          /* same, restricted to s!=0 */
    u64 cov_plain, m_plain;              /* covered / Phi match on covered */
    u64 cov_tauC, m_tauC, probe_tauC[16];/* probe index used (1=clean first) */
    u64 cov_tauA, m_tauA;
    u64 cov_ideal, m_ideal;
    u64 relabel_any;
} g3w_t;

static inline int phi_sel(const u8 *D16, const int *idx, u8 *fp){
    u8 t[14];
    for(int k = 0; k < 14; k++) t[k] = D16[idx[k]];
    return phi14(t, fp);
}

static void *cov_worker(void *arg){
    g3w_t *w = (g3w_t*)arg;
    rng_t R; rng_seed2(&R, w->seed, 9901ULL * (u64)(w->tid + 1) + 3);
    for(long long t = w->tid; t < w->n; t += w->nthr){
        u8 key[16], P[16], C[16];
        for(int i = 0; i < 16; i++){ key[i] = rng8(&R); P[i] = rng8(&R); }
        aesctx ctx; aes_set_key(&ctx, key, 7);
        trace_t T; aes_enc_trace(&ctx, P, C, &T);
        refs24_t Rr;
        for(int i = 0; i < 4; i++) Rr.x3c[i] = T.x[3][i];
        for(int i = 0; i < 16; i++) Rr.x4[i] = T.x[4][i];
        for(int i = 0; i < 4; i++) Rr.x5d[i] = T.x[5][DIAG0[i]];
        u8 e_true[256]; prop2_all(&Rr, e_true);
        u8 s = T.x[6][0];
        u8 Doff[16], Don[16];
        int bad[16]; int nb14 = 0, nb16 = 0;
        for(int j = 1; j <= 16; j++){
            u8 d = e_true[j];
            Doff[j-1] = GF_inv[d];
            int b = 0;
            if(d == 0){ b = 1; w->lab_d0[j-1]++; }
            else if(d == s){ b = 1; w->lab_ds[j-1]++; }   /* s != 0 here since d != 0 */
            bad[j-1] = b;
            if(b){ nb16++; if(j <= 14) nb14++; }
        }
        w->inst++;
        w->hb14[nb14]++; w->hb16[nb16]++;
        if(s == 0) w->s_zero++; else { w->hb14_s[nb14]++; w->hb16_s[nb16]++; }
        /* online g_v for labels 1..16 through the true key bytes */
        u8 km1[4], k7a[4];
        for(int i = 0; i < 4; i++) km1[i] = ctx.rk[0][DIAG0[i]];
        for(int i = 0; i < 4; i++) k7a[i] = ctx.rk[7][ADIAG[i]];
        u8 z1r[4], w1r[4];
        for(int i = 0; i < 4; i++) z1r[i] = SBOX[(u8)(P[DIAG0[i]] ^ km1[i])];
        for(int r = 0; r < 4; r++) w1r[r] = (u8)(gmul(MCc[r][0], z1r[0]) ^ gmul(MCc[r][1], z1r[1]) ^ gmul(MCc[r][2], z1r[2]) ^ gmul(MCc[r][3], z1r[3]));
        u8 xi = (u8)(w1r[0] ^ ctx.rk[1][0]);
        if(xi != T.x[2][0]){ fprintf(stderr, "FATAL: xi mismatch\n"); exit(9); }
        nik_t N; nik_set(&N, &ctx);
        u8 x7ref[4];
        for(int i = 0; i < 4; i++) x7ref[i] = iSBOX[(u8)(C[ADIAG[i]] ^ k7a[i])];
        for(int j = 1; j <= 16; j++){
            u8 v = (u8)(iSBOX[(u8)(SBOX[xi] ^ j)] ^ xi);
            u8 w1[4] = { (u8)(w1r[0] ^ v), w1r[1], w1r[2], w1r[3] };
            u8 z1[4];
            for(int r = 0; r < 4; r++) z1[r] = (u8)(gmul(iMCc[r][0], w1[0]) ^ gmul(iMCc[r][1], w1[1]) ^ gmul(iMCc[r][2], w1[2]) ^ gmul(iMCc[r][3], w1[3]));
            u8 PP[16]; memcpy(PP, P, 16);
            for(int i = 0; i < 4; i++) PP[DIAG0[i]] = (u8)(iSBOX[z1[i]] ^ km1[i]);
            u8 Cv[16];
            _mm_storeu_si128((__m128i*)Cv, nik_enc(&N, _mm_loadu_si128((const __m128i*)PP)));
            u8 dy6 = 0;
            for(int i = 0; i < 4; i++) dy6 ^= gmul(iMCc[0][i], (u8)(iSBOX[(u8)(Cv[ADIAG[i]] ^ k7a[i])] ^ x7ref[i]));
            Don[j-1] = GF_inv[Linv[dy6]];
        }
        u8 fon[12], foff[12]; int mon, moff; int idx[14];
        /* ---- plain: labels 1..14 ---- */
        if(s != 0 && nb14 == 0){
            w->cov_plain++;
            for(int k = 0; k < 14; k++) idx[k] = k;
            mon = phi_sel(Don, idx, fon); moff = phi_sel(Doff, idx, foff);
            if(mon == moff && !memcmp(fon, foff, 12)) w->m_plain++;
            if(moff != 1) w->relabel_any++;
        }
        if(s != 0){
            /* ---- tauC: adjacent-pair deletion, probe order t=15 then 1..14 ---- */
            static const int order[15] = {15,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
            for(int q = 0; q < 15; q++){
                int td = order[q];            /* delete labels td, td+1 (1-based) */
                int clean = 1;
                for(int j = 1; j <= 16 && clean; j++) if(j != td && j != td+1 && bad[j-1]) clean = 0;
                if(!clean) continue;
                int kk = 0;
                for(int j = 1; j <= 16; j++) if(j != td && j != td+1) idx[kk++] = j-1;
                mon = phi_sel(Don, idx, fon); moff = phi_sel(Doff, idx, foff);
                w->cov_tauC++; w->probe_tauC[q]++;
                if(mon == moff && !memcmp(fon, foff, 12)) w->m_tauC++;
                break;
            }
            /* ---- tauA: W_0={1..14}; W_j = {1..14}\{j} U {15} ---- */
            {
                int found = -1;
                if(nb14 == 0) found = 0;
                else if(nb14 == 1 && !bad[14]){
                    for(int j = 1; j <= 14; j++) if(bad[j-1]){ found = j; break; }
                }
                if(found >= 0){
                    int kk = 0;
                    for(int j = 1; j <= 15; j++){ if(j == found) continue; if(found == 0 && j == 15) continue; idx[kk++] = j-1; }
                    mon = phi_sel(Don, idx, fon); moff = phi_sel(Doff, idx, foff);
                    w->cov_tauA++;
                    if(mon == moff && !memcmp(fon, foff, 12)) w->m_tauA++;
                }
            }
            /* ---- idealized: any 14 of 16 bad-free  <=>  nb16 <= 2 ---- */
            if(nb16 <= 2){
                int kk = 0;
                for(int j = 1; j <= 16 && kk < 14; j++) if(!bad[j-1]) idx[kk++] = j-1;
                mon = phi_sel(Don, idx, fon); moff = phi_sel(Doff, idx, foff);
                w->cov_ideal++;
                if(mon == moff && !memcmp(fon, foff, 12)) w->m_ideal++;
            }
        }
    }
    return NULL;
}

static int mode_cov(int argc, char **argv){
    long long n = argll(argc, argv, "n", 100000);
    int nthr = (int)argll(argc, argv, "threads", 8);
    u64 seed = (u64)argll(argc, argv, "seed", 0xc0be7a6e);
    g3w_t *ws = (g3w_t*)calloc((size_t)nthr, sizeof(g3w_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    for(int i = 0; i < nthr; i++){ ws[i].tid = i; ws[i].nthr = nthr; ws[i].n = n; ws[i].seed = seed; pthread_create(&th[i], NULL, cov_worker, &ws[i]); }
    g3w_t S; memset(&S, 0, sizeof S);
    for(int i = 0; i < nthr; i++){
        pthread_join(th[i], NULL);
        S.inst += ws[i].inst; S.s_zero += ws[i].s_zero;
        for(int j = 0; j < 16; j++){ S.lab_d0[j] += ws[i].lab_d0[j]; S.lab_ds[j] += ws[i].lab_ds[j]; S.probe_tauC[j] += ws[i].probe_tauC[j]; }
        for(int k = 0; k < 15; k++){ S.hb14[k] += ws[i].hb14[k]; S.hb14_s[k] += ws[i].hb14_s[k]; }
        for(int k = 0; k < 17; k++){ S.hb16[k] += ws[i].hb16[k]; S.hb16_s[k] += ws[i].hb16_s[k]; }
        S.cov_plain += ws[i].cov_plain; S.m_plain += ws[i].m_plain;
        S.cov_tauC += ws[i].cov_tauC;   S.m_tauC += ws[i].m_tauC;
        S.cov_tauA += ws[i].cov_tauA;   S.m_tauA += ws[i].m_tauA;
        S.cov_ideal += ws[i].cov_ideal; S.m_ideal += ws[i].m_ideal;
        S.relabel_any += ws[i].relabel_any;
    }
    double pn = (double)S.inst;
    #define SE(c) sqrt(((double)(c)/pn) * (1.0 - (double)(c)/pn) / pn)
    printf("[cov] instances=%llu seed=0x%llx (real AES-128 7r, true key, delta-set at x_2[0], match byte x_6[0], labels eta=1..16)\n",
           (unsigned long long)S.inst, (unsigned long long)seed);
    printf("[cov] s=0 instances: %llu  rate=%.6f  se=%.6f  (ideal 1/256=%.6f)\n",
           (unsigned long long)S.s_zero, S.s_zero / pn, SE(S.s_zero), 1.0/256);
    u64 td0 = 0, tds = 0;
    for(int j = 0; j < 16; j++){ td0 += S.lab_d0[j]; tds += S.lab_ds[j]; }
    u64 td0_14 = 0, tds_14 = 0;
    for(int j = 0; j < 14; j++){ td0_14 += S.lab_d0[j]; tds_14 += S.lab_ds[j]; }
    printf("[cov] per-label rates over labels 1..16 (N=16n=%llu): d=0: %.6f  d=s(s!=0): %.6f  bad: %.6f  se(bad)=%.6f  (ideal 2/256-1/65536=%.6f ; 2/256=%.6f)\n",
           (unsigned long long)(16*S.inst), td0 / (16.0*pn), tds / (16.0*pn), (td0+tds) / (16.0*pn),
           sqrt(((td0+tds)/(16.0*pn))*(1-(td0+tds)/(16.0*pn))/(16.0*pn)), 2.0/256 - 1.0/65536, 2.0/256);
    printf("[cov] per-label rates over labels 1..14 only: d=0: %.6f  d=s: %.6f  bad: %.6f\n",
           td0_14 / (14.0*pn), tds_14 / (14.0*pn), (td0_14+tds_14) / (14.0*pn));
    for(int j = 0; j < 16; j++)
        printf("[cov] label %2d: d0=%llu ds=%llu rate=%.6f\n", j+1, (unsigned long long)S.lab_d0[j],
               (unsigned long long)S.lab_ds[j], (S.lab_d0[j]+S.lab_ds[j]) / pn);
    double mean14 = 0, mean16 = 0;
    printf("[cov] #bad among labels 1..14 (all instances):");
    for(int k = 0; k < 15; k++){ printf(" %d:%llu", k, (unsigned long long)S.hb14[k]); mean14 += k * (double)S.hb14[k]; }
    printf("  mean=%.6f (14*2/256=%.6f)\n", mean14 / pn, 14.0*(2.0/256-1.0/65536));
    printf("[cov] #bad among labels 1..16 (all instances):");
    for(int k = 0; k < 17; k++){ printf(" %d:%llu", k, (unsigned long long)S.hb16[k]); mean16 += k * (double)S.hb16[k]; }
    printf("  mean=%.6f (16*2/256=%.6f)\n", mean16 / pn, 16.0*(2.0/256-1.0/65536));
    printf("[cov] #bad among 1..14 | s!=0:");
    for(int k = 0; k < 15; k++) printf(" %d:%llu", k, (unsigned long long)S.hb14_s[k]);
    printf("\n[cov] #bad among 1..16 | s!=0:");
    for(int k = 0; k < 17; k++) printf(" %d:%llu", k, (unsigned long long)S.hb16_s[k]);
    printf("\n");
    double eta_plain_an = pow(254.0/256, 14) * (255.0/256);
    double q = 254.0/256, p = 2.0/256;
    double eta_tauC_an = (pow(q,16) + 16*p*pow(q,15) + 15*p*p*pow(q,14)) * (255.0/256);
    double eta_tauA_an = (pow(q,14) + 14*p*pow(q,14)) * (255.0/256);   /* W_j needs 13 clean + label15 clean */
    double eta_ideal_an = (pow(q,16) + 16*p*pow(q,15) + 120*p*p*pow(q,14)) * (255.0/256);
    printf("[cov] PLAIN  eta = %llu/%llu = %.6f  se=%.6f   analytic (254/256)^14*(255/256) = %.6f\n",
           (unsigned long long)S.cov_plain, (unsigned long long)S.inst, S.cov_plain / pn, SE(S.cov_plain), eta_plain_an);
    printf("[cov] TAU-C  eta = %llu/%llu = %.6f  se=%.6f   analytic (q^16+16pq^15+15p^2q^14)*(255/256) = %.6f\n",
           (unsigned long long)S.cov_tauC, (unsigned long long)S.inst, S.cov_tauC / pn, SE(S.cov_tauC), eta_tauC_an);
    printf("[cov] TAU-A  eta = %llu/%llu = %.6f  se=%.6f   analytic (q^14+14pq^14)*(255/256) = %.6f\n",
           (unsigned long long)S.cov_tauA, (unsigned long long)S.inst, S.cov_tauA / pn, SE(S.cov_tauA), eta_tauA_an);
    printf("[cov] IDEAL  eta = %llu/%llu = %.6f  se=%.6f   analytic P(Bin(16,2/256)<=2)*(255/256) = %.6f\n",
           (unsigned long long)S.cov_ideal, (unsigned long long)S.inst, S.cov_ideal / pn, SE(S.cov_ideal), eta_ideal_an);
    printf("[cov] GATE (online Phi == offline Phi on covered instances):\n");
    printf("[cov]   plain: %llu/%llu %s\n", (unsigned long long)S.m_plain, (unsigned long long)S.cov_plain, S.m_plain == S.cov_plain ? "PASS" : "FAIL");
    printf("[cov]   tauC : %llu/%llu %s\n", (unsigned long long)S.m_tauC, (unsigned long long)S.cov_tauC, S.m_tauC == S.cov_tauC ? "PASS" : "FAIL");
    printf("[cov]   tauA : %llu/%llu %s\n", (unsigned long long)S.m_tauA, (unsigned long long)S.cov_tauA, S.m_tauA == S.cov_tauA ? "PASS" : "FAIL");
    printf("[cov]   ideal: %llu/%llu %s\n", (unsigned long long)S.m_ideal, (unsigned long long)S.cov_ideal, S.m_ideal == S.cov_ideal ? "PASS" : "FAIL");
    printf("[cov] tauC probe index used (1=clean window {1..14}, then t=1..14 -> probes 2..15):");
    for(int k = 0; k < 15; k++) printf(" %d:%llu", k+1, (unsigned long long)S.probe_tauC[k]);
    printf("\n[cov] relabel (m!=2) events on plain window: %llu\n", (unsigned long long)S.relabel_any);
    int ok = (S.m_plain == S.cov_plain) && (S.m_tauC == S.cov_tauC) && (S.m_tauA == S.cov_tauA) && (S.m_ideal == S.cov_ideal);
    printf("[cov] %s\n", ok ? "PASS (0 mismatches on covered instances, all rules)" : "FAIL");
    free(ws); free(th);
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* WRONG-CANDIDATE generation                                            */
/* ------------------------------------------------------------------ */
typedef struct __attribute__((packed)) { u8 fp[12]; u32 id; } wrec_t;   /* 16 B */

static u64 G_nwrong = 0, G_seed = 0;
static volatile u64 G_next = 0, G_done = 0;
static gbuck_t G_bset;
static int G_chunklog = 16;
static char G_dir[512];
static volatile int G_stop_mon = 0;
static int G_mon_recsz = 16;

static void *mon_fn(void *arg){
    const char *tag = (const char*)arg;
    double t0 = now_sec();
    while(!G_stop_mon){
        for(int i = 0; i < 40 && !G_stop_mon; i++) usleep(250000);
        double t = now_sec() - t0; u64 d = G_done;
        fprintf(stderr, "[%s mon] t=%.0fs done=%llu (2^%.2f) rate=%.2f M/s  io=%.1f GiB\n", tag, t,
                (unsigned long long)d, log2((double)d + 1), d / t / 1e6, (double)d * G_mon_recsz / (1024.0*1024*1024));
    }
    return NULL;
}

typedef struct { int tid; stats_t st; u64 g_zero, lab_eq_ref; } wrw_t;

static void *wrong_worker(void *arg){
    wrw_t *w = (wrw_t*)arg;
    stats_init(&w->st);
    gtl_t tl; gtl_init(&tl, &G_bset, 4096);
    const u64 chunk = 1ULL << G_chunklog;
    for(;;){
        u64 c = __sync_fetch_and_add(&G_next, 1);
        u64 start = c * chunk;
        if(start >= G_nwrong) break;
        u64 end = start + chunk; if(end > G_nwrong) end = G_nwrong;
        rng_t R; rng_seed2(&R, G_seed, 0x57a7000000ULL + c);
        for(u64 id = start; id < end; id++){
            u8 key[16], P[16], km1[4], k7a[4];
            for(int i = 0; i < 16; i++){ key[i] = rng8(&R); P[i] = rng8(&R); }
            for(int i = 0; i < 4; i++){ km1[i] = rng8(&R); k7a[i] = rng8(&R); }
            u8 k10g = rng8(&R);                    /* wrong guess of k_1[0] -> xi */
            aesctx ctx; aes_set_key(&ctx, key, 7);
            nik_t N; nik_set(&N, &ctx);
            /* build delta-set reference column w1r from P and the (wrong) km1 guess */
            u8 z1r[4], w1r[4];
            for(int i = 0; i < 4; i++) z1r[i] = SBOX[(u8)(P[DIAG0[i]] ^ km1[i])];
            for(int r = 0; r < 4; r++) w1r[r] = (u8)(gmul(MCc[r][0], z1r[0]) ^ gmul(MCc[r][1], z1r[1]) ^ gmul(MCc[r][2], z1r[2]) ^ gmul(MCc[r][3], z1r[3]));
            u8 xi = (u8)(w1r[0] ^ k10g);
            /* reference ciphertext (dv=0): PP = P */
            u8 C0[16]; _mm_storeu_si128((__m128i*)C0, nik_enc(&N, _mm_loadu_si128((const __m128i*)P)));
            u8 x7ref[4];
            for(int i = 0; i < 4; i++) x7ref[i] = iSBOX[(u8)(C0[ADIAG[i]] ^ k7a[i])];
            u8 D[14];
            for(int j = 1; j <= 14; j++){
                u8 v = (u8)(iSBOX[(u8)(SBOX[xi] ^ j)] ^ xi);
                u8 w1[4] = { (u8)(w1r[0] ^ v), w1r[1], w1r[2], w1r[3] };
                u8 z1[4];
                for(int r = 0; r < 4; r++) z1[r] = (u8)(gmul(iMCc[r][0], w1[0]) ^ gmul(iMCc[r][1], w1[1]) ^ gmul(iMCc[r][2], w1[2]) ^ gmul(iMCc[r][3], w1[3]));
                u8 PP[16]; memcpy(PP, P, 16);
                for(int i = 0; i < 4; i++) PP[DIAG0[i]] = (u8)(iSBOX[z1[i]] ^ km1[i]);
                u8 Cv[16];
                _mm_storeu_si128((__m128i*)Cv, nik_enc(&N, _mm_loadu_si128((const __m128i*)PP)));
                u8 dy6 = 0;
                for(int i = 0; i < 4; i++) dy6 ^= gmul(iMCc[0][i], (u8)(iSBOX[(u8)(Cv[ADIAG[i]] ^ k7a[i])] ^ x7ref[i]));
                u8 g = GF_inv[Linv[dy6]];
                if(g == 0) w->g_zero++;
                D[j-1] = g;
            }
            wrec_t rec;
            int m = phi14(D, rec.fp);
            if(m < 0) w->st.degenerate++;
            else if(m != 1) w->st.relabel++;
            rec.id = (u32)id;
            stats_add(&w->st, rec.fp);
            gtl_put(&tl, rec.fp[0], &rec);
        }
        __sync_fetch_and_add(&G_done, end - start);
    }
    gtl_flush(&tl);
    return NULL;
}

static int mode_wrong(int argc, char **argv){
    int logn = (int)argll(argc, argv, "logn", 20);
    int nthr = (int)argll(argc, argv, "threads", 64);
    G_seed = (u64)argll(argc, argv, "seed", 0x57a7e1);
    snprintf(G_dir, sizeof G_dir, "%s", argstr(argc, argv, "dir", "./wrong"));
    G_nwrong = 1ULL << logn;
    if(!nik_selftest()){ fprintf(stderr, "FATAL: AES-NI self-test failed\n"); return 2; }
    fprintf(stderr, "[wrong] AES-NI 7-round == scalar aes_enc on 2000 random (key,P): PASS\n");
    if(gbuck_open(&G_bset, G_dir, "w", sizeof(wrec_t))) return 1;
    fprintf(stderr, "[wrong] N=2^%d wrong-candidate fingerprints, %d threads, seed=%llx -> %s\n",
            logn, nthr, (unsigned long long)G_seed, G_dir);
    wrw_t *ws = (wrw_t*)calloc((size_t)nthr, sizeof(wrw_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    pthread_t mon; double t0 = now_sec();
    pthread_create(&mon, NULL, mon_fn, "wrong");
    for(int i = 0; i < nthr; i++){ ws[i].tid = i; pthread_create(&th[i], NULL, wrong_worker, &ws[i]); }
    for(int i = 0; i < nthr; i++) pthread_join(th[i], NULL);
    G_stop_mon = 1; pthread_join(mon, NULL);
    gbuck_close(&G_bset);
    double t1 = now_sec();
    /* merge stats */
    stats_t S; stats_init(&S); u64 gz = 0;
    for(int i = 0; i < nthr; i++){ stats_merge(&S, &ws[i].st); gz += ws[i].g_zero; }
    char path[700]; snprintf(path, sizeof path, "%s/stats.txt", G_dir);
    stats_write(&S, path);
    u64 totc = 0; for(int b = 0; b < 256; b++) totc += G_bset.count[b];
    fprintf(stderr, "[wrong] done: %llu records in %.1fs (%.2f M/s). relabel(m!=2)=%llu degenerate=%llu g=0 labels=%llu\n",
            (unsigned long long)totc, t1 - t0, totc / (t1 - t0) / 1e6,
            (unsigned long long)S.relabel, (unsigned long long)S.degenerate, (unsigned long long)gz);
    printf("[wrong] records=%llu relabel=%llu degenerate=%llu g_zero_labels=%llu (of %llu labels)\n",
           (unsigned long long)totc, (unsigned long long)S.relabel, (unsigned long long)S.degenerate,
           (unsigned long long)gz, (unsigned long long)(14 * totc));
    return 0;
}

/* ------------------------------------------------------------------ */
/* IDEAL baseline: Phi of 14 iid uniform bytes (same record format)     */
/* ------------------------------------------------------------------ */
static void *ideal_worker(void *arg){
    wrw_t *w = (wrw_t*)arg;
    stats_init(&w->st);
    gtl_t tl; gtl_init(&tl, &G_bset, 4096);
    const u64 chunk = 1ULL << G_chunklog;
    for(;;){
        u64 c = __sync_fetch_and_add(&G_next, 1);
        u64 start = c * chunk;
        if(start >= G_nwrong) break;
        u64 end = start + chunk; if(end > G_nwrong) end = G_nwrong;
        rng_t R; rng_seed2(&R, G_seed, 0x1dea1000000ULL + c);
        for(u64 id = start; id < end; id++){
            u8 D[14];
            u64 r = rng_next(&R); memcpy(D, &r, 8);
            r = rng_next(&R); memcpy(D + 8, &r, 6);
            wrec_t rec;
            int m = phi14(D, rec.fp);
            if(m < 0) w->st.degenerate++; else if(m != 1) w->st.relabel++;
            rec.id = (u32)id;
            stats_add(&w->st, rec.fp);
            gtl_put(&tl, rec.fp[0], &rec);
        }
        __sync_fetch_and_add(&G_done, end - start);
    }
    gtl_flush(&tl);
    return NULL;
}

static int mode_ideal(int argc, char **argv){
    int logn = (int)argll(argc, argv, "logn", 20);
    int nthr = (int)argll(argc, argv, "threads", 64);
    G_seed = (u64)argll(argc, argv, "seed", 0x1dea1);
    snprintf(G_dir, sizeof G_dir, "%s", argstr(argc, argv, "dir", "./ideal"));
    G_nwrong = 1ULL << logn;
    if(gbuck_open(&G_bset, G_dir, "w", sizeof(wrec_t))) return 1;
    fprintf(stderr, "[ideal] N=2^%d Phi(14 iid uniform bytes), %d threads, seed=%llx -> %s\n",
            logn, nthr, (unsigned long long)G_seed, G_dir);
    wrw_t *ws = (wrw_t*)calloc((size_t)nthr, sizeof(wrw_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    pthread_t mon; double t0 = now_sec();
    pthread_create(&mon, NULL, mon_fn, "ideal");
    for(int i = 0; i < nthr; i++){ ws[i].tid = i; pthread_create(&th[i], NULL, ideal_worker, &ws[i]); }
    for(int i = 0; i < nthr; i++) pthread_join(th[i], NULL);
    G_stop_mon = 1; pthread_join(mon, NULL);
    gbuck_close(&G_bset);
    double t1 = now_sec();
    stats_t S; stats_init(&S);
    for(int i = 0; i < nthr; i++) stats_merge(&S, &ws[i].st);
    char path[700]; snprintf(path, sizeof path, "%s/stats.txt", G_dir);
    stats_write(&S, path);
    u64 totc = 0; for(int b = 0; b < 256; b++) totc += G_bset.count[b];
    fprintf(stderr, "[ideal] done: %llu records in %.1fs (%.2f M/s). relabel(m!=2)=%llu degenerate=%llu\n",
            (unsigned long long)totc, t1 - t0, totc / (t1 - t0) / 1e6,
            (unsigned long long)S.relabel, (unsigned long long)S.degenerate);
    printf("[ideal] records=%llu relabel=%llu degenerate=%llu\n",
           (unsigned long long)totc, (unsigned long long)S.relabel, (unsigned long long)S.degenerate);
    return 0;
}

/* ------------------------------------------------------------------ */
/* TABLE-ENTRY generation                                                */
/* ------------------------------------------------------------------ */
typedef struct __attribute__((packed)) {
    u8 fp[12]; u8 d[14]; u16 flags; u64 seqh; u32 base; u32 idx;
} trec_t;                                                   /* 44 B */

static int G_logbases = 10, G_walkbits = 12;

static void *table_worker(void *arg){
    wrw_t *w = (wrw_t*)arg;
    stats_init(&w->st);
    gtl_t tl; gtl_init(&tl, &G_bset, 2048);
    walk_t *W = (walk_t*)aligned_alloc(64, sizeof(walk_t));
    const u64 nb = 1ULL << G_logbases;
    const u32 NW = 1u << G_walkbits;
    for(;;){
        u64 b = __sync_fetch_and_add(&G_next, 1);
        if(b >= nb) break;
        rng_t R; rng_seed2(&R, G_seed, 0xba5e00000ULL + b);
        base_t B; base_sample(&B, &R);
        walk_init(W, &B);
        for(u32 i = 0; i < NW; i++){
            if(i) walk_step(W, __builtin_ctz(i));
            trec_t rec;
            u8 D[14];
            for(int j = 1; j <= 14; j++){ rec.d[j-1] = W->E[j]; D[j-1] = GF_inv[W->E[j]]; }
            int m = phi14(D, rec.fp);
            rec.flags = (u16)((m < 0) ? 0xffff : m);
            if(m < 0) w->st.degenerate++; else if(m != 1) w->st.relabel++;
            int zw = 0; for(int j = 0; j < 14; j++) if(rec.d[j] == 0) zw++;
            if(zw) w->st.zero_in_window++;
            rec.seqh = seq_hash(W->E);
            rec.base = (u32)b; rec.idx = i;
            stats_add(&w->st, rec.fp);
            gtl_put(&tl, rec.fp[0], &rec);
        }
        __sync_fetch_and_add(&G_done, NW);
    }
    gtl_flush(&tl);
    free(W);
    return NULL;
}

static int mode_table(int argc, char **argv){
    G_logbases = (int)argll(argc, argv, "logbases", 10);
    G_walkbits = (int)argll(argc, argv, "walkbits", 12);
    int nthr = (int)argll(argc, argv, "threads", 64);
    G_seed = (u64)argll(argc, argv, "seed", 0x7ab1e5);
    snprintf(G_dir, sizeof G_dir, "%s", argstr(argc, argv, "dir", "./table"));
    G_mon_recsz = (int)sizeof(trec_t);
    if(gbuck_open(&G_bset, G_dir, "t", sizeof(trec_t))) return 1;
    /* self-test: walk vs cold on a few bases */
    {
        rng_t R; rng_seed2(&R, G_seed ^ 0x5151, 7);
        walk_t *W = (walk_t*)aligned_alloc(64, sizeof(walk_t));
        int ok = 0, tot = 0;
        for(int t = 0; t < 6; t++){
            base_t B; base_sample(&B, &R);
            if(!base_check(&B)){ fprintf(stderr, "FATAL base_check\n"); return 2; }
            walk_init(W, &B);
            for(int i = 0; i < 2000; i++){
                if(i) walk_step(W, __builtin_ctz(i));
                if(i % 131) continue;
                u8 Ec[256], Ep[256]; cold_E(&W->R, Ec);
                refs24_t P; ref24_to_refs24(&W->R, &P); prop2_all(&P, Ep);
                tot++;
                if(!memcmp(Ec, W->E, 256) && !memcmp(Ep, W->E, 256)) ok++;
            }
        }
        free(W);
        fprintf(stderr, "[table] selftest walk==cold_E==prop2_all %d/%d\n", ok, tot);
        if(ok != tot){ fprintf(stderr, "FATAL selftest\n"); return 3; }
    }
    fprintf(stderr, "[table] N=2^%d entries = 2^%d bases x 2^%d Gray-walk steps, %d threads, seed=%llx -> %s\n",
            G_logbases + G_walkbits, G_logbases, G_walkbits, nthr, (unsigned long long)G_seed, G_dir);
    wrw_t *ws = (wrw_t*)calloc((size_t)nthr, sizeof(wrw_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    pthread_t mon; double t0 = now_sec();
    pthread_create(&mon, NULL, mon_fn, "table");
    for(int i = 0; i < nthr; i++){ ws[i].tid = i; pthread_create(&th[i], NULL, table_worker, &ws[i]); }
    for(int i = 0; i < nthr; i++) pthread_join(th[i], NULL);
    G_stop_mon = 1; pthread_join(mon, NULL);
    gbuck_close(&G_bset);
    double t1 = now_sec();
    stats_t S; stats_init(&S);
    for(int i = 0; i < nthr; i++) stats_merge(&S, &ws[i].st);
    char path[700]; snprintf(path, sizeof path, "%s/stats.txt", G_dir);
    stats_write(&S, path);
    u64 totc = 0; for(int b = 0; b < 256; b++) totc += G_bset.count[b];
    fprintf(stderr, "[table] done: %llu records in %.1fs (%.2f M/s). relabel=%llu degenerate=%llu zero_in_window=%llu\n",
            (unsigned long long)totc, t1 - t0, totc / (t1 - t0) / 1e6,
            (unsigned long long)S.relabel, (unsigned long long)S.degenerate, (unsigned long long)S.zero_in_window);
    printf("[table] records=%llu relabel=%llu degenerate=%llu zero_in_window=%llu\n",
           (unsigned long long)totc, (unsigned long long)S.relabel, (unsigned long long)S.degenerate,
           (unsigned long long)S.zero_in_window);
    return 0;
}

/* ------------------------------------------------------------------ */
/* ANALYSIS: per-bucket sort, prefix-collision counts, sorted key export  */
/* ------------------------------------------------------------------ */
/* counts of pairs sharing a k-byte prefix, k = 4..12 */
typedef struct { u64 pairs[13]; u64 dupgroups, duprecs; } coll_t;

static int cmp12(const void *a, const void *b){ return memcmp(a, b, 12); }
static int cmp_trec(const void *a, const void *b){ return memcmp(a, b, 12); }

static u8 *read_file(const char *path, size_t *len){
    FILE *f = fopen(path, "rb"); if(!f){ perror(path); return NULL; }
    fseeko(f, 0, SEEK_END); off_t sz = ftello(f); fseeko(f, 0, SEEK_SET);
    u8 *buf = (u8*)malloc((size_t)sz ? (size_t)sz : 1);
    size_t got = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    *len = got;
    return buf;
}

/* after sorting records by the 12-byte key (stride given), accumulate
 * prefix-collision pair counts for k=4..12 into C */
static void scan_prefix_coll(const u8 *base, u64 n, size_t stride, coll_t *C){
    /* for each k, runs of equal k-prefix: count pairs in one pass using the
     * common-prefix length between consecutive sorted records. */
    if(!n) return;
    u64 run[13];
    for(int k = 0; k < 13; k++) run[k] = 1;
    for(u64 i = 1; i < n; i++){
        const u8 *a = base + (i - 1) * stride, *b = base + i * stride;
        int cpl = 0;
        while(cpl < 12 && a[cpl] == b[cpl]) cpl++;
        for(int k = cpl + 1; k <= 12; k++){ C->pairs[k] += run[k] * (run[k] - 1) / 2; run[k] = 1; }
        for(int k = 1; k <= cpl; k++) run[k]++;
    }
    for(int k = 1; k <= 12; k++) C->pairs[k] += run[k] * (run[k] - 1) / 2;
}

typedef struct {
    int tid, nthr; const char *dir; const char *pfx; int recsz; int write_keys;
    coll_t C; u64 total;
    /* table-side duplicate classification */
    u64 dg_total, dg_samebase, dg_crossbase, dg_window_ident, dg_window_affine, dg_seq_ident, dg_other;
    FILE *ex; int nex;
} anw_t;

static pthread_mutex_t g_ex_mtx = PTHREAD_MUTEX_INITIALIZER;
static volatile int g_next_bucket = 0;

static void *ana_worker(void *arg){
    anw_t *w = (anw_t*)arg;
    char path[700];
    for(;;){
        int b = __sync_fetch_and_add(&g_next_bucket, 1);
        if(b >= 256) break;
        snprintf(path, sizeof path, "%s/%s%02x.bin", w->dir, w->pfx, b);
        size_t len = 0; u8 *buf = read_file(path, &len);
        if(!buf) continue;
        u64 n = len / (u64)w->recsz;
        w->total += n;
        qsort(buf, (size_t)n, (size_t)w->recsz, w->recsz == 16 ? cmp12 : cmp_trec);
        scan_prefix_coll(buf, n, (size_t)w->recsz, &w->C);
        /* full-duplicate groups */
        u64 i = 0;
        while(i < n){
            u64 j = i + 1;
            while(j < n && !memcmp(buf + i * w->recsz, buf + j * w->recsz, 12)) j++;
            if(j - i >= 2){
                w->C.dupgroups++; w->C.duprecs += j - i;
                if(w->recsz == (int)sizeof(trec_t)){
                    /* classify every pair against the first member */
                    const trec_t *r0 = (const trec_t*)(buf + i * w->recsz);
                    for(u64 q = i + 1; q < j; q++){
                        const trec_t *r1 = (const trec_t*)(buf + q * w->recsz);
                        w->dg_total++;
                        if(r0->base == r1->base) w->dg_samebase++; else w->dg_crossbase++;
                        if(!memcmp(r0->d, r1->d, 14)){
                            if(r0->seqh == r1->seqh) w->dg_seq_ident++; else w->dg_window_ident++;
                        } else {
                            /* verify affine equivalence of the two D-windows */
                            u8 D0[14], D1[14];
                            for(int t = 0; t < 14; t++){ D0[t] = GF_inv[r0->d[t]]; D1[t] = GF_inv[r1->d[t]]; }
                            int m = (r0->flags == 0xffff) ? -1 : r0->flags;
                            int affine_ok = 0; u8 alpha = 0, beta = 0;
                            if(m >= 0 && (u8)(D0[m] ^ D0[0])){
                                alpha = gmul((u8)(D1[m] ^ D1[0]), GF_inv[(u8)(D0[m] ^ D0[0])]);
                                beta = (u8)(D1[0] ^ gmul(alpha, D0[0]));
                                affine_ok = (alpha != 0);
                                for(int t = 0; t < 14 && affine_ok; t++) if(D1[t] != (u8)(gmul(alpha, D0[t]) ^ beta)) affine_ok = 0;
                            }
                            if(affine_ok) w->dg_window_affine++; else w->dg_other++;
                            pthread_mutex_lock(&g_ex_mtx);
                            if(w->nex < 40 && w->ex){
                                w->nex++;
                                fprintf(w->ex, "DUP fp=");
                                for(int t = 0; t < 12; t++) fprintf(w->ex, "%02x", r0->fp[t]);
                                fprintf(w->ex, " base=%u idx=%u d=", r0->base, r0->idx);
                                for(int t = 0; t < 14; t++) fprintf(w->ex, "%02x", r0->d[t]);
                                fprintf(w->ex, " seqh=%016llx | base=%u idx=%u d=", (unsigned long long)r0->seqh, r1->base, r1->idx);
                                for(int t = 0; t < 14; t++) fprintf(w->ex, "%02x", r1->d[t]);
                                fprintf(w->ex, " seqh=%016llx affine=%d alpha=%02x beta=%02x m=%d\n",
                                        (unsigned long long)r1->seqh, affine_ok, alpha, beta, m);
                                fflush(w->ex);
                            }
                            pthread_mutex_unlock(&g_ex_mtx);
                        }
                    }
                }
            }
            i = j;
        }
        /* export sorted 12-byte keys for cross-matching */
        if(w->write_keys){
            snprintf(path, sizeof path, "%s/k%02x.key", w->dir, b);
            FILE *f = fopen(path, "wb");
            if(f){
                if(w->recsz == 12){ fwrite(buf, 12, (size_t)n, f); }
                else {
                    u8 *k = (u8*)malloc((size_t)n * 12);
                    for(u64 q = 0; q < n; q++) memcpy(k + q * 12, buf + q * w->recsz, 12);
                    fwrite(k, 12, (size_t)n, f); free(k);
                }
                fclose(f);
            }
        }
        free(buf);
        fprintf(stderr, "[ana] bucket %02x n=%llu done\n", b, (unsigned long long)n);
    }
    return NULL;
}

static int mode_analyze(int argc, char **argv, int table){
    const char *dir = argstr(argc, argv, "dir", table ? "./table" : "./wrong");
    int nthr = (int)argll(argc, argv, "threads", 32);
    int write_keys = (int)argll(argc, argv, "keys", 1);
    anw_t *ws = (anw_t*)calloc((size_t)nthr, sizeof(anw_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    FILE *ex = NULL;
    if(table){
        char path[700]; snprintf(path, sizeof path, "%s/dup_examples.txt", dir);
        ex = fopen(path, "w");
    }
    g_next_bucket = 0;
    double t0 = now_sec();
    for(int i = 0; i < nthr; i++){
        ws[i].tid = i; ws[i].nthr = nthr; ws[i].dir = dir; ws[i].pfx = table ? "t" : "w";
        ws[i].recsz = table ? (int)sizeof(trec_t) : (int)sizeof(wrec_t); ws[i].write_keys = write_keys;
        ws[i].ex = ex;
        pthread_create(&th[i], NULL, ana_worker, &ws[i]);
    }
    coll_t C; memset(&C, 0, sizeof C);
    u64 N = 0, dg_total = 0, dg_same = 0, dg_cross = 0, dg_wi = 0, dg_wa = 0, dg_seq = 0, dg_other = 0;
    for(int i = 0; i < nthr; i++){
        pthread_join(th[i], NULL);
        for(int k = 0; k < 13; k++) C.pairs[k] += ws[i].C.pairs[k];
        C.dupgroups += ws[i].C.dupgroups; C.duprecs += ws[i].C.duprecs; N += ws[i].total;
        dg_total += ws[i].dg_total; dg_same += ws[i].dg_samebase; dg_cross += ws[i].dg_crossbase;
        dg_wi += ws[i].dg_window_ident; dg_wa += ws[i].dg_window_affine; dg_seq += ws[i].dg_seq_ident; dg_other += ws[i].dg_other;
    }
    if(ex) fclose(ex);
    double t1 = now_sec();
    double pairs = (double)N * ((double)N - 1) / 2.0;
    printf("[%s] N=%llu records (2^%.3f), total pairs=2^%.3f, analysis %.1fs\n", table ? "anat" : "anaw",
           (unsigned long long)N, log2((double)N), log2(pairs), t1 - t0);
    printf("[%s] prefix-collision pair counts (observed vs ideal C(N,2)*2^-8k, and -log2 empirical collision prob):\n",
           table ? "anat" : "anaw");
    for(int k = 4; k <= 12; k++){
        double ideal = pairs * pow(2.0, -8.0 * k);
        double phat = (double)C.pairs[k] / pairs;
        printf("  k=%2d bytes (%3d bits): observed=%llu  ideal=%.4g  ratio=%.4f  -log2(p_coll)=%s%.2f\n",
               k, 8 * k, (unsigned long long)C.pairs[k], ideal,
               ideal > 0 ? (double)C.pairs[k] / ideal : 0.0,
               C.pairs[k] ? "" : ">", C.pairs[k] ? -log2(phat) : -log2(3.0 / pairs));
    }
    printf("[%s] full 96-bit duplicate groups=%llu records-in-dup-groups=%llu\n", table ? "anat" : "anaw",
           (unsigned long long)C.dupgroups, (unsigned long long)C.duprecs);
    if(table){
        printf("[anat] structural duplicate pairs classified (vs first group member): total=%llu samebase=%llu crossbase=%llu; "
               "identical-seq=%llu identical-window(diff seq)=%llu affine-equivalent-window=%llu other=%llu\n",
               (unsigned long long)dg_total, (unsigned long long)dg_same, (unsigned long long)dg_cross,
               (unsigned long long)dg_seq, (unsigned long long)dg_wi, (unsigned long long)dg_wa, (unsigned long long)dg_other);
        printf("[anat] examples in %s/dup_examples.txt\n", dir);
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* CROSS-MATCH: sorted wrong keys vs sorted table keys                    */
/* ------------------------------------------------------------------ */
typedef struct {
    const char *wdir, *tdir;
    u64 probes, tentries;
    u64 m32, m64, m96;
} crw_t;

/* count pairs (i in A, j in B) with equal pfx bytes; A,B sorted, 12-byte keys */
static u64 count_prefix_pairs(const u8 *A, u64 na, const u8 *B, u64 nb, int pfx){
    u64 i = 0, j = 0, total = 0;
    while(i < na && j < nb){
        int c = memcmp(A + i * 12, B + j * 12, (size_t)pfx);
        if(c < 0){ i++; continue; }
        if(c > 0){ j++; continue; }
        u64 i2 = i + 1, j2 = j + 1;
        while(i2 < na && !memcmp(A + i2 * 12, A + i * 12, (size_t)pfx)) i2++;
        while(j2 < nb && !memcmp(B + j2 * 12, B + j * 12, (size_t)pfx)) j2++;
        total += (i2 - i) * (j2 - j);
        i = i2; j = j2;
    }
    return total;
}

static void *cross_worker(void *arg){
    crw_t *w = (crw_t*)arg;
    char path[700];
    for(;;){
        int b = __sync_fetch_and_add(&g_next_bucket, 1);
        if(b >= 256) break;
        size_t la = 0, lb = 0;
        snprintf(path, sizeof path, "%s/k%02x.key", w->wdir, b);
        u8 *A = read_file(path, &la);
        snprintf(path, sizeof path, "%s/k%02x.key", w->tdir, b);
        u8 *B = read_file(path, &lb);
        if(!A || !B){ free(A); free(B); continue; }
        u64 na = la / 12, nb = lb / 12;
        w->probes += na; w->tentries += nb;
        w->m32 += count_prefix_pairs(A, na, B, nb, 4);
        w->m64 += count_prefix_pairs(A, na, B, nb, 8);
        w->m96 += count_prefix_pairs(A, na, B, nb, 12);
        free(A); free(B);
        fprintf(stderr, "[cross] bucket %02x done (na=%llu nb=%llu)\n", b, (unsigned long long)na, (unsigned long long)nb);
    }
    return NULL;
}

static int mode_cross(int argc, char **argv){
    const char *wdir = argstr(argc, argv, "wdir", "./wrong");
    const char *tdir = argstr(argc, argv, "tdir", "./table");
    int nthr = (int)argll(argc, argv, "threads", 32);
    crw_t *ws = (crw_t*)calloc((size_t)nthr, sizeof(crw_t));
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * (size_t)nthr);
    g_next_bucket = 0;
    double t0 = now_sec();
    for(int i = 0; i < nthr; i++){ ws[i].wdir = wdir; ws[i].tdir = tdir; pthread_create(&th[i], NULL, cross_worker, &ws[i]); }
    u64 P = 0, T = 0, m32 = 0, m64 = 0, m96 = 0;
    for(int i = 0; i < nthr; i++){
        pthread_join(th[i], NULL);
        P += ws[i].probes; T += ws[i].tentries; m32 += ws[i].m32; m64 += ws[i].m64; m96 += ws[i].m96;
    }
    double t1 = now_sec();
    double pairs = (double)P * (double)T;
    printf("[cross] probes(wrong)=%llu (2^%.2f) table entries=%llu (2^%.2f) probe-pairs=2^%.2f  time %.1fs\n",
           (unsigned long long)P, log2((double)P), (unsigned long long)T, log2((double)T), log2(pairs), t1 - t0);
    printf("[cross] 32-bit prefix matches: observed=%llu ideal=%.4g ratio=%.4f\n",
           (unsigned long long)m32, pairs * pow(2, -32), m32 / (pairs * pow(2, -32)));
    printf("[cross] 64-bit prefix matches: observed=%llu ideal=%.4g\n",
           (unsigned long long)m64, pairs * pow(2, -64));
    printf("[cross] 96-bit full matches:   observed=%llu ideal=%.4g\n",
           (unsigned long long)m96, pairs * pow(2, -96));
    return 0;
}

/* ------------------------------------------------------------------ */
int main(int argc, char **argv){
    const char *mode = argc > 1 ? argv[1] : "help";
    aesbench_init();
    if(!strcmp(mode, "gate1")) return mode_gate1(argc, argv);
    if(!strcmp(mode, "gate2")) return mode_gate2(argc, argv);
    if(!strcmp(mode, "cov"))   return mode_cov(argc, argv);
    if(!strcmp(mode, "wrong")) return mode_wrong(argc, argv);
    if(!strcmp(mode, "ideal")) return mode_ideal(argc, argv);
    if(!strcmp(mode, "table")) return mode_table(argc, argv);
    if(!strcmp(mode, "anaw"))  return mode_analyze(argc, argv, 0);
    if(!strcmp(mode, "anat"))  return mode_analyze(argc, argv, 1);
    if(!strcmp(mode, "cross")) return mode_cross(argc, argv);
    fprintf(stderr, "modes: gate1 gate2 wrong table anaw anat cross\n");
    return 1;
}
