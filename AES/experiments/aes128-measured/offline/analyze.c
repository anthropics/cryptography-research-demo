// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* analyze.c -- sort/dedup/statistics over the staged bucket files, the
 * wrong-candidate probe test, and the known-answer coverage test.
 *
 * Modes (run `analyze all ...` for the whole pipeline in one process):
 *   sort      : sort each bucket by (fp0,fp1) in RAM, rewrite it sorted,
 *               per-bucket stats, duplicate pairs, fp0-word collisions,
 *               top-20-bit histogram.
 *   fp1stats  : scan sorted buckets, partition fp1 values by top byte in
 *               RAM, sort, count fp1-word collisions.
 *   probe     : generate NP genuine online wrong-candidate fingerprints
 *               (real 7-round AES delta-sets, honest peel with wrong k7a),
 *               join against the table (both words, all prefix levels).
 *   ka        : known-answer presence + genuine online recovery.
 *   all       : sort, fp1stats, probe, ka in sequence.
 *
 * Build: gcc -O3 -march=native -pthread -o analyze analyze.c -lm
 */
#include "common.h"
#include <sys/stat.h>

static int    g_threads = 190;
static char   g_dir[512] = "./buckets";
static int    g_sortpar = 24;          /* buckets sorted concurrently */
static u64    g_probe_inst = 1u << 18, g_probe_g = 1u << 16;   /* 2^34 probes */
static u64    g_seed = 0xA5E5BE5C0FFEEULL;
static int    g_nka = 128, g_nb = 65536;
static double g_tsc_hz;
static FILE  *g_log;                    /* analysis log (findings) */

#define LOGF(...) do { fprintf(stdout, __VA_ARGS__); if(g_log){ fprintf(g_log, __VA_ARGS__); fflush(g_log); } fflush(stdout); } while(0)

/* ------------------------------------------------------------------ */
/* sorting primitives                                                    */
/* ------------------------------------------------------------------ */
static inline int rec_lt(const rec_t *a, const rec_t *b){
    return a->fp0 < b->fp0 || (a->fp0 == b->fp0 && a->fp1 < b->fp1);
}
static void isort16(rec_t *a, long n){
    for(long i = 1; i < n; i++){ rec_t k = a[i]; long j = i - 1; while(j >= 0 && rec_lt(&k, &a[j])){ a[j+1] = a[j]; j--; } a[j+1] = k; }
}
static void qsort16(rec_t *a, long n){
    /* iterative quicksort, median of three, insertion for small */
    typedef struct { long lo, hi; } range_t;
    range_t stack[128]; int sp = 0;
    stack[sp].lo = 0; stack[sp].hi = n - 1; sp++;
    while(sp){
        long lo = stack[--sp].lo, hi = stack[sp].hi;
        while(hi - lo > 24){
            long mid = lo + (hi - lo) / 2;
            if(rec_lt(&a[mid], &a[lo])) { rec_t t = a[mid]; a[mid] = a[lo]; a[lo] = t; }
            if(rec_lt(&a[hi], &a[lo]))  { rec_t t = a[hi]; a[hi] = a[lo]; a[lo] = t; }
            if(rec_lt(&a[hi], &a[mid])) { rec_t t = a[hi]; a[hi] = a[mid]; a[mid] = t; }
            rec_t p = a[mid];
            long i = lo, j = hi;
            for(;;){
                while(rec_lt(&a[i], &p)) i++;
                while(rec_lt(&p, &a[j])) j--;
                if(i >= j) break;
                rec_t t = a[i]; a[i] = a[j]; a[j] = t; i++; j--;
            }
            /* push larger side, loop on smaller */
            if(j - lo < hi - j){ if(j + 1 < hi){ stack[sp].lo = j+1; stack[sp].hi = hi; sp++; } hi = j; }
            else { if(lo < j){ stack[sp].lo = lo; stack[sp].hi = j; sp++; } lo = j + 1; }
            if(sp > 126){ isort16(a + lo, hi - lo + 1); break; }
        }
        if(hi - lo <= 24 && hi > lo) isort16(a + lo, hi - lo + 1);
    }
}
/* full sort of a bucket: split by byte 6 of fp0, then sort each segment */
static void sort_bucket_records(rec_t *a, long n, rec_t *scratch){
    long cnt[257]; memset(cnt, 0, sizeof cnt);
    for(long i = 0; i < n; i++) cnt[((a[i].fp0 >> 48) & 0xff) + 1]++;
    for(int i = 0; i < 256; i++) cnt[i+1] += cnt[i];
    long pos[256]; for(int i = 0; i < 256; i++) pos[i] = cnt[i];
    for(long i = 0; i < n; i++) scratch[pos[(a[i].fp0 >> 48) & 0xff]++] = a[i];
    for(int s = 0; s < 256; s++){ long lo = cnt[s], hi = cnt[s+1]; if(hi - lo > 1) qsort16(scratch + lo, hi - lo); }
    memcpy(a, scratch, (size_t)n * sizeof(rec_t));
}
/* u64 quicksort */
static void isort8(u64 *a, long n){
    for(long i = 1; i < n; i++){ u64 k = a[i]; long j = i - 1; while(j >= 0 && a[j] > k){ a[j+1] = a[j]; j--; } a[j+1] = k; }
}
static void qsort8(u64 *a, long n){
    typedef struct { long lo, hi; } range_t;
    range_t stack[128]; int sp = 0;
    stack[sp].lo = 0; stack[sp].hi = n - 1; sp++;
    while(sp){
        long lo = stack[--sp].lo, hi = stack[sp].hi;
        while(hi - lo > 24){
            long mid = lo + (hi - lo) / 2;
            if(a[mid] < a[lo]) { u64 t = a[mid]; a[mid] = a[lo]; a[lo] = t; }
            if(a[hi] < a[lo])  { u64 t = a[hi]; a[hi] = a[lo]; a[lo] = t; }
            if(a[hi] < a[mid]) { u64 t = a[hi]; a[hi] = a[mid]; a[mid] = t; }
            u64 p = a[mid];
            long i = lo, j = hi;
            for(;;){
                while(a[i] < p) i++;
                while(p < a[j]) j--;
                if(i >= j) break;
                u64 t = a[i]; a[i] = a[j]; a[j] = t; i++; j--;
            }
            if(j - lo < hi - j){ if(j + 1 < hi){ stack[sp].lo = j+1; stack[sp].hi = hi; sp++; } hi = j; }
            else { if(lo < j){ stack[sp].lo = lo; stack[sp].hi = j; sp++; } lo = j + 1; }
            if(sp > 126){ isort8(a + lo, hi - lo + 1); break; }
        }
        if(hi - lo <= 24 && hi > lo) isort8(a + lo, hi - lo + 1);
    }
}
/* sort u64 array of values sharing the same top byte: split on byte 6 */
static void sort_u64_part(u64 *a, long n, u64 *scratch){
    long cnt[257]; memset(cnt, 0, sizeof cnt);
    for(long i = 0; i < n; i++) cnt[((a[i] >> 48) & 0xff) + 1]++;
    for(int i = 0; i < 256; i++) cnt[i+1] += cnt[i];
    long pos[256]; for(int i = 0; i < 256; i++) pos[i] = cnt[i];
    for(long i = 0; i < n; i++) scratch[pos[(a[i] >> 48) & 0xff]++] = a[i];
    for(int s = 0; s < 256; s++){ long lo = cnt[s], hi = cnt[s+1]; if(hi - lo > 1) qsort8(scratch + lo, hi - lo); }
    memcpy(a, scratch, (size_t)n * sizeof(u64));
}

/* ------------------------------------------------------------------ */
/* file helpers                                                          */
/* ------------------------------------------------------------------ */
static void bucket_path(char *p, size_t cap, int b){ snprintf(p, cap, "%s/b%02x.bin", g_dir, b); }
static off_t file_size(const char *p){ struct stat st; if(stat(p, &st)) return -1; return st.st_size; }
static int read_all(const char *p, void *buf, size_t len){
    int fd = open(p, O_RDONLY); if(fd < 0){ perror(p); return -1; }
    size_t off = 0;
    while(off < len){ ssize_t r = pread(fd, (char*)buf + off, len - off, (off_t)off); if(r <= 0){ perror("pread"); close(fd); return -1; } off += (size_t)r; }
    close(fd); return 0;
}
static int write_all(const char *p, const void *buf, size_t len){
    int fd = open(p, O_WRONLY | O_CREAT | O_TRUNC, 0644); if(fd < 0){ perror(p); return -1; }
    size_t off = 0;
    while(off < len){ ssize_t w = pwrite(fd, (const char*)buf + off, len - off, (off_t)off); if(w <= 0){ perror("pwrite"); close(fd); return -1; } off += (size_t)w; }
    fsync(fd); close(fd); return 0;
}
/* on-disk binary search for fp0 in a sorted bucket file: returns index of
 * first record with rec.fp0 >= v (and reads a window of records around). */
static long disk_lower_bound_fp0(int fd, long nrec, u64 v){
    long lo = 0, hi = nrec;
    rec_t r;
    while(lo < hi){
        long mid = lo + (hi - lo) / 2;
        if(pread(fd, &r, sizeof r, (off_t)mid * (off_t)sizeof r) != (ssize_t)sizeof r){ perror("pread"); return -1; }
        if(r.fp0 < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* ------------------------------------------------------------------ */
/* PHASE SORT                                                            */
/* ------------------------------------------------------------------ */
typedef struct {
    u64 n;            /* records */
    u64 dup_pairs;    /* records equal to predecessor in both words */
    u64 distinct;     /* distinct (fp0,fp1) pairs */
    u64 coll0;        /* pairs of records sharing fp0 but differing fp1 (counted per excess record) */
    u64 fp0_distinct;
} bstat_t;
static bstat_t g_bstat[NBUCKET];
static u64 *g_hist20 = NULL;          /* 2^20 counters: top 20 bits of fp0 */
static u64 g_fp1top[NBUCKET];         /* histogram of fp1 top byte (filled in sort phase) */
static int g_have_fp1top = 0;
static pthread_mutex_t g_list_mtx = PTHREAD_MUTEX_INITIALIZER;
static FILE *g_fdups = NULL, *g_fcoll0 = NULL;
static volatile int g_next_bucket = 0;

static void *sort_worker(void *arg){
    (void)arg;
    for(;;){
        int b = __sync_fetch_and_add(&g_next_bucket, 1);
        if(b >= NBUCKET) break;
        char p[700], ps[720]; bucket_path(p, sizeof p, b);
        snprintf(ps, sizeof ps, "%s.sorting", p);
        off_t sz = file_size(p);
        if(sz < 0){ fprintf(stderr, "missing %s\n", p); continue; }
        long n = (long)(sz / (off_t)sizeof(rec_t));
        rec_t *a = (rec_t*)malloc((size_t)sz + 16);
        rec_t *scr = (rec_t*)malloc((size_t)sz + 16);
        if(!a || !scr){ fprintf(stderr, "malloc fail bucket %d\n", b); exit(9); }
        double t0 = now_sec();
        if(read_all(p, a, (size_t)sz)){ free(a); free(scr); continue; }
        double t1 = now_sec();
        sort_bucket_records(a, n, scr);
        double t2 = now_sec();
        /* stats pass */
        bstat_t S; memset(&S, 0, sizeof S);
        S.n = (u64)n;
        u64 *h = g_hist20 + ((u64)b << 12);
        u64 t1h[NBUCKET]; memset(t1h, 0, sizeof t1h);
        for(long i = 0; i < n; i++){
            h[(a[i].fp0 >> 44) & 0xfff]++;
            t1h[a[i].fp1 >> 56]++;
            if(i == 0){ S.distinct = 1; S.fp0_distinct = 1; continue; }
            int same0 = (a[i].fp0 == a[i-1].fp0);
            int same1 = (a[i].fp1 == a[i-1].fp1);
            if(same0 && same1){
                S.dup_pairs++;
                /* log the duplicated value once per run */
                if(i < 2 || !(a[i-2].fp0 == a[i].fp0 && a[i-2].fp1 == a[i].fp1)){
                    pthread_mutex_lock(&g_list_mtx);
                    fprintf(g_fdups, "%016llx %016llx\n", (unsigned long long)a[i].fp0, (unsigned long long)a[i].fp1);
                    pthread_mutex_unlock(&g_list_mtx);
                }
            } else {
                S.distinct++;
                if(same0){
                    S.coll0++;
                    pthread_mutex_lock(&g_list_mtx);
                    fprintf(g_fcoll0, "%016llx %016llx %016llx\n", (unsigned long long)a[i].fp0,
                            (unsigned long long)a[i-1].fp1, (unsigned long long)a[i].fp1);
                    pthread_mutex_unlock(&g_list_mtx);
                } else S.fp0_distinct++;
            }
        }
        g_bstat[b] = S;
        pthread_mutex_lock(&g_list_mtx);
        for(int q = 0; q < NBUCKET; q++) g_fp1top[q] += t1h[q];
        pthread_mutex_unlock(&g_list_mtx);
        /* write sorted back: write to temp then rename */
        if(write_all(ps, a, (size_t)sz) == 0) rename(ps, p);
        double t3 = now_sec();
        free(a); free(scr);
        fprintf(stderr, "[sort] bucket %02x n=%ld read %.1fs sort %.1fs write %.1fs  dups=%llu coll0=%llu\n",
                b, n, t1 - t0, t2 - t1, t3 - t2,
                (unsigned long long)S.dup_pairs, (unsigned long long)S.coll0);
    }
    return NULL;
}

static void phase_sort(void){
    char p[700];
    snprintf(p, sizeof p, "%s/../dups.txt", g_dir);  g_fdups  = fopen(p, "w");
    snprintf(p, sizeof p, "%s/../coll0.txt", g_dir); g_fcoll0 = fopen(p, "w");
    g_hist20 = (u64*)calloc(1u << 20, sizeof(u64));
    g_next_bucket = 0;
    double t0 = now_sec();
    pthread_t th[256];
    int P = g_sortpar;
    for(int i = 0; i < P; i++) pthread_create(&th[i], NULL, sort_worker, NULL);
    for(int i = 0; i < P; i++) pthread_join(th[i], NULL);
    double t1 = now_sec();
    fclose(g_fdups); fclose(g_fcoll0);

    /* aggregate */
    u64 N = 0, dups = 0, distinct = 0, coll0 = 0;
    u64 bmin = ~0ULL, bmax = 0;
    for(int b = 0; b < NBUCKET; b++){
        N += g_bstat[b].n; dups += g_bstat[b].dup_pairs; distinct += g_bstat[b].distinct; coll0 += g_bstat[b].coll0;
        if(g_bstat[b].n < bmin) bmin = g_bstat[b].n;
        if(g_bstat[b].n > bmax) bmax = g_bstat[b].n;
    }
    double ebuck = (double)N / NBUCKET;
    double chi2_b = 0;
    for(int b = 0; b < NBUCKET; b++){ double d = (double)g_bstat[b].n - ebuck; chi2_b += d * d / ebuck; }
    double e20 = (double)N / (double)(1u << 20);
    double chi2_20 = 0; u64 hmin = ~0ULL, hmax = 0;
    for(u64 i = 0; i < (1u << 20); i++){
        double d = (double)g_hist20[i] - e20; chi2_20 += d * d / e20;
        if(g_hist20[i] < hmin) hmin = g_hist20[i];
        if(g_hist20[i] > hmax) hmax = g_hist20[i];
    }
    double npairs = (double)N * (double)(N - 1) / 2.0;
    double exp_coll = npairs / 18446744073709551616.0;   /* /2^64 */
    LOGF("\n=== SORT/DEDUP RESULTS (%.0f s wall) ===\n", t1 - t0);
    LOGF("total records N = %llu (2^%.4f)\n", (unsigned long long)N, log2((double)N));
    LOGF("distinct (fp0,fp1) pairs = %llu ; duplicate records (full 128-bit pair equal) = %llu ; distinct duplicated values listed in dups.txt\n",
         (unsigned long long)distinct, (unsigned long long)dups);
    LOGF("fp0-word (64-bit) collisions with different fp1 = %llu   (random expectation C(N,2)/2^64 = %.1f)\n",
         (unsigned long long)coll0, exp_coll);
    LOGF("empirical fp0-word pairwise collision prob = %.3e  (ideal 2^-64 = %.3e)\n",
         (double)coll0 / npairs, 1.0 / 18446744073709551616.0);
    LOGF("bucket size (top byte of fp0): min %llu max %llu mean %.0f ; chi2 = %.1f (df 255, z = %.2f)\n",
         (unsigned long long)bmin, (unsigned long long)bmax, ebuck, chi2_b, (chi2_b - 255) / sqrt(2.0 * 255));
    LOGF("top-20-bit prefix histogram: min %llu max %llu mean %.1f ; chi2 = %.1f (df 2^20-1=1048575, z = %.2f)\n",
         (unsigned long long)hmin, (unsigned long long)hmax, e20, chi2_20, (chi2_20 - 1048575.0) / sqrt(2.0 * 1048575.0));
    /* fp1 top-byte distribution (fp1 is a MIN over candidate hashes -> skewed low) */
    {
        g_have_fp1top = 1;
        u64 t1min = ~0ULL, t1max = 0; double chi2 = 0;
        for(int q = 0; q < NBUCKET; q++){
            if(g_fp1top[q] < t1min) t1min = g_fp1top[q];
            if(g_fp1top[q] > t1max) t1max = g_fp1top[q];
            double d = (double)g_fp1top[q] - ebuck; chi2 += d * d / ebuck;
        }
        LOGF("fp1 top-byte histogram (fp1 = min of candidate hashes => NOT uniform): min %llu max %llu mean %.0f ; chi2 = %.1f (df 255)\n",
             (unsigned long long)t1min, (unsigned long long)t1max, ebuck, chi2);
        LOGF("fp1 top-byte counts by byte value: q=0x00:%llu 0x10:%llu 0x20:%llu 0x40:%llu 0x80:%llu 0xc0:%llu 0xe0:%llu 0xff:%llu\n",
             (unsigned long long)g_fp1top[0x00], (unsigned long long)g_fp1top[0x10], (unsigned long long)g_fp1top[0x20],
             (unsigned long long)g_fp1top[0x40], (unsigned long long)g_fp1top[0x80], (unsigned long long)g_fp1top[0xc0],
             (unsigned long long)g_fp1top[0xe0], (unsigned long long)g_fp1top[0xff]);
    }
    free(g_hist20); g_hist20 = NULL;
}

/* ------------------------------------------------------------------ */
/* PHASE FP1 stats: partition all fp1 values by top byte, sort, count  */
/* ------------------------------------------------------------------ */
static u64   *g_p1[NBUCKET];          /* partitions of fp1 values */
static long   g_p1n[NBUCKET], g_p1cap[NBUCKET], g_p1capmax;
static pthread_mutex_t g_p1mtx[NBUCKET];
static int    g_p1_countonly = 0;
static u64    g_p1cnt_tmp[NBUCKET];

static void *fp1_load_worker(void *arg){
    (void)arg;
    const long CH = 1 << 22;   /* 4M records = 64 MiB chunk */
    rec_t *buf = (rec_t*)malloc((size_t)CH * sizeof(rec_t));
    u64 *lb[NBUCKET]; int ln[NBUCKET];
    u64 lc[NBUCKET]; memset(lc, 0, sizeof lc);
    const int LCAP = 4096;
    if(!g_p1_countonly) for(int q = 0; q < NBUCKET; q++){ lb[q] = (u64*)malloc(LCAP * sizeof(u64)); ln[q] = 0; }
    for(;;){
        int b = __sync_fetch_and_add(&g_next_bucket, 1);
        if(b >= NBUCKET) break;
        char p[700]; bucket_path(p, sizeof p, b);
        int fd = open(p, O_RDONLY); if(fd < 0){ perror(p); continue; }
        off_t off = 0;
        for(;;){
            ssize_t r = pread(fd, buf, (size_t)CH * sizeof(rec_t), off);
            if(r <= 0) break;
            long cnt = (long)(r / (ssize_t)sizeof(rec_t));
            if(g_p1_countonly){
                for(long i = 0; i < cnt; i++) lc[buf[i].fp1 >> 56]++;
            } else {
                for(long i = 0; i < cnt; i++){
                    u64 v = buf[i].fp1; int q = (int)(v >> 56);
                    lb[q][ln[q]++] = v;
                    if(ln[q] == LCAP){
                        pthread_mutex_lock(&g_p1mtx[q]);
                        if(g_p1n[q] + LCAP > g_p1cap[q]){ fprintf(stderr, "FATAL fp1 partition %02x overflow\n", q); exit(9); }
                        memcpy(g_p1[q] + g_p1n[q], lb[q], LCAP * sizeof(u64)); g_p1n[q] += LCAP;
                        pthread_mutex_unlock(&g_p1mtx[q]);
                        ln[q] = 0;
                    }
                }
            }
            off += r;
            if(cnt < CH) break;
        }
        close(fd);
    }
    if(g_p1_countonly){
        pthread_mutex_lock(&g_list_mtx);
        for(int q = 0; q < NBUCKET; q++) g_p1cnt_tmp[q] += lc[q];
        pthread_mutex_unlock(&g_list_mtx);
    } else {
        for(int q = 0; q < NBUCKET; q++){
            if(ln[q]){
                pthread_mutex_lock(&g_p1mtx[q]);
                if(g_p1n[q] + ln[q] > g_p1cap[q]){ fprintf(stderr, "FATAL fp1 partition %02x overflow (tail)\n", q); exit(9); }
                memcpy(g_p1[q] + g_p1n[q], lb[q], (size_t)ln[q] * sizeof(u64)); g_p1n[q] += ln[q];
                pthread_mutex_unlock(&g_p1mtx[q]);
            }
            free(lb[q]);
        }
    }
    free(buf);
    return NULL;
}
static u64 g_coll1 = 0;
static FILE *g_fcoll1 = NULL;
static void *fp1_sort_worker(void *arg){
    (void)arg;
    u64 *scr = (u64*)malloc((size_t)g_p1capmax * sizeof(u64));
    u64 local = 0;
    for(;;){
        int q = __sync_fetch_and_add(&g_next_bucket, 1);
        if(q >= NBUCKET) break;
        sort_u64_part(g_p1[q], g_p1n[q], scr);
        for(long i = 1; i < g_p1n[q]; i++) if(g_p1[q][i] == g_p1[q][i-1]){
            local++;
            pthread_mutex_lock(&g_list_mtx);
            fprintf(g_fcoll1, "%016llx\n", (unsigned long long)g_p1[q][i]);
            pthread_mutex_unlock(&g_list_mtx);
        }
    }
    __sync_fetch_and_add(&g_coll1, local);
    free(scr);
    return NULL;
}
static int g_have_p1 = 0;
static void run_fp1_load(int countonly, int P){
    g_p1_countonly = countonly;
    g_next_bucket = 0;
    pthread_t th[256];
    if(P > 64) P = 64;
    for(int i = 0; i < P; i++) pthread_create(&th[i], NULL, fp1_load_worker, NULL);
    for(int i = 0; i < P; i++) pthread_join(th[i], NULL);
}
static void phase_fp1stats(u64 Ntotal){
    LOGF("\n=== FP1 WORD ANALYSIS (%llu values) ===\n", (unsigned long long)Ntotal);
    /* exact per-partition counts: from the sort phase, or a counting scan */
    if(!g_have_fp1top){
        memset(g_p1cnt_tmp, 0, sizeof g_p1cnt_tmp);
        double tc = now_sec();
        run_fp1_load(1, g_threads);
        for(int q = 0; q < NBUCKET; q++) g_fp1top[q] = g_p1cnt_tmp[q];
        g_have_fp1top = 1;
        LOGF("fp1 top-byte counting scan: %.0fs\n", now_sec() - tc);
        {
            u64 t1min = ~0ULL, t1max = 0; double mean = (double)Ntotal / NBUCKET, chi2 = 0;
            for(int q = 0; q < NBUCKET; q++){
                if(g_fp1top[q] < t1min) t1min = g_fp1top[q];
                if(g_fp1top[q] > t1max) t1max = g_fp1top[q];
                double d = (double)g_fp1top[q] - mean; chi2 += d * d / mean;
            }
            LOGF("fp1 top-byte histogram (fp1 = min of candidate hashes => NOT uniform): min %llu max %llu mean %.0f ; chi2 = %.1f (df 255)\n",
                 (unsigned long long)t1min, (unsigned long long)t1max, mean, chi2);
            LOGF("fp1 top-byte counts by byte value: q=0x00:%llu 0x10:%llu 0x20:%llu 0x40:%llu 0x80:%llu 0xc0:%llu 0xe0:%llu 0xff:%llu\n",
                 (unsigned long long)g_fp1top[0x00], (unsigned long long)g_fp1top[0x10], (unsigned long long)g_fp1top[0x20],
                 (unsigned long long)g_fp1top[0x40], (unsigned long long)g_fp1top[0x80], (unsigned long long)g_fp1top[0xc0],
                 (unsigned long long)g_fp1top[0xe0], (unsigned long long)g_fp1top[0xff]);
        }
    }
    g_p1capmax = 0;
    double gib = 0;
    for(int q = 0; q < NBUCKET; q++){
        g_p1cap[q] = (long)g_fp1top[q] + 16;
        if(g_p1cap[q] > g_p1capmax) g_p1capmax = g_p1cap[q];
        gib += (double)g_p1cap[q] * 8.0 / (1024.0*1024*1024);
        g_p1[q] = (u64*)malloc((size_t)g_p1cap[q] * sizeof(u64)); g_p1n[q] = 0;
        pthread_mutex_init(&g_p1mtx[q], NULL);
        if(!g_p1[q]){ fprintf(stderr, "malloc fail fp1 part\n"); exit(9); }
    }
    LOGF("allocated 256 exact-size fp1 partitions, %.1f GiB total\n", gib);
    double t0 = now_sec();
    run_fp1_load(0, g_threads);
    double t1 = now_sec();
    u64 n1 = 0; for(int q = 0; q < NBUCKET; q++) n1 += (u64)g_p1n[q];
    LOGF("fp1 values loaded: %llu in %.0fs\n", (unsigned long long)n1, t1 - t0);
    char p[700]; snprintf(p, sizeof p, "%s/../coll1.txt", g_dir); g_fcoll1 = fopen(p, "w");
    g_next_bucket = 0;
    {
        int P = g_threads;
        pthread_t th[512];
        for(int i = 0; i < P; i++) pthread_create(&th[i], NULL, fp1_sort_worker, NULL);
        for(int i = 0; i < P; i++) pthread_join(th[i], NULL);
    }
    fclose(g_fcoll1);
    double t2 = now_sec();
    double npairs = (double)n1 * (double)(n1 - 1) / 2.0;
    /* skew-adjusted expectation: sum over top-byte partitions of C(n_q,2)/2^56 */
    double exp_skew = 0;
    for(int q = 0; q < NBUCKET; q++){ double nq = (double)g_p1n[q]; exp_skew += nq * (nq - 1) / 2.0 / 72057594037927936.0; }
    LOGF("fp1-word (64-bit) collisions = %llu  (uniform-64-bit expectation %.1f; expectation given observed top-byte skew [uniform within partition] %.1f)  [sorted in %.0fs]\n",
         (unsigned long long)g_coll1, npairs / 18446744073709551616.0, exp_skew, t2 - t1);
    LOGF("empirical fp1-word pairwise collision prob = %.3e  (ideal 2^-64 = %.3e)\n",
         (double)g_coll1 / npairs, 1.0 / 18446744073709551616.0);
    g_have_p1 = 1;
}

/* ------------------------------------------------------------------ */
/* PHASE PROBE                                                           */
/* ------------------------------------------------------------------ */
static u64 *g_pr0 = NULL, *g_pr1 = NULL;          /* probe fH0, fH1 arrays, index = inst*G+g */
static u64 *g_pr0s = NULL, *g_pr1s = NULL;        /* partition-sorted copies */
static long g_pr_start[NBUCKET+1], g_pr1_start[NBUCKET+1];
static volatile u64 g_next_inst = 0;

/* one probe instance: real AES, genuine delta-set, wrong k7a guesses */
typedef struct {
    aesctx ctx;
    u8 P[16], km1[4], k7a[4];
    u8 Cd[256][16];
} pinst_t;

static void pinst_make(pinst_t *I, rng_t *rng){
    u8 key[16];
    for(int i = 0; i < 16; i++) key[i] = rng8(rng);
    aes_set_key(&I->ctx, key, NR);
    for(int i = 0; i < 16; i++) I->P[i] = rng8(rng);
    for(int i = 0; i < 4; i++) I->km1[i] = I->ctx.rk[0][DIAG0[i]];
    for(int i = 0; i < 4; i++) I->k7a[i] = I->ctx.rk[7][ADIAG[i]];
    /* genuine delta-set (true km1) */
    u8 z1r[4]; for(int i = 0; i < 4; i++) z1r[i] = SBOX[I->P[DIAG0[i]] ^ I->km1[i]];
    u8 w1r[4];
    for(int r = 0; r < 4; r++) w1r[r] = (u8)(gmul(MCc[r][0], z1r[0]) ^ gmul(MCc[r][1], z1r[1]) ^ gmul(MCc[r][2], z1r[2]) ^ gmul(MCc[r][3], z1r[3]));
    for(int dv = 0; dv < 256; dv++){
        u8 w1[4] = { (u8)(w1r[0] ^ dv), w1r[1], w1r[2], w1r[3] };
        u8 z1[4];
        for(int r = 0; r < 4; r++) z1[r] = (u8)(gmul(iMCc[r][0], w1[0]) ^ gmul(iMCc[r][1], w1[1]) ^ gmul(iMCc[r][2], w1[2]) ^ gmul(iMCc[r][3], w1[3]));
        u8 PP[16]; memcpy(PP, I->P, 16);
        for(int i = 0; i < 4; i++) PP[DIAG0[i]] = (u8)(iSBOX[z1[i]] ^ I->km1[i]);
        aes_enc(&I->ctx, PP, I->Cd[dv]);
    }
}
/* online fingerprint for a k7a guess: honest peel -> eH -> kernel */
static inline void probe_fp(const pinst_t *I, const u8 k7a[4], u64 *f0, u64 *f1){
    u8 eH[256] __attribute__((aligned(64)));
    eH[0] = 0;
    /* per byte-position tables: t_i[dv] = iMC[0][i] * iSBOX[Cd[dv][adiag i] ^ k7a[i]] */
    const u8 *am0 = MUL[iMCc[0][0]], *am1 = MUL[iMCc[0][1]], *am2 = MUL[iMCc[0][2]], *am3 = MUL[iMCc[0][3]];
    u8 k0 = k7a[0], k1 = k7a[1], k2 = k7a[2], k3 = k7a[3];
    const int p0 = ADIAG[0], p1 = ADIAG[1], p2 = ADIAG[2], p3 = ADIAG[3];
    u8 ref = (u8)(am0[iSBOX[I->Cd[0][p0] ^ k0]] ^ am1[iSBOX[I->Cd[0][p1] ^ k1]] ^ am2[iSBOX[I->Cd[0][p2] ^ k2]] ^ am3[iSBOX[I->Cd[0][p3] ^ k3]]);
    for(int dv = 1; dv < 256; dv++){
        const u8 *C = I->Cd[dv];
        u8 v = (u8)(am0[iSBOX[C[p0] ^ k0]] ^ am1[iSBOX[C[p1] ^ k1]] ^ am2[iSBOX[C[p2] ^ k2]] ^ am3[iSBOX[C[p3] ^ k3]]);
        eH[dv] = Linv[(u8)(v ^ ref)];
    }
    kernel_v2_vec(eH, ZROWS, ZROWS, ZROWS, f0, f1);
}
static void probe_guess_k7a(rng_t *rng, const u8 truek[4], u8 out[4]){
    do { for(int i = 0; i < 4; i++) out[i] = rng8(rng); } while(!memcmp(out, truek, 4));
}

typedef struct { int dummy; } probegen_arg_t;
static volatile u64 g_probe_done = 0;
/* resolve mode: list of target values to identify */
static u64 *g_tv0 = NULL, *g_tv1 = NULL; static long g_ntv0 = 0, g_ntv1 = 0;
static int tv_has(const u64 *A, long n, u64 v){ long lo = 0, hi = n - 1; while(lo <= hi){ long m = (lo+hi)>>1; if(A[m] == v) return 1; if(A[m] < v) lo = m+1; else hi = m-1; } return 0; }
static FILE *g_fprobe_resolve = NULL;
static int g_probe_mode = 0;  /* 0 = fill arrays, 1 = resolve targets */

static void *probegen_worker(void *arg){
    (void)arg;
    pinst_t I; rng_t rng, rg;
    for(;;){
        u64 j = __sync_fetch_and_add(&g_next_inst, 1);
        if(j >= g_probe_inst) break;
        rng_seed2(&rng, g_seed ^ 0x50524f4245ULL, j);
        pinst_make(&I, &rng);
        rng_seed2(&rg, g_seed ^ 0x4755455353ULL, j);
        u64 base_idx = j * g_probe_g;
        for(u64 g = 0; g < g_probe_g; g++){
            u8 k7a[4]; probe_guess_k7a(&rg, I.k7a, k7a);
            u64 f0, f1; probe_fp(&I, k7a, &f0, &f1);
            if(g_probe_mode == 0){
                g_pr0[base_idx + g] = f0; g_pr1[base_idx + g] = f1;
            } else {
                if(tv_has(g_tv0, g_ntv0, f0) || tv_has(g_tv1, g_ntv1, f1)){
                    pthread_mutex_lock(&g_list_mtx);
                    fprintf(g_fprobe_resolve, "PROBE inst=%llu g=%llu k7a=%02x%02x%02x%02x truek7a=%02x%02x%02x%02x f0=%016llx f1=%016llx\n",
                            (unsigned long long)j, (unsigned long long)g, k7a[0], k7a[1], k7a[2], k7a[3],
                            I.k7a[0], I.k7a[1], I.k7a[2], I.k7a[3], (unsigned long long)f0, (unsigned long long)f1);
                    fflush(g_fprobe_resolve);
                    pthread_mutex_unlock(&g_list_mtx);
                }
            }
        }
        __sync_fetch_and_add(&g_probe_done, g_probe_g);
    }
    return NULL;
}
/* partition+sort probe arrays */
static u64 *g_psrc, *g_pdst; static long *g_pstart;
static void *probesort_part_worker(void *arg){
    (void)arg;
    long cap = 0; for(int q = 0; q < NBUCKET; q++){ long len = g_pstart[q+1] - g_pstart[q]; if(len > cap) cap = len; }
    u64 *scr = (u64*)malloc((size_t)cap * sizeof(u64));
    for(;;){
        int q = __sync_fetch_and_add(&g_next_bucket, 1);
        if(q >= NBUCKET) break;
        sort_u64_part(g_pdst + g_pstart[q], g_pstart[q+1] - g_pstart[q], scr);
    }
    free(scr);
    return NULL;
}
static void probe_partition_sort(u64 *src, u64 *dst, long n, long *start){
    /* counting pass */
    long cnt[NBUCKET+1]; memset(cnt, 0, sizeof cnt);
    for(long i = 0; i < n; i++) cnt[(src[i] >> 56) + 1]++;
    for(int q = 0; q < NBUCKET; q++) cnt[q+1] += cnt[q];
    memcpy(start, cnt, sizeof(long) * (NBUCKET + 1));
    long *pos = (long*)malloc(sizeof(long) * NBUCKET);
    memcpy(pos, cnt, sizeof(long) * NBUCKET);
    for(long i = 0; i < n; i++) dst[pos[src[i] >> 56]++] = src[i];
    free(pos);
    g_psrc = src; g_pdst = dst; g_pstart = start;
    g_next_bucket = 0;
    int P = g_threads;
    pthread_t th[512];
    for(int i = 0; i < P; i++) pthread_create(&th[i], NULL, probesort_part_worker, NULL);
    for(int i = 0; i < P; i++) pthread_join(th[i], NULL);
}

/* prefix-level pair counting between two sorted arrays sharing top byte */
static const int LEV[5] = {32, 40, 48, 56, 64};
static void count_prefix_pairs(const u64 *A, long na, const u64 *B, long nb, double out[5], u64 *exact_vals, int *nexact, int cap_exact){
    for(int L = 0; L < 5; L++){
        int k = LEV[L];
        u64 mask = (k == 64) ? ~0ULL : ~((1ULL << (64 - k)) - 1ULL);
        long i = 0, j = 0; double tot = 0;
        while(i < na && j < nb){
            u64 ga = A[i] & mask, gb = B[j] & mask;
            if(ga < gb){ do i++; while(i < na && (A[i] & mask) == ga); continue; }
            if(gb < ga){ do j++; while(j < nb && (B[j] & mask) == gb); continue; }
            long i2 = i; while(i2 < na && (A[i2] & mask) == ga) i2++;
            long j2 = j; while(j2 < nb && (B[j2] & mask) == gb) j2++;
            tot += (double)(i2 - i) * (double)(j2 - j);
            if(k == 64 && exact_vals && *nexact < cap_exact) exact_vals[(*nexact)++] = ga;
            else if(k == 64 && exact_vals) (*nexact)++;
            i = i2; j = j2;
        }
        out[L] = tot;
    }
}

/* J1: table fp1 partitions (RAM) vs probe arr1 partitions */
static double g_j1[5]; static u64 g_j1_vals[10000]; static int g_j1_nv = 0;
static void *j1_worker(void *arg){
    (void)arg;
    double acc[5] = {0,0,0,0,0};
    for(;;){
        int q = __sync_fetch_and_add(&g_next_bucket, 1);
        if(q >= NBUCKET) break;
        double out[5]; u64 vals[256]; int nv = 0;
        count_prefix_pairs(g_p1[q], g_p1n[q], g_pr1s + g_pr1_start[q], g_pr1_start[q+1] - g_pr1_start[q], out, vals, &nv, 256);
        for(int L = 0; L < 5; L++) acc[L] += out[L];
        if(nv){
            pthread_mutex_lock(&g_list_mtx);
            for(int t = 0; t < nv && t < 256; t++) if(g_j1_nv < 10000) g_j1_vals[g_j1_nv++] = vals[t];
            pthread_mutex_unlock(&g_list_mtx);
        }
    }
    pthread_mutex_lock(&g_list_mtx);
    for(int L = 0; L < 5; L++) g_j1[L] += acc[L];
    pthread_mutex_unlock(&g_list_mtx);
    return NULL;
}
/* J2: table sorted bucket files (fp0) vs probe arr0 partitions; also capture
 * records whose fp1 is in the J1 exact-match value set */
static double g_j2[5];
typedef struct { u64 fp0, fp1; int bucket; long pos; int kind; } cap_t;   /* kind 0=fp0 match, 1=fp1 match */
static cap_t g_caps[20000]; static int g_ncaps = 0;
static u64 g_j1sorted[10000]; static int g_j1n_sorted = 0;
static void *j2_worker(void *arg){
    (void)arg;
    const long CH = 1 << 22;
    rec_t *buf = (rec_t*)malloc((size_t)CH * sizeof(rec_t));
    double acc[5] = {0,0,0,0,0};
    for(;;){
        int b = __sync_fetch_and_add(&g_next_bucket, 1);
        if(b >= NBUCKET) break;
        char p[700]; bucket_path(p, sizeof p, b);
        off_t sz = file_size(p); long n = (long)(sz / (off_t)sizeof(rec_t));
        rec_t *a = (rec_t*)malloc((size_t)sz + 16);
        if(read_all(p, a, (size_t)sz)){ free(a); continue; }
        const u64 *B = g_pr0s + g_pr_start[b]; long nb = g_pr_start[b+1] - g_pr_start[b];
        /* prefix-level counts over fp0 (A is a->fp0 sorted) */
        for(int L = 0; L < 5; L++){
            int k = LEV[L];
            u64 mask = (k == 64) ? ~0ULL : ~((1ULL << (64 - k)) - 1ULL);
            long i = 0, j = 0; double tot = 0;
            while(i < n && j < nb){
                u64 ga = a[i].fp0 & mask, gb = B[j] & mask;
                if(ga < gb){ do i++; while(i < n && (a[i].fp0 & mask) == ga); continue; }
                if(gb < ga){ do j++; while(j < nb && (B[j] & mask) == gb); continue; }
                long i2 = i; while(i2 < n && (a[i2].fp0 & mask) == ga) i2++;
                long j2 = j; while(j2 < nb && (B[j2] & mask) == gb) j2++;
                tot += (double)(i2 - i) * (double)(j2 - j);
                if(k == 64){
                    pthread_mutex_lock(&g_list_mtx);
                    for(long t = i; t < i2 && g_ncaps < 20000; t++){
                        g_caps[g_ncaps].fp0 = a[t].fp0; g_caps[g_ncaps].fp1 = a[t].fp1;
                        g_caps[g_ncaps].bucket = b; g_caps[g_ncaps].pos = t; g_caps[g_ncaps].kind = 0; g_ncaps++;
                    }
                    pthread_mutex_unlock(&g_list_mtx);
                }
                i = i2; j = j2;
            }
            acc[L] += tot;
        }
        /* capture fp1-side exact matches */
        if(g_j1n_sorted){
            for(long i = 0; i < n; i++){
                if(tv_has(g_j1sorted, g_j1n_sorted, a[i].fp1)){
                    pthread_mutex_lock(&g_list_mtx);
                    if(g_ncaps < 20000){
                        g_caps[g_ncaps].fp0 = a[i].fp0; g_caps[g_ncaps].fp1 = a[i].fp1;
                        g_caps[g_ncaps].bucket = b; g_caps[g_ncaps].pos = i; g_caps[g_ncaps].kind = 1; g_ncaps++;
                    }
                    pthread_mutex_unlock(&g_list_mtx);
                }
            }
        }
        free(a);
    }
    pthread_mutex_lock(&g_list_mtx);
    for(int L = 0; L < 5; L++) g_j2[L] += acc[L];
    pthread_mutex_unlock(&g_list_mtx);
    free(buf);
    return NULL;
}
static int cmp_u64q(const void *a, const void *b){ u64 x = *(const u64*)a, y = *(const u64*)b; return x < y ? -1 : x > y; }

static void phase_probe(u64 Ntable){
    u64 NP = g_probe_inst * g_probe_g;
    LOGF("\n=== PROBE TEST: %llu genuine online wrong-candidate fingerprints (%llu instances x %llu wrong k7a guesses) ===\n",
         (unsigned long long)NP, (unsigned long long)g_probe_inst, (unsigned long long)g_probe_g);
    g_pr0 = (u64*)malloc(NP * sizeof(u64)); g_pr1 = (u64*)malloc(NP * sizeof(u64));
    g_pr0s = (u64*)malloc(NP * sizeof(u64)); g_pr1s = (u64*)malloc(NP * sizeof(u64));
    if(!g_pr0 || !g_pr1 || !g_pr0s || !g_pr1s){ fprintf(stderr, "malloc probes\n"); exit(9); }
    double t0 = now_sec();
    g_next_inst = 0; g_probe_done = 0; g_probe_mode = 0;
    {
        pthread_t th[512];
        for(int i = 0; i < g_threads; i++) pthread_create(&th[i], NULL, probegen_worker, NULL);
        for(int i = 0; i < g_threads; i++) pthread_join(th[i], NULL);
    }
    double t1 = now_sec();
    LOGF("probe generation: %llu probes in %.0fs (%.1f Mprobe/s aggregate, %.3f Mprobe/s/core; ~%.0f TSC-ticks/probe incl. online peel + kernel)\n",
         (unsigned long long)NP, t1 - t0, NP / (t1 - t0) / 1e6, NP / (t1 - t0) / 1e6 / g_threads,
         g_tsc_hz * (t1 - t0) * g_threads / (double)NP);
    /* partition + sort both arrays */
    probe_partition_sort(g_pr0, g_pr0s, (long)NP, g_pr_start);
    probe_partition_sort(g_pr1, g_pr1s, (long)NP, g_pr1_start);
    free(g_pr0); g_pr0 = NULL; free(g_pr1); g_pr1 = NULL;
    double t2 = now_sec();
    LOGF("probe arrays partitioned+sorted in %.0fs\n", t2 - t1);
    /* within-probe-set word self-collision sanity (both sorted): */
    {
        u64 d0 = 0, d1 = 0;
        for(int q = 0; q < NBUCKET; q++){
            const u64 *A = g_pr0s + g_pr_start[q]; long na = g_pr_start[q+1] - g_pr_start[q];
            for(long i = 1; i < na; i++) if(A[i] == A[i-1]) d0++;
            const u64 *C = g_pr1s + g_pr1_start[q]; long nc = g_pr1_start[q+1] - g_pr1_start[q];
            for(long i = 1; i < nc; i++) if(C[i] == C[i-1]) d1++;
        }
        double np = (double)NP * (double)(NP - 1) / 2.0;
        double skew1 = 0;
        u64 pmin = ~0ULL, pmax = 0;
        for(int q = 0; q < NBUCKET; q++){
            double nq = (double)(g_pr1_start[q+1] - g_pr1_start[q]);
            skew1 += nq * (nq - 1) / 2.0 / 72057594037927936.0;
            u64 n0q = (u64)(g_pr_start[q+1] - g_pr_start[q]);
            if(n0q < pmin) pmin = n0q;
            if(n0q > pmax) pmax = n0q;
        }
        LOGF("probe fH0 top-byte partition sizes: min %llu max %llu (uniform %.0f); fH1 top byte skewed like fp1 (min-of-hashes)\n",
             (unsigned long long)pmin, (unsigned long long)pmax, (double)NP / NBUCKET);
        LOGF("within-probe 64-bit duplicates: fH0 %llu (uniform expectation %.1f), fH1 %llu (uniform expectation %.1f, topbyte-skew expectation %.1f)\n",
             (unsigned long long)d0, np / 18446744073709551616.0, (unsigned long long)d1, np / 18446744073709551616.0, skew1);
    }
    double exp64 = (double)Ntable * (double)NP / 18446744073709551616.0;
    /* J1: fp1 side (needs table fp1 partitions in RAM) */
    if(g_have_p1){
        memset(g_j1, 0, sizeof g_j1); g_j1_nv = 0;
        g_next_bucket = 0;
        pthread_t th[512];
        for(int i = 0; i < g_threads; i++) pthread_create(&th[i], NULL, j1_worker, NULL);
        for(int i = 0; i < g_threads; i++) pthread_join(th[i], NULL);
        /* skew-adjusted expectations: both fp1 and fH1 are minima of candidate
         * hashes, so their top bytes are non-uniform; model uniform below the
         * top byte: pairs with >=k common bits ~ sum_q T_q*P_q / 2^(k-8). */
        double skewsum = 0;
        for(int q = 0; q < NBUCKET; q++) skewsum += (double)g_p1n[q] * (double)(g_pr1_start[q+1] - g_pr1_start[q]);
        LOGF("J1 (table fp1 vs probe fH1): prefix-level matching pairs  [fp1/fH1 are MIN-of-candidates => top byte skewed; two expectation models]\n");
        for(int L = 0; L < 5; L++)
            LOGF("   >=%2d common bits: observed %.4g   uniform-model expectation %.4g (ratio %.3f)   topbyte-skew-model expectation %.4g (ratio %.3f)\n",
                 LEV[L], g_j1[L],
                 (double)Ntable * (double)NP / ldexp(1.0, LEV[L]), g_j1[L] / ((double)Ntable * (double)NP / ldexp(1.0, LEV[L])),
                 skewsum / ldexp(1.0, LEV[L] - 8), g_j1[L] / (skewsum / ldexp(1.0, LEV[L] - 8)));
        LOGF("   exact 64-bit fp1-word matches: %d distinct values captured\n", g_j1_nv);
        for(int i = 0; i < g_j1_nv; i++) g_j1sorted[i] = g_j1_vals[i];
        g_j1n_sorted = g_j1_nv > 10000 ? 10000 : g_j1_nv;
        qsort(g_j1sorted, (size_t)g_j1n_sorted, sizeof(u64), cmp_u64q);
        /* free table fp1 partitions */
        for(int q = 0; q < NBUCKET; q++){ free(g_p1[q]); g_p1[q] = NULL; }
        g_have_p1 = 0;
    } else LOGF("J1 skipped (no fp1 partitions in RAM)\n");
    /* J2: fp0 side against sorted bucket files */
    double t3 = now_sec();
    memset(g_j2, 0, sizeof g_j2); g_ncaps = 0;
    g_next_bucket = 0;
    {
        int P = g_threads < 48 ? g_threads : 48;
        pthread_t th[512];
        for(int i = 0; i < P; i++) pthread_create(&th[i], NULL, j2_worker, NULL);
        for(int i = 0; i < P; i++) pthread_join(th[i], NULL);
    }
    double t4 = now_sec();
    LOGF("J2 (table fp0 [sorted bucket files, binary/merge search] vs probe fH0): %.0fs\n", t4 - t3);
    for(int L = 0; L < 5; L++)
        LOGF("   >=%2d common bits: observed %.4g   expected %.4g   ratio %.3f\n", LEV[L], g_j2[L],
             (double)Ntable * (double)NP / ldexp(1.0, LEV[L]), g_j2[L] / ((double)Ntable * (double)NP / ldexp(1.0, LEV[L])));
    LOGF("exact 64-bit matches: fp0-word %.0f, fp1-word %.0f  (each expected %.1f for an ideal uniform 64-bit word)\n",
         g_j2[4], g_j1[4], exp64);
    LOGF("captured matching table records: %d\n", g_ncaps);
    /* J3: resolve probe provenance for all matching values */
    {
        long n0 = 0, n1 = 0;
        u64 *tv0 = (u64*)malloc(sizeof(u64) * (size_t)(g_ncaps + 1));
        u64 *tv1 = (u64*)malloc(sizeof(u64) * (size_t)(g_ncaps + 1));
        for(int i = 0; i < g_ncaps; i++){ if(g_caps[i].kind == 0) tv0[n0++] = g_caps[i].fp0; else tv1[n1++] = g_caps[i].fp1; }
        qsort(tv0, (size_t)n0, sizeof(u64), cmp_u64q); qsort(tv1, (size_t)n1, sizeof(u64), cmp_u64q);
        /* dedup */
        long m0 = 0, m1 = 0;
        for(long i = 0; i < n0; i++) if(!i || tv0[i] != tv0[i-1]) tv0[m0++] = tv0[i];
        for(long i = 0; i < n1; i++) if(!i || tv1[i] != tv1[i-1]) tv1[m1++] = tv1[i];
        g_tv0 = tv0; g_ntv0 = m0; g_tv1 = tv1; g_ntv1 = m1;
        char pp[700]; snprintf(pp, sizeof pp, "%s/../probe_resolve.txt", g_dir);
        g_fprobe_resolve = fopen(pp, "w");
        /* also dump the captured table records */
        for(int i = 0; i < g_ncaps; i++)
            fprintf(g_fprobe_resolve, "TABREC kind=%d bucket=%02x pos=%ld fp0=%016llx fp1=%016llx\n",
                    g_caps[i].kind, g_caps[i].bucket, g_caps[i].pos,
                    (unsigned long long)g_caps[i].fp0, (unsigned long long)g_caps[i].fp1);
        double t5 = now_sec();
        g_next_inst = 0; g_probe_mode = 1;
        pthread_t th[512];
        for(int i = 0; i < g_threads; i++) pthread_create(&th[i], NULL, probegen_worker, NULL);
        for(int i = 0; i < g_threads; i++) pthread_join(th[i], NULL);
        fclose(g_fprobe_resolve); g_fprobe_resolve = NULL;
        double t6 = now_sec();
        LOGF("J3 probe provenance regeneration for %ld fp0-values + %ld fp1-values: %.0fs -> probe_resolve.txt\n", m0, m1, t6 - t5);
        /* read back probe lines to classify joint matches */
        FILE *f = fopen(pp, "r");
        char line[512];
        typedef struct { u64 f0, f1; unsigned long long inst, g; } prl_t;
        prl_t *prl = (prl_t*)malloc(sizeof(prl_t) * 4096); int npr = 0;
        while(f && fgets(line, sizeof line, f)){
            if(strncmp(line, "PROBE", 5)) continue;
            prl_t x; unsigned long long a, b2;
            char *q;
            if((q = strstr(line, "inst=")) == NULL) continue; x.inst = strtoull(q + 5, 0, 10);
            if((q = strstr(line, " g=")) == NULL) continue;    x.g = strtoull(q + 3, 0, 10);
            if((q = strstr(line, "f0=")) == NULL) continue;    a = strtoull(q + 3, 0, 16); x.f0 = a;
            if((q = strstr(line, "f1=")) == NULL) continue;    b2 = strtoull(q + 3, 0, 16); x.f1 = b2;
            if(npr < 4096) prl[npr++] = x;
        }
        if(f) fclose(f);
        int joint = 0, w0 = 0, w1 = 0;
        for(int i = 0; i < g_ncaps; i++){
            for(int t = 0; t < npr; t++){
                int m0x = (g_caps[i].fp0 == prl[t].f0), m1x = (g_caps[i].fp1 == prl[t].f1);
                if(m0x && m1x) joint++;
                if(g_caps[i].kind == 0 && m0x) w0++;
                if(g_caps[i].kind == 1 && m1x) w1++;
            }
        }
        LOGF("probe/table word-level matches resolved: fp0-side %d, fp1-side %d;  JOINT (same probe matches BOTH words of the same entry) = %d\n",
             w0, w1, joint);
        LOGF("EITHER-rule hits (the real attack's acceptance rule with 8-byte words): %.0f (fp0) + %.0f (fp1) candidate/entry pairs; 128-bit joint hits = %d\n",
             g_j2[4], g_j1[4], joint);
        LOGF("model: joint 128-bit hit expectation |T|*|P|/2^128 = %.3e;  96-bit-key model |T|*|P|/2^96 = %.3e;  64-bit-word model |T|*|P|/2^64 = %.1f per word\n",
             (double)Ntable * (double)NP / ldexp(1.0, 128), (double)Ntable * (double)NP / ldexp(1.0, 96), exp64);
        free(prl);
    }
    free(g_pr0s); g_pr0s = NULL; free(g_pr1s); g_pr1s = NULL;
}

/* ------------------------------------------------------------------ */
/* PHASE KA: known-answer coverage                                       */
/* ------------------------------------------------------------------ */
static void phase_ka(void){
    char p[700]; snprintf(p, sizeof p, "%s/../ka.bin", g_dir);
    FILE *f = fopen(p, "rb"); if(!f){ LOGF("KA: cannot open %s\n", p); return; }
    int n; if(fread(&n, sizeof(int), 1, f) != 1){ fclose(f); return; }
    ka_inst_t *K = (ka_inst_t*)malloc(sizeof(ka_inst_t) * (size_t)n);
    if(fread(K, sizeof(ka_inst_t), (size_t)n, f) != (size_t)n){ LOGF("KA read short\n"); fclose(f); return; }
    fclose(f);
    LOGF("\n=== KNOWN-ANSWER COVERAGE: %d planted true instances (real AES-128, 7 rounds) ===\n", n);
    int present = 0, either = 0, m0 = 0, m1 = 0, recovered = 0, azero = 0, found_by_online = 0, br = 0;
    for(int i = 0; i < n; i++){
        ka_inst_t *k = &K[i];
        if(k->a_ref == 0) azero++;
        if(k->bridge_ok >= 253) br++;
        /* recompute online fingerprint fresh from the stored key/P (independent recomputation) */
        rng_t dummy; (void)dummy;
        ka_inst_t R2 = *k;
        /* presence: binary search (fD0,fD1) in sorted bucket file */
        int b = (int)(k->fD0 >> 56);
        char bp[700]; bucket_path(bp, sizeof bp, b);
        int fd = open(bp, O_RDONLY);
        long nrec = (long)(file_size(bp) / (off_t)sizeof(rec_t));
        int pres = 0; long posfound = -1;
        if(fd >= 0){
            long lo = disk_lower_bound_fp0(fd, nrec, k->fD0);
            rec_t r;
            for(long t = lo; t < nrec && t < lo + 64; t++){
                if(pread(fd, &r, sizeof r, (off_t)t * (off_t)sizeof r) != (ssize_t)sizeof r) break;
                if(r.fp0 != k->fD0) break;
                if(r.fp1 == k->fD1){ pres = 1; posfound = t; break; }
            }
        }
        if(pres) present++;
        int e0 = (k->fD0 == k->fH0), e1 = (k->fD1 == k->fH1);
        if(e0) m0++;
        if(e1) m1++;
        if(e0 || e1) either++;
        /* attack-style lookup: online computes (fH0,fH1); look up fH0 in the fp0 index */
        int hit = 0;
        if(fd >= 0 && e0){
            long lo = disk_lower_bound_fp0(fd, nrec, k->fH0);
            rec_t r;
            for(long t = lo; t < nrec && t < lo + 64; t++){
                if(pread(fd, &r, sizeof r, (off_t)t * (off_t)sizeof r) != (ssize_t)sizeof r) break;
                if(r.fp0 != k->fH0) break;
                if(r.fp0 == k->fD0 && r.fp1 == k->fD1){ hit = 1; break; }
            }
        }
        if(hit) found_by_online++;
        if(pres && (e0 || e1)) recovered++;
        if(fd >= 0) close(fd);
        if(i < 8 || !(pres && (e0||e1)))
            LOGF("  ka[%3d] a=%02x bridge=%3d/255 fD0=%016llx fH0=%016llx fD1=%016llx fH1=%016llx present=%d(pos %ld) match0=%d match1=%d lookup_hit=%d\n",
                 i, k->a_ref, k->bridge_ok, (unsigned long long)k->fD0, (unsigned long long)k->fH0,
                 (unsigned long long)k->fD1, (unsigned long long)k->fH1, pres, posfound, e0, e1, hit);
        (void)R2;
    }
    LOGF("KA summary: instances=%d  a_ref==0 (bridge-degenerate)=%d  bridge>=253/255=%d\n", n, azero, br);
    LOGF("KA summary: true entry present in sorted table=%d/%d ; genuine online fingerprint matches its entry: fp0-word %d, fp1-word %d, EITHER %d/%d\n",
         present, n, m0, m1, either, n);
    LOGF("KA summary: RECOVERED (present AND genuine online EITHER-match) = %d/%d ; found by fp0-index binary search (raw-parity case) = %d\n",
         recovered, n, found_by_online);
    free(K);
}

/* ------------------------------------------------------------------ */
int main(int argc, char **argv){
    const char *mode = argc > 1 ? argv[1] : "all";
    for(int i = 2; i < argc; i++){
        if(!strncmp(argv[i], "threads=", 8)) g_threads = atoi(argv[i]+8);
        else if(!strncmp(argv[i], "dir=", 4)) snprintf(g_dir, sizeof g_dir, "%s", argv[i]+4);
        else if(!strncmp(argv[i], "sortpar=", 8)) g_sortpar = atoi(argv[i]+8);
        else if(!strncmp(argv[i], "pinst=", 6)) g_probe_inst = strtoull(argv[i]+6, 0, 0);
        else if(!strncmp(argv[i], "pg=", 3)) g_probe_g = strtoull(argv[i]+3, 0, 0);
        else if(!strncmp(argv[i], "seed=", 5)) g_seed = strtoull(argv[i]+5, 0, 0);
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }
    double ti = now_sec();
    aesbench_init();
    g_tsc_hz = tsc_hz();
    char lp[700]; snprintf(lp, sizeof lp, "%s/../analysis.log", g_dir);
    g_log = fopen(lp, "a");
    LOGF("[analyze] mode=%s init %.1fs TSC %.3fGHz threads=%d dir=%s\n", mode, now_sec() - ti, g_tsc_hz / 1e9, g_threads, g_dir);

    /* table size from bucket files */
    u64 Ntable = 0;
    for(int b = 0; b < NBUCKET; b++){ char p[700]; bucket_path(p, sizeof p, b); off_t s = file_size(p); if(s > 0) Ntable += (u64)s / sizeof(rec_t); }
    LOGF("[analyze] table records on disk: %llu (2^%.4f), %.1f GiB\n", (unsigned long long)Ntable,
         log2((double)(Ntable ? Ntable : 1)), (double)Ntable * 16.0 / (1024.0*1024*1024));

    int do_sort = !strcmp(mode, "sort") || !strcmp(mode, "all");
    int do_fp1  = !strcmp(mode, "fp1stats") || !strcmp(mode, "all") || !strcmp(mode, "probe");
    int do_probe= !strcmp(mode, "probe") || !strcmp(mode, "all");
    int do_ka   = !strcmp(mode, "ka") || !strcmp(mode, "all");
    if(do_sort) phase_sort();
    if(do_fp1) phase_fp1stats(Ntable);
    if(do_probe) phase_probe(Ntable);
    if(do_ka) phase_ka();
    LOGF("[analyze] done.\n");
    return 0;
}
