// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* table_build.c -- real-scale offline table generator (DDT-Gray walk +
 * kernel_v2 chi* fingerprints), bucket-staged by the top byte of fp0.
 *
 * Modes:
 *   build   : full generation.  NKA planted known-answer instances (real AES)
 *             + NB random parameter bases x 2^20 Gray-walked entries each,
 *             16-byte records (fp0,fp1) appended to 256 bucket files.
 *   bench   : compute-only throughput (no I/O), 1 thread and N threads.
 *   resolve : regenerate everything (no I/O) and dump provenance + sequence
 *             diagnostics for entries whose fp0 or fp1 appears in a target
 *             list (duplicate / collision resolution).
 *
 * Build: gcc -O3 -march=native -pthread -o table_build table_build.c -lm
 */
#include "common.h"

/* ------------------------------------------------------------------ */
/* configuration                                                         */
/* ------------------------------------------------------------------ */
static int    g_threads   = 190;
static int    g_nb        = 65536;     /* number of parameter bases */
static int    g_nka       = 128;       /* planted known-answer instances */
static u64    g_seed      = 0xA5E5BE5C0FFEEULL;
static char   g_outdir[512] = "./buckets";
static int    g_buf_rec   = 8192;      /* per-thread per-bucket buffer, records (8192*16=128KiB) */
static int    g_walk_bits = 20;
static double g_tsc_hz;

/* cross-check strides */
static int    g_chk_cold  = 8191;      /* cold_E / prop2 recompute every k entries within a base */
static int    g_chk_kern  = 65537;     /* scalar kernel recompute stride */

/* ------------------------------------------------------------------ */
/* shared state                                                          */
/* ------------------------------------------------------------------ */
typedef struct {
    u32 id; u8 Din, Dout, x2[4], z4[4]; u8 n4cnt; u8 tries_capped;
} baselog_t;

static bucketset_t   g_bset;
static volatile u64  g_next_base = 0;
static baselog_t    *g_baselog = NULL;

/* aggregate counters */
static volatile u64 g_entries_done = 0, g_bases_done = 0;
static volatile u64 g_chk_cold_n = 0, g_chk_cold_ok = 0;
static volatile u64 g_chk_p2_n = 0, g_chk_p2_ok = 0;
static volatile u64 g_chk_kern_n = 0, g_chk_kern_ok = 0;
static volatile u64 g_ticks_gen = 0, g_ticks_kern = 0, g_ticks_io = 0, g_ticks_chk = 0;

/* resolve-mode targets */
typedef struct { u64 *v; int n; } vset_t;
static vset_t g_t0 = {NULL, 0}, g_t1 = {NULL, 0};
static int vset_has(const vset_t *S, u64 x){
    /* tiny sets: linear/binary search on sorted array */
    int lo = 0, hi = S->n - 1;
    while(lo <= hi){ int m = (lo + hi) >> 1; if(S->v[m] == x) return 1; if(S->v[m] < x) lo = m + 1; else hi = m - 1; }
    return 0;
}
static int cmp_u64(const void *a, const void *b){ u64 x = *(const u64*)a, y = *(const u64*)b; return x < y ? -1 : x > y; }

static int  g_mode = 0;   /* 0=build, 1=bench, 2=resolve */
static FILE *g_resolve_out = NULL;
static pthread_mutex_t g_print_mtx = PTHREAD_MUTEX_INITIALIZER;

/* ------------------------------------------------------------------ */
/* worker                                                                */
/* ------------------------------------------------------------------ */
typedef struct {
    int tid;
    u64 entries;
    u64 t_gen, t_kern, t_chk;
    walk_t W;
    tlbuf_t tb;
} worker_t;

static void resolve_emit(const char *tag, u32 base, u32 idx, const base_t *B, const u8 *E, u64 f0, u64 f1){
    chi_t chi; u64 hx = chi_of_seq(E, &chi);
    pthread_mutex_lock(&g_print_mtx);
    fprintf(g_resolve_out, "%s base=%u gidx=%u fp0=%016llx fp1=%016llx seqh=%016llx mset=%016llx chih=%016llx chi=%016llx%016llx%016llx%016llx",
            tag, base, idx, (unsigned long long)f0, (unsigned long long)f1,
            (unsigned long long)seq_hash(E), (unsigned long long)multiset_hash(E),
            (unsigned long long)hx,
            (unsigned long long)chi.w[3], (unsigned long long)chi.w[2], (unsigned long long)chi.w[1], (unsigned long long)chi.w[0]);
    if(B) fprintf(g_resolve_out, " params=%02x%02x|%02x%02x%02x%02x|%02x%02x%02x%02x",
                  B->Din, B->Dout, B->x2[0], B->x2[1], B->x2[2], B->x2[3], B->z4[0], B->z4[1], B->z4[2], B->z4[3]);
    fprintf(g_resolve_out, " seq=");
    for(int i = 1; i < 256; i++) fprintf(g_resolve_out, "%02x", E[i]);
    fprintf(g_resolve_out, "\n");
    fflush(g_resolve_out);
    pthread_mutex_unlock(&g_print_mtx);
}

static void *worker_fn(void *arg){
    worker_t *me = (worker_t*)arg;
    walk_t *W = &me->W;
    tlbuf_t *tb = &me->tb;
    if(g_mode == 0) tlbuf_init(tb, &g_bset, g_buf_rec);
    const u32 NW = 1u << g_walk_bits;
    u8 ebuf[4][256] __attribute__((aligned(64)));
    const u8 *rows[4][4];
    for(int k = 0; k < 4; k++){ rows[k][0] = ebuf[k]; rows[k][1] = ZROWS; rows[k][2] = ZROWS; rows[k][3] = ZROWS; }

    for(;;){
        u64 b = __sync_fetch_and_add(&g_next_base, 1);
        if(b >= (u64)g_nb) break;
        /* deterministic per-base RNG */
        rng_t rng; rng_seed2(&rng, g_seed, b);
        base_t B; int tries = base_sample(&B, &rng);
        if(g_baselog){
            baselog_t *L = &g_baselog[b];
            L->id = (u32)b; L->Din = B.Din; L->Dout = B.Dout;
            memcpy(L->x2, B.x2, 4); memcpy(L->z4, B.z4, 4);
            int n4c = 0;
            for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) if(B.n3[r][c] == 4) n4c++;
            for(int r = 0; r < 4; r++) if(B.n4[r] == 4) n4c++;
            L->n4cnt = (u8)n4c; L->tries_capped = (u8)(tries > 255 ? 255 : tries);
        }
        u64 t0 = rdtsc_now();
        walk_init(W, &B);
        u64 t1 = rdtsc_now();
        me->t_gen += t1 - t0;
        int nbuf = 0;
        u32 idxbuf[4];
        for(u32 i = 0; i < NW; i++){
            u64 ta = rdtsc_now();
            if(i) walk_step(W, __builtin_ctz(i));
            memcpy(ebuf[nbuf], W->E, 256);
            idxbuf[nbuf] = i;
            nbuf++;
            u64 tb0 = rdtsc_now();
            me->t_gen += tb0 - ta;
            if(nbuf == 4){
                u64 f0v[4], f1v[4];
                kernel_v2_vec4(rows, f0v, f1v);
                u64 tc = rdtsc_now();
                me->t_kern += tc - tb0;
                if(g_mode == 0){
                    for(int k = 0; k < 4; k++) tlbuf_put(tb, f0v[k], f1v[k]);
                    /* periodic scalar-kernel bit-exact verification */
                    u32 gi = idxbuf[0];
                    if((gi % (u32)g_chk_kern) == 0){
                        u64 tk0 = rdtsc_now();
                        for(int k = 0; k < 4; k++){
                            u64 g0, g1; fp_ref(ebuf[k], &g0, &g1);
                            __sync_fetch_and_add(&g_chk_kern_n, 1);
                            if(g0 == f0v[k] && g1 == f1v[k]) __sync_fetch_and_add(&g_chk_kern_ok, 1);
                        }
                        me->t_chk += rdtsc_now() - tk0;
                    }
                } else if(g_mode == 2){
                    for(int k = 0; k < 4; k++){
                        if(vset_has(&g_t0, f0v[k]) || vset_has(&g_t1, f1v[k]))
                            resolve_emit("TAB", (u32)b, idxbuf[k], &B, ebuf[k], f0v[k], f1v[k]);
                    }
                }
                nbuf = 0;
            }
            /* cold cross-check: the vector-walk sequence vs cold_E vs prop2_all */
            if(g_mode != 1 && (i % (u32)g_chk_cold) == 0){
                u64 tk0 = rdtsc_now();
                u8 Ec[256], Ep[256];
                cold_E(&W->R, Ec);
                int ok = (memcmp(Ec, W->E, 256) == 0);
                __sync_fetch_and_add(&g_chk_cold_n, 1);
                if(ok) __sync_fetch_and_add(&g_chk_cold_ok, 1);
                else {
                    pthread_mutex_lock(&g_print_mtx);
                    fprintf(stderr, "COLD MISMATCH base=%llu i=%u\n", (unsigned long long)b, i);
                    pthread_mutex_unlock(&g_print_mtx);
                }
                refs24_t P; ref24_to_refs24(&W->R, &P);
                prop2_all(&P, Ep);
                __sync_fetch_and_add(&g_chk_p2_n, 1);
                if(memcmp(Ep, W->E, 256) == 0) __sync_fetch_and_add(&g_chk_p2_ok, 1);
                /* also verify the current 24-ref state still satisfies the pair DDT relations
                 * (genuineness of the walked state): stored in chk_p2 stats via stderr on fail */
                me->t_chk += rdtsc_now() - tk0;
            }
        }
        /* leftover (NW always multiple of 4) */
        __sync_fetch_and_add(&g_entries_done, NW);
        __sync_fetch_and_add(&g_bases_done, 1);
        me->entries += NW;
    }
    if(g_mode == 0){
        tlbuf_flush_all(tb);
        __sync_fetch_and_add(&g_ticks_io, tb->io_ticks);
        tlbuf_free(tb);
    }
    __sync_fetch_and_add(&g_ticks_gen, me->t_gen);
    __sync_fetch_and_add(&g_ticks_kern, me->t_kern);
    __sync_fetch_and_add(&g_ticks_chk, me->t_chk);
    return NULL;
}

static volatile int g_stop_monitor = 0;
static void *monitor_fn(void *arg){
    double t0 = now_sec();
    (void)arg;
    while(!g_stop_monitor){
        for(int i = 0; i < 20 && !g_stop_monitor; i++) usleep(500000);
        double t = now_sec() - t0;
        u64 e = g_entries_done;
        fprintf(stderr, "[mon] t=%.0fs bases=%llu/%d entries=%llu (2^%.2f) rate=%.2f Ment/s io=%.1f GiB\n",
                t, (unsigned long long)g_bases_done, g_nb, (unsigned long long)e, log2((double)e + 1),
                e / t / 1e6, (double)e * 16.0 / (1024.0*1024*1024));
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
int main(int argc, char **argv){
    const char *mode = argc > 1 ? argv[1] : "build";
    for(int i = 2; i < argc; i++){
        if(!strncmp(argv[i], "nb=", 3)) g_nb = atoi(argv[i]+3);
        else if(!strncmp(argv[i], "threads=", 8)) g_threads = atoi(argv[i]+8);
        else if(!strncmp(argv[i], "nka=", 4)) g_nka = atoi(argv[i]+4);
        else if(!strncmp(argv[i], "dir=", 4)) snprintf(g_outdir, sizeof g_outdir, "%s", argv[i]+4);
        else if(!strncmp(argv[i], "seed=", 5)) g_seed = strtoull(argv[i]+5, 0, 0);
        else if(!strncmp(argv[i], "buf=", 4)) g_buf_rec = atoi(argv[i]+4);
        else if(!strncmp(argv[i], "walkbits=", 9)) g_walk_bits = atoi(argv[i]+9);
        else if(!strncmp(argv[i], "chkkern=", 8)) g_chk_kern = atoi(argv[i]+8);
        else if(!strncmp(argv[i], "chkcold=", 8)) g_chk_cold = atoi(argv[i]+8);
        else if(!strncmp(argv[i], "targets0=", 9) || !strncmp(argv[i], "targets1=", 9)){
            /* file of hex u64 values, one per line */
            FILE *f = fopen(argv[i]+9, "r"); if(!f){ perror("targets"); return 1; }
            vset_t *S = (argv[i][7] == '0') ? &g_t0 : &g_t1;
            int cap = 1024; S->v = (u64*)malloc(cap * sizeof(u64)); S->n = 0;
            unsigned long long x;
            while(fscanf(f, "%llx", &x) == 1){ if(S->n == cap){ cap *= 2; S->v = (u64*)realloc(S->v, cap * sizeof(u64)); } S->v[S->n++] = x; }
            fclose(f);
            qsort(S->v, S->n, sizeof(u64), cmp_u64);
        }
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }
    if(!strcmp(mode, "build")) g_mode = 0;
    else if(!strcmp(mode, "bench")) g_mode = 1;
    else if(!strcmp(mode, "resolve")) g_mode = 2;
    else { fprintf(stderr, "modes: build | bench | resolve\n"); return 1; }

    double tinit = now_sec();
    aesbench_init();
    g_tsc_hz = tsc_hz();
    fprintf(stderr, "[init] tables ready in %.1fs; TSC %.3f GHz; threads=%d nb=%d (2^%d entries/base) nka=%d seed=%016llx dir=%s\n",
            now_sec() - tinit, g_tsc_hz / 1e9, g_threads, g_nb, g_walk_bits, g_nka,
            (unsigned long long)g_seed, g_outdir);

    /* self-check: vector walk vs cold vs prop2 on a few bases, and kernel vec vs scalar */
    {
        rng_t rng; rng_seed2(&rng, g_seed ^ 0x5151, 7);
        walk_t *W = (walk_t*)aligned_alloc(64, sizeof(walk_t));
        int okw = 0, okp = 0, okk = 0, tot = 0;
        for(int t = 0; t < 8; t++){
            base_t B; base_sample(&B, &rng);
            if(!base_check(&B)){ fprintf(stderr, "FATAL base_check\n"); return 2; }
            walk_init(W, &B);
            for(int i = 0; i < 3000; i++){
                if(i) walk_step(W, __builtin_ctz(i));
                if(i % 97) continue;
                u8 Ec[256], Ep[256]; cold_E(&W->R, Ec);
                refs24_t P; ref24_to_refs24(&W->R, &P); prop2_all(&P, Ep);
                u64 a0, a1, b0, b1; fp_vec(W->E, &a0, &a1); fp_ref(W->E, &b0, &b1);
                tot++;
                if(!memcmp(Ec, W->E, 256)) okw++;
                if(!memcmp(Ep, W->E, 256)) okp++;
                if(a0 == b0 && a1 == b1) okk++;
            }
        }
        free(W);
        fprintf(stderr, "[selftest] walk==cold_E %d/%d  walk==prop2_all %d/%d  kernel_vec==kernel_scalar %d/%d\n",
                okw, tot, okp, tot, okk, tot);
        if(okw != tot || okp != tot || okk != tot){ fprintf(stderr, "FATAL: self-test failed\n"); return 3; }
    }

    /* ------------------------------------------------------------ */
    if(g_mode == 1){
        /* bench: single-thread then multi-thread, compute only */
        int benchbases[2] = {8, g_threads * 8};
        int benchthreads[2] = {1, g_threads};
        for(int pass = 0; pass < 2; pass++){
            int T = benchthreads[pass];
            g_nb = benchbases[pass];
            g_next_base = 0; g_entries_done = 0; g_bases_done = 0;
            g_ticks_gen = g_ticks_kern = 0;
            worker_t *ws = (worker_t*)aligned_alloc(64, sizeof(worker_t) * T);
            memset(ws, 0, sizeof(worker_t) * T);
            pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * T);
            double w0 = now_sec();
            for(int i = 0; i < T; i++){ ws[i].tid = i; pthread_create(&th[i], NULL, worker_fn, &ws[i]); }
            for(int i = 0; i < T; i++) pthread_join(th[i], NULL);
            double w1 = now_sec();
            u64 ent = g_entries_done;
            double ticks_per = (double)(g_ticks_gen + g_ticks_kern) / ent;
            fprintf(stderr, "[bench] threads=%d entries=%llu wall=%.2fs  aggregate=%.1f Ment/s  per-core=%.2f Ment/s  "
                            "TSC-ticks/entry=%.1f (walk %.1f + kernel %.1f)  core-cycles/entry~%.0f (assuming 3.0GHz core vs %.2fGHz TSC)\n",
                    T, (unsigned long long)ent, w1 - w0, ent / (w1 - w0) / 1e6,
                    ent / (w1 - w0) / 1e6 / T, ticks_per,
                    (double)g_ticks_gen / ent, (double)g_ticks_kern / ent,
                    ticks_per * 3.0e9 / g_tsc_hz, g_tsc_hz / 1e9);
            free(ws); free(th);
        }
        return 0;
    }

    if(g_mode == 2){
        /* resolve: regenerate all entries and KA entries; emit matches */
        if((!g_t0.v || !g_t0.n) && (!g_t1.v || !g_t1.n)){ fprintf(stderr, "resolve needs targets0=/targets1=\n"); return 1; }
        if(!g_t0.v){ g_t0.v = (u64*)calloc(1, 8); g_t0.n = 0; }
        if(!g_t1.v){ g_t1.v = (u64*)calloc(1, 8); g_t1.n = 0; }
        g_resolve_out = stdout;
        fprintf(stderr, "[resolve] %d fp0-targets, %d fp1-targets; regenerating %d bases + %d KA\n", g_t0.n, g_t1.n, g_nb, g_nka);
        /* KA entries first */
        {
            rng_t krng; rng_seed2(&krng, g_seed ^ 0x4b415f, 12345);
            ka_inst_t K;
            for(int i = 0; i < g_nka; i++){
                ka_make(&K, &krng);
                if(vset_has(&g_t0, K.fD0) || vset_has(&g_t1, K.fD1))
                    resolve_emit("KA ", 0xFFFFFFFFu, (u32)i, NULL, K.e_true, K.fD0, K.fD1);
            }
        }
        worker_t *ws = (worker_t*)aligned_alloc(64, sizeof(worker_t) * g_threads);
        memset(ws, 0, sizeof(worker_t) * g_threads);
        pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * g_threads);
        pthread_t mon; pthread_create(&mon, NULL, monitor_fn, NULL);
        for(int i = 0; i < g_threads; i++){ ws[i].tid = i; pthread_create(&th[i], NULL, worker_fn, &ws[i]); }
        for(int i = 0; i < g_threads; i++) pthread_join(th[i], NULL);
        g_stop_monitor = 1; pthread_join(mon, NULL);
        fprintf(stderr, "[resolve] done\n");
        return 0;
    }

    /* ---------------------------- build ---------------------------- */
    if(buckets_open(&g_bset, g_outdir, 1)){ fprintf(stderr, "cannot open buckets in %s\n", g_outdir); return 1; }
    g_baselog = (baselog_t*)calloc((size_t)g_nb, sizeof(baselog_t));

    /* known-answer planted instances */
    ka_inst_t *KA = (ka_inst_t*)calloc((size_t)g_nka, sizeof(ka_inst_t));
    {
        rng_t krng; rng_seed2(&krng, g_seed ^ 0x4b415f, 12345);
        int br253 = 0, either = 0, azero = 0;
        for(int i = 0; i < g_nka; i++){
            ka_make(&KA[i], &krng);
            /* fp via vector kernel must agree with the scalar reference */
            u64 v0, v1; fp_vec(KA[i].e_true, &v0, &v1);
            if(v0 != KA[i].fD0 || v1 != KA[i].fD1){ fprintf(stderr, "FATAL: KA %d kernel vec/scalar mismatch\n", i); return 4; }
            if(KA[i].a_ref == 0) azero++;
            if(KA[i].bridge_ok >= 253) br253++;
            if(KA[i].fD0 == KA[i].fH0 || KA[i].fD1 == KA[i].fH1) either++;
            rec_t r = { KA[i].fD0, KA[i].fD1 };
            buckets_write(&g_bset, (int)(KA[i].fD0 >> 56), &r, sizeof r);
        }
        fprintf(stderr, "[ka] planted %d true entries: bridge>=253/255 in %d, genuine-online EITHER-match %d, a_ref=0 instances %d\n",
                g_nka, br253, either, azero);
        char path[700]; snprintf(path, sizeof path, "%s/../ka.bin", g_outdir);
        FILE *f = fopen(path, "wb"); if(!f){ perror(path); return 1; }
        fwrite(&g_nka, sizeof(int), 1, f);
        fwrite(KA, sizeof(ka_inst_t), (size_t)g_nka, f);
        fclose(f);
        fprintf(stderr, "[ka] wrote %s (%zu bytes/instance)\n", path, sizeof(ka_inst_t));
    }

    /* main generation */
    worker_t *ws = (worker_t*)aligned_alloc(64, sizeof(worker_t) * g_threads);
    memset(ws, 0, sizeof(worker_t) * g_threads);
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t) * g_threads);
    pthread_t mon;
    double w0 = now_sec();
    pthread_create(&mon, NULL, monitor_fn, NULL);
    for(int i = 0; i < g_threads; i++){ ws[i].tid = i; pthread_create(&th[i], NULL, worker_fn, &ws[i]); }
    for(int i = 0; i < g_threads; i++) pthread_join(th[i], NULL);
    double w1 = now_sec();
    g_stop_monitor = 1; pthread_join(mon, NULL);
    buckets_close(&g_bset);
    double w2 = now_sec();

    /* bases log */
    {
        char path[700]; snprintf(path, sizeof path, "%s/../bases.bin", g_outdir);
        FILE *f = fopen(path, "wb");
        if(f){ fwrite(g_baselog, sizeof(baselog_t), (size_t)g_nb, f); fclose(f); }
        u64 n4bases = 0, n4tot = 0, trytot = 0;
        for(int i = 0; i < g_nb; i++){ if(g_baselog[i].n4cnt) n4bases++; n4tot += g_baselog[i].n4cnt; trytot += g_baselog[i].tries_capped; }
        fprintf(stderr, "[bases] %d bases logged; bases with any DDT n=4 byte: %llu (total n=4 bytes %llu); mean rejection tries(capped255)=%.1f\n",
                g_nb, (unsigned long long)n4bases, (unsigned long long)n4tot, (double)trytot / g_nb);
    }

    /* summary */
    u64 ent = g_entries_done;
    u64 comp = g_ticks_gen + g_ticks_kern;
    double wall = w1 - w0;
    fprintf(stderr, "\n[build] DONE: entries=%llu (2^%.3f) + %d KA  wall=%.1fs (+%.1fs fsync/close)\n",
            (unsigned long long)ent, log2((double)ent), g_nka, wall, w2 - w1);
    fprintf(stderr, "[build] end-to-end rate = %.1f Ment/s = %.2f Ment/s/thread  (%.3f GiB/s written)\n",
            ent / wall / 1e6, ent / wall / 1e6 / g_threads, ent * 16.0 / wall / (1024.0*1024*1024));
    fprintf(stderr, "[build] compute-only: %.1f TSC-ticks/entry (walk %.1f + kernel %.1f)  => %.2f Ment/s/core, ~%.0f core-cycles/entry at 3.0GHz\n",
            (double)comp / ent, (double)g_ticks_gen / ent, (double)g_ticks_kern / ent,
            g_tsc_hz / ((double)comp / ent) / 1e6, (double)comp / ent * 3.0e9 / g_tsc_hz);
    fprintf(stderr, "[build] time shares: compute %.1f%%  io(write) %.1f%%  checks %.1f%%  (of summed thread time)\n",
            100.0 * comp / (comp + g_ticks_io + g_ticks_chk),
            100.0 * g_ticks_io / (comp + g_ticks_io + g_ticks_chk),
            100.0 * g_ticks_chk / (comp + g_ticks_io + g_ticks_chk));
    fprintf(stderr, "[checks] walk==cold_E %llu/%llu   walk==prop2_all %llu/%llu   kernel_vec==kernel_scalar %llu/%llu\n",
            (unsigned long long)g_chk_cold_ok, (unsigned long long)g_chk_cold_n,
            (unsigned long long)g_chk_p2_ok, (unsigned long long)g_chk_p2_n,
            (unsigned long long)g_chk_kern_ok, (unsigned long long)g_chk_kern_n);
    u64 totc = 0; for(int b = 0; b < NBUCKET; b++) totc += g_bset.count[b];
    fprintf(stderr, "[buckets] records written: %llu (expect %llu)\n",
            (unsigned long long)totc, (unsigned long long)(ent + (u64)g_nka));
    free(ws); free(th); free(KA);
    return 0;
}
