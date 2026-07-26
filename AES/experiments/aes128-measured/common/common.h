// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#ifndef AESBENCH_COMMON_H
#define AESBENCH_COMMON_H
/* common.h -- shared machinery for the real-scale offline table build.
 *
 * GENERATOR GEOMETRY (DDT-Gray, as in AES/c/ddt-gray/ddt_gray_e2e.c):
 *   Parameter base (10 bytes, DFJ'13 Prop-2 style):
 *        Din = Delta z_1[0]   x2[0..3] = x_2[col0]   Dout = Delta w_4[0]   z4[0..3] = z_4[col0]
 *   rebound -> Delta x_3[4][4], Delta y_3[4][4], Delta x_4[diag], Delta z_4[col0] and a base
 *   24-ref state (x2[4], x3[4][4], x4[diag]).  The 20 DDT branch bits (4 for
 *   x_4[diag], 16 for x_3) are walked in binary-reflected Gray order; every
 *   one of the 2^20 states is a valid 24-ref state and yields a genuine
 *   255-element difference sequence  d_w = Delta x_5[0]  (w = 1..255).
 *   This is the SAME 3-round map as prop2.h's refs24_t (x3c,x4,x5d) ->
 *   Delta x_6[0] one round later; cross-checked numerically at runtime.
 *
 * FINGERPRINT: kernel_v2 (brute-free):
 *   input  e[256] (= the d_w sequence; the kernel inverts internally),
 *   output (fp0, fp1) = raw-parity and Add-0-parity chi* canonical FNV64s.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>
#include <x86intrin.h>
#include <immintrin.h>

#include "aes_core.h"
#include "prop2.h"
#include "chi256.h"
#include "chi256_fast.h"
#include "kernel_fast.h"
#include "kernel_v2_scalar.h"
#include "kernel_v2_vec.h"

/* ------------------------------------------------------------------ */
/* RNG: splitmix64 / xoroshiro-ish, per-thread state                     */
/* ------------------------------------------------------------------ */
typedef struct { u64 s0, s1; } rng_t;
static inline u64 splitmix64(u64 *x){
    u64 z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline void rng_seed2(rng_t *r, u64 a, u64 b){
    u64 x = a * 0x9e3779b97f4a7c15ULL + b * 0xc2b2ae3d27d4eb4fULL + 0x165667b19e3779f9ULL;
    r->s0 = splitmix64(&x); r->s1 = splitmix64(&x);
    if(!r->s0 && !r->s1) r->s0 = 1;
}
static inline u64 rng_next(rng_t *r){
    u64 s0 = r->s0, s1 = r->s1, res = s0 + s1;
    s1 ^= s0;
    r->s0 = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
    r->s1 = (s1 << 36) | (s1 >> 28);
    return res;
}
static inline u8 rng8(rng_t *r){ return (u8)(rng_next(r) >> 33); }
static inline u8 rng8nz(rng_t *r){ u8 v; do { v = rng8(r); } while(!v); return v; }

static inline double now_sec(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static inline u64 rdtsc_now(void){ unsigned aux; return __rdtscp(&aux); }

/* TSC frequency calibration (ticks per second) */
static double tsc_hz(void){
    double t0 = now_sec(); u64 c0 = rdtsc_now();
    struct timespec req = {0, 300 * 1000 * 1000};
    nanosleep(&req, NULL);
    double t1 = now_sec(); u64 c1 = rdtsc_now();
    return (double)(c1 - c0) / (t1 - t0);
}

/* ------------------------------------------------------------------ */
/* Global init                                                           */
/* ------------------------------------------------------------------ */
static u8 ZROWS[256] __attribute__((aligned(64)));
static __m512i VSBOX_MAT;   /* AES affine matrix for GF2P8AFFINEINVQB */

static void aesbench_init(void){
    aes_core_init(); ddt_init(); chifast_init();
    kernel_v2_vec_init();           /* kernel_init + canon_even_init + vec tables */
    memset(ZROWS, 0, sizeof ZROWS);
    /* Build the AES affine matrix in GFNI encoding and verify against SBOX.
     * GFNI: out.bit[i] = parity(A.byte[7-i] & in) ^ imm.bit[i], where
     * A.byte[k] = bits 8k..8k+7 of the qword operand.
     * AES: out_i = in_i ^ in_{i+4} ^ in_{i+5} ^ in_{i+6} ^ in_{i+7} (mod 8) ^ c_i, c=0x63. */
    u64 M = 0;
    for(int i = 0; i < 8; i++){
        u8 m = (u8)((1u << i) | (1u << ((i+4)&7)) | (1u << ((i+5)&7)) | (1u << ((i+6)&7)) | (1u << ((i+7)&7)));
        M |= (u64)m << (8 * (7 - i));
    }
    VSBOX_MAT = _mm512_set1_epi64((long long)M);
    /* self-test */
    u8 in[256] __attribute__((aligned(64))), out[256] __attribute__((aligned(64)));
    for(int i = 0; i < 256; i++) in[i] = (u8)i;
    for(int i = 0; i < 4; i++){
        __m512i v = _mm512_load_si512(in + 64*i);
        __m512i s = _mm512_gf2p8affineinv_epi64_epi8(v, VSBOX_MAT, 0x63);
        _mm512_store_si512(out + 64*i, s);
    }
    int bad = 0;
    for(int i = 0; i < 256; i++) if(out[i] != SBOX[i]) bad++;
    if(bad){
        /* try the widely-quoted constant as a fallback */
        VSBOX_MAT = _mm512_set1_epi64((long long)0xF1E3C78F1F3E7CF8ULL);
        bad = 0;
        for(int i = 0; i < 4; i++){
            __m512i v = _mm512_load_si512(in + 64*i);
            __m512i s = _mm512_gf2p8affineinv_epi64_epi8(v, VSBOX_MAT, 0x63);
            _mm512_store_si512(out + 64*i, s);
        }
        for(int i = 0; i < 256; i++) if(out[i] != SBOX[i]) bad++;
        if(bad){ fprintf(stderr, "FATAL: vector SBOX self-test failed (%d mismatches)\n", bad); abort(); }
    }
}
static inline __m512i vsbox(__m512i x){ return _mm512_gf2p8affineinv_epi64_epi8(x, VSBOX_MAT, 0x63); }
static inline __m512i vgmul(__m512i x, u8 c){ return _mm512_gf2p8mul_epi8(x, _mm512_set1_epi8((char)c)); }

/* ------------------------------------------------------------------ */
/* Parameter base, rebound (ddt_gray geometry, aes_core tables)         */
/* ------------------------------------------------------------------ */
typedef struct {
    u8 Din, Dout;
    u8 x2[4], z4[4];
    /* derived pair differences */
    u8 Dx3[4][4], Dy3[4][4], Dx4[4], Dz4[4];
    /* base 24-ref state */
    u8 bx2[4], bx3[4][4], bx4[4];
    /* DDT multiplicities (2 or 4) */
    u8 n3[4][4], n4[4];
} base_t;

/* Fills derived fields; returns 1 if valid (all DDT entries nonzero). */
static int base_rebound(base_t *B){
    u8 Dy2[4];
    for(int r = 0; r < 4; r++){
        u8 Dx2 = gmul(MCc[r][0], B->Din);
        Dy2[r] = (u8)(SBOX[B->x2[r] ^ Dx2] ^ SBOX[B->x2[r]]);
        if(Dy2[r] == 0) return 0;
    }
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){
        int rc = (4 - c) & 3;
        B->Dx3[r][c] = gmul(MCc[r][rc], Dy2[rc]);
    }
    for(int r = 0; r < 4; r++){
        B->Dz4[r] = gmul(iMCc[r][0], B->Dout);
        u8 x4 = iSBOX[B->z4[r]];
        B->Dx4[r] = (u8)(iSBOX[B->z4[r] ^ B->Dz4[r]] ^ x4);
        B->bx4[r] = x4;
        if(B->Dx4[r] == 0) return 0;
    }
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){
        int cc = ((c - r) + 4) & 3;
        B->Dy3[r][c] = gmul(iMCc[r][cc], B->Dx4[cc]);
    }
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){
        u8 n = DDTn[B->Dx3[r][c]][B->Dy3[r][c]];
        if(n == 0) return 0;
        B->n3[r][c] = n;
        B->bx3[r][c] = DDTsol[B->Dx3[r][c]][B->Dy3[r][c]][0];
    }
    for(int r = 0; r < 4; r++){
        /* x_4[diag] pair diff check: SB(x4)^SB(x4^Dx4) must equal Dz4 (it does by construction) */
        B->n4[r] = DDTn[B->Dx4[r]][B->Dz4[r]];
        if(B->n4[r] == 0) return 0; /* cannot happen */
    }
    for(int r = 0; r < 4; r++) B->bx2[r] = B->x2[r];
    return 1;
}
/* verify base state satisfies the DDT relations */
static int base_check(const base_t *B){
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++)
        if((u8)(SBOX[B->bx3[r][c]] ^ SBOX[B->bx3[r][c] ^ B->Dx3[r][c]]) != B->Dy3[r][c]) return 0;
    for(int r = 0; r < 4; r++)
        if((u8)(SBOX[B->bx4[r]] ^ SBOX[B->bx4[r] ^ B->Dx4[r]]) != B->Dz4[r]) return 0;
    return 1;
}
/* rejection-sample a valid base from rng; returns number of tries */
static int base_sample(base_t *B, rng_t *R){
    int tries = 0;
    for(;;){
        tries++;
        B->Din = rng8nz(R); B->Dout = rng8nz(R);
        for(int r = 0; r < 4; r++){ B->x2[r] = rng8(R); B->z4[r] = rng8(R); }
        if(base_rebound(B)) return tries;
    }
}

/* ------------------------------------------------------------------ */
/* 24-ref state + scalar reference sequence (cold_E, from ddt_gray)      */
/* ------------------------------------------------------------------ */
typedef struct { u8 x2[4]; u8 x3[4][4]; u8 x4[4]; } ref24_t;

static void cold_E(const ref24_t *R, u8 E[256]){
    u8 sb_x2[4], sb_x3[4][4], sb_x4[4];
    for(int r = 0; r < 4; r++) sb_x2[r] = SBOX[R->x2[r]];
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) sb_x3[r][c] = SBOX[R->x3[r][c]];
    for(int r = 0; r < 4; r++) sb_x4[r] = SBOX[R->x4[r]];
    E[0] = 0;
    for(int w = 1; w < 256; w++){
        u8 dy2[4];
        for(int r = 0; r < 4; r++) dy2[r] = (u8)(SBOX[R->x2[r] ^ gmul(MCc[r][0], (u8)w)] ^ sb_x2[r]);
        u8 dy3[4][4];
        for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){
            int rc = (4 - c) & 3;
            u8 dx3 = gmul(MCc[r][rc], dy2[rc]);
            dy3[r][c] = (u8)(SBOX[R->x3[r][c] ^ dx3] ^ sb_x3[r][c]);
        }
        u8 dy4[4];
        for(int r = 0; r < 4; r++){
            u8 dx4 = 0;
            for(int rp = 0; rp < 4; rp++) dx4 ^= gmul(MCc[r][rp], dy3[rp][(r + rp) & 3]);
            dy4[r] = (u8)(SBOX[R->x4[r] ^ dx4] ^ sb_x4[r]);
        }
        E[w] = (u8)(gmul(2, dy4[0]) ^ gmul(3, dy4[1]) ^ dy4[2] ^ dy4[3]);
    }
}
/* map to prop2.h refs24_t convention (the next-round labelling) */
static void ref24_to_refs24(const ref24_t *R, refs24_t *P){
    for(int i = 0; i < 4; i++) P->x3c[i] = R->x2[i];
    for(int c = 0; c < 4; c++) for(int r = 0; r < 4; r++) P->x4[4*c + r] = R->x3[r][c];
    for(int i = 0; i < 4; i++) P->x5d[i] = R->x4[i];
}

/* ------------------------------------------------------------------ */
/* Vectorized incremental Gray-walk state (per thread)                   */
/* ------------------------------------------------------------------ */
typedef struct {
    /* per-omega caches, omega = 0..255 (omega=0 lane is all-zero diffs) */
    u8 dx3a[4][4][256];   /* static per base (depends on x2 only) */
    u8 dy3a[4][4][256];
    u8 dx4a[4][256];
    u8 dy4a[4][256];
    u8 E[256];            /* the current difference sequence, E[0]=0 */
    ref24_t R;            /* scalar refs (for cross-checks / provenance) */
    u8 sbx3[4][4], sbx4[4];
    const base_t *B;
} __attribute__((aligned(64))) walk_t;

/* cache init from base: scalar, once per base (cost ~2*10^4 SB) */
static void walk_init(walk_t *W, const base_t *B){
    W->B = B;
    for(int r = 0; r < 4; r++) W->R.x2[r] = B->bx2[r];
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) W->R.x3[r][c] = B->bx3[r][c];
    for(int r = 0; r < 4; r++) W->R.x4[r] = B->bx4[r];
    u8 sbx2[4];
    for(int r = 0; r < 4; r++) sbx2[r] = SBOX[W->R.x2[r]];
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++) W->sbx3[r][c] = SBOX[W->R.x3[r][c]];
    for(int r = 0; r < 4; r++) W->sbx4[r] = SBOX[W->R.x4[r]];
    /* omega = 0 lane: all zeros */
    for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){ W->dx3a[r][c][0] = 0; W->dy3a[r][c][0] = 0; }
    for(int r = 0; r < 4; r++){ W->dx4a[r][0] = 0; W->dy4a[r][0] = 0; }
    W->E[0] = 0;
    for(int w = 1; w < 256; w++){
        u8 dy2[4];
        for(int r = 0; r < 4; r++) dy2[r] = (u8)(SBOX[W->R.x2[r] ^ gmul(MCc[r][0], (u8)w)] ^ sbx2[r]);
        for(int r = 0; r < 4; r++) for(int c = 0; c < 4; c++){
            int rc = (4 - c) & 3;
            u8 dx3 = gmul(MCc[r][rc], dy2[rc]);
            W->dx3a[r][c][w] = dx3;
            W->dy3a[r][c][w] = (u8)(SBOX[W->R.x3[r][c] ^ dx3] ^ W->sbx3[r][c]);
        }
        u8 e = 0;
        for(int r = 0; r < 4; r++){
            u8 dx4 = 0;
            for(int rp = 0; rp < 4; rp++) dx4 ^= gmul(MCc[r][rp], W->dy3a[rp][(r + rp) & 3][w]);
            W->dx4a[r][w] = dx4;
            u8 dy4 = (u8)(SBOX[W->R.x4[r] ^ dx4] ^ W->sbx4[r]);
            W->dy4a[r][w] = dy4;
            e ^= gmul(MCc[0][r], dy4);
        }
        W->E[w] = e;
    }
}

/* flip x_4[diag r] by Dx4[r]  (Gray bit r in 0..3) */
static inline void walk_flip_z4(walk_t *W, int r){
    W->R.x4[r] ^= W->B->Dx4[r];
    W->sbx4[r] = SBOX[W->R.x4[r]];
    const __m512i x4b = _mm512_set1_epi8((char)W->R.x4[r]);
    const __m512i sbb = _mm512_set1_epi8((char)W->sbx4[r]);
    const __m512i mc  = _mm512_set1_epi8((char)MCc[0][r]);
    u8 *dx4 = W->dx4a[r], *dy4 = W->dy4a[r], *E = W->E;
    for(int i = 0; i < 4; i++){
        __m512i vdx4 = _mm512_load_si512(dx4 + 64*i);
        __m512i nd   = _mm512_xor_si512(vsbox(_mm512_xor_si512(x4b, vdx4)), sbb);
        __m512i ody  = _mm512_load_si512(dy4 + 64*i);
        __m512i de   = _mm512_gf2p8mul_epi8(_mm512_xor_si512(nd, ody), mc);
        __m512i ve   = _mm512_load_si512(E + 64*i);
        _mm512_store_si512(E + 64*i, _mm512_xor_si512(ve, de));
        _mm512_store_si512(dy4 + 64*i, nd);
    }
}
/* flip x_3[rp][cp] by Dx3[rp][cp]  (Gray bit 4+p, p = rp + 4*cp) */
static inline void walk_flip_x3(walk_t *W, int rp, int cp){
    W->R.x3[rp][cp] ^= W->B->Dx3[rp][cp];
    W->sbx3[rp][cp] = SBOX[W->R.x3[rp][cp]];
    int rs = ((cp - rp) + 4) & 3;                 /* affected x_4 diag index */
    const __m512i x3b = _mm512_set1_epi8((char)W->R.x3[rp][cp]);
    const __m512i sb3 = _mm512_set1_epi8((char)W->sbx3[rp][cp]);
    const __m512i mc3 = _mm512_set1_epi8((char)MCc[rs][rp]);
    const __m512i x4b = _mm512_set1_epi8((char)W->R.x4[rs]);
    const __m512i sb4 = _mm512_set1_epi8((char)W->sbx4[rs]);
    const __m512i mc4 = _mm512_set1_epi8((char)MCc[0][rs]);
    u8 *dx3 = W->dx3a[rp][cp], *dy3 = W->dy3a[rp][cp];
    u8 *dx4 = W->dx4a[rs], *dy4 = W->dy4a[rs], *E = W->E;
    for(int i = 0; i < 4; i++){
        __m512i vdx3 = _mm512_load_si512(dx3 + 64*i);
        __m512i nd3  = _mm512_xor_si512(vsbox(_mm512_xor_si512(x3b, vdx3)), sb3);
        __m512i ody3 = _mm512_load_si512(dy3 + 64*i);
        __m512i vdx4 = _mm512_load_si512(dx4 + 64*i);
        vdx4 = _mm512_xor_si512(vdx4, _mm512_gf2p8mul_epi8(_mm512_xor_si512(nd3, ody3), mc3));
        _mm512_store_si512(dy3 + 64*i, nd3);
        _mm512_store_si512(dx4 + 64*i, vdx4);
        __m512i nd4  = _mm512_xor_si512(vsbox(_mm512_xor_si512(x4b, vdx4)), sb4);
        __m512i ody4 = _mm512_load_si512(dy4 + 64*i);
        __m512i ve   = _mm512_load_si512(E + 64*i);
        ve = _mm512_xor_si512(ve, _mm512_gf2p8mul_epi8(_mm512_xor_si512(nd4, ody4), mc4));
        _mm512_store_si512(E + 64*i, ve);
        _mm512_store_si512(dy4 + 64*i, nd4);
    }
}
/* Gray step i (i>=1): bit j = ctz(i) */
static inline void walk_step(walk_t *W, int j){
    if(j < 4) walk_flip_z4(W, j);
    else { int p = j - 4; walk_flip_x3(W, p & 3, p >> 2); }
}

/* ------------------------------------------------------------------ */
/* Fingerprints                                                          */
/* ------------------------------------------------------------------ */
static inline void fp_vec(const u8 *E, u64 *f0, u64 *f1){
    kernel_v2_vec(E, ZROWS, ZROWS, ZROWS, f0, f1);
}
static inline void fp_ref(const u8 *E, u64 *f0, u64 *f1){
    kernel_full_v2(E, ZROWS, ZROWS, ZROWS, f0, f1);
}
/* FNV64 of the ordered sequence E[1..255] */
static inline u64 seq_hash(const u8 *E){
    u64 h = 0xcbf29ce484222325ULL;
    for(int i = 1; i < 256; i++){ h ^= E[i]; h *= 0x100000001b3ULL; }
    return h;
}
/* FNV64 of the multiset (sorted counts) of d = inv(E[1..255]) */
static inline u64 multiset_hash(const u8 *E){
    u8 cnt[256]; memset(cnt, 0, 256);
    for(int i = 1; i < 256; i++) cnt[GF_inv[E[i]]]++;
    u64 h = 0xcbf29ce484222325ULL;
    for(int v = 0; v < 256; v++){ h ^= cnt[v]; h *= 0x100000001b3ULL; }
    return h;
}
/* parity set chi(D) of d = inv(E) as 4 qwords; returns FNV of it */
static inline u64 chi_of_seq(const u8 *E, chi_t *out){
    chi_t c; chi_clr(&c);
    for(int i = 1; i < 256; i++) chi_flip(&c, GF_inv[E[i]]);
    if(out) *out = c;
    return chi_hash(&c);
}

/* ------------------------------------------------------------------ */
/* Record format / bucket writer                                         */
/* ------------------------------------------------------------------ */
typedef struct { u64 fp0, fp1; } rec_t;   /* 16 bytes, little-endian */

#define NBUCKET 256
typedef struct {
    int fd[NBUCKET];
    pthread_mutex_t mtx[NBUCKET];
    u64 count[NBUCKET];
    char dir[512];
} bucketset_t;

static int buckets_open(bucketset_t *S, const char *dir, int truncate){
    snprintf(S->dir, sizeof S->dir, "%s", dir);
    char path[640];
    for(int b = 0; b < NBUCKET; b++){
        snprintf(path, sizeof path, "%s/b%02x.bin", dir, b);
        S->fd[b] = open(path, O_WRONLY | O_CREAT | (truncate ? O_TRUNC : O_APPEND), 0644);
        if(S->fd[b] < 0){ perror(path); return -1; }
        pthread_mutex_init(&S->mtx[b], NULL);
        S->count[b] = 0;
    }
    return 0;
}
static void buckets_close(bucketset_t *S){
    for(int b = 0; b < NBUCKET; b++){ fsync(S->fd[b]); close(S->fd[b]); }
}
static void buckets_write(bucketset_t *S, int b, const void *buf, size_t len){
    pthread_mutex_lock(&S->mtx[b]);
    const char *p = (const char*)buf; size_t left = len;
    while(left){
        ssize_t w = write(S->fd[b], p, left);
        if(w <= 0){ perror("write"); abort(); }
        p += w; left -= (size_t)w;
    }
    S->count[b] += len / sizeof(rec_t);
    pthread_mutex_unlock(&S->mtx[b]);
}

/* thread-local record buffering */
typedef struct {
    rec_t *buf[NBUCKET];
    int    n[NBUCKET];
    int    cap;
    bucketset_t *S;
    u64 io_ticks;
} tlbuf_t;
static void tlbuf_init(tlbuf_t *T, bucketset_t *S, int cap_records){
    T->S = S; T->cap = cap_records; T->io_ticks = 0;
    for(int b = 0; b < NBUCKET; b++){
        T->buf[b] = (rec_t*)malloc(sizeof(rec_t) * cap_records);
        T->n[b] = 0;
    }
}
static inline void tlbuf_flush_one(tlbuf_t *T, int b){
    if(!T->n[b]) return;
    u64 t0 = rdtsc_now();
    if(T->S) buckets_write(T->S, b, T->buf[b], (size_t)T->n[b] * sizeof(rec_t));
    T->io_ticks += rdtsc_now() - t0;
    T->n[b] = 0;
}
static inline void tlbuf_put(tlbuf_t *T, u64 fp0, u64 fp1){
    int b = (int)(fp0 >> 56);
    T->buf[b][T->n[b]].fp0 = fp0;
    T->buf[b][T->n[b]].fp1 = fp1;
    if(++T->n[b] == T->cap) tlbuf_flush_one(T, b);
}
static void tlbuf_flush_all(tlbuf_t *T){
    for(int b = 0; b < NBUCKET; b++) tlbuf_flush_one(T, b);
}
static void tlbuf_free(tlbuf_t *T){
    for(int b = 0; b < NBUCKET; b++) free(T->buf[b]);
}

/* ------------------------------------------------------------------ */
/* Known-answer (planted) instance: real AES-128 7-round, a1_e2e geometry  */
/*   refs at rounds 3,4,5 (x_3[col0], x_4, x_5[diag]) -> e = Delta x_6[0]  */
/* ------------------------------------------------------------------ */
#define NR 7
typedef struct {
    u8 key[16], P[16];
    u8 x3c[4], x4[16], x5d[4];
    u8 a_ref;        /* x_6[0,0] */
    u8 km1[4];       /* k_0[diag0]   (true) */
    u8 k7a[4];       /* k_7[anti-diag] (true) */
    u8 e_true[256];
    u64 fD0, fD1;    /* offline fingerprint of the true entry */
    u64 fH0, fH1;    /* genuine online fingerprint with the true key bytes */
    int bridge_ok;   /* bridge identity count / 255 */
} ka_inst_t;

static void ka_make(ka_inst_t *K, rng_t *rng){
    aesctx ctx;
    for(int i = 0; i < 16; i++) K->key[i] = rng8(rng);
    aes_set_key(&ctx, K->key, NR);
    for(int i = 0; i < 16; i++) K->P[i] = rng8(rng);
    u8 C[16]; trace_t T;
    aes_enc_trace(&ctx, K->P, C, &T);
    refs24_t R;
    for(int i = 0; i < 4; i++) R.x3c[i] = T.x[3][i];
    for(int i = 0; i < 16; i++) R.x4[i] = T.x[4][i];
    for(int i = 0; i < 4; i++) R.x5d[i] = T.x[5][DIAG0[i]];
    memcpy(K->x3c, R.x3c, 4); memcpy(K->x4, R.x4, 16); memcpy(K->x5d, R.x5d, 4);
    K->a_ref = T.x[6][0];
    for(int i = 0; i < 4; i++) K->km1[i] = ctx.rk[0][DIAG0[i]];
    for(int i = 0; i < 4; i++) K->k7a[i] = ctx.rk[7][ADIAG[i]];
    prop2_all(&R, K->e_true);
    fp_ref(K->e_true, &K->fD0, &K->fD1);
    /* genuine online computation with TRUE km1 (delta-set) and TRUE k7a (peel) */
    /* w_1[col0] of the reference, rebuilt from P and km1 exactly as the attacker does */
    u8 z1r[4]; for(int i = 0; i < 4; i++) z1r[i] = SBOX[K->P[DIAG0[i]] ^ K->km1[i]];
    u8 w1r[4];
    for(int r = 0; r < 4; r++) w1r[r] = (u8)(gmul(MCc[r][0], z1r[0]) ^ gmul(MCc[r][1], z1r[1]) ^ gmul(MCc[r][2], z1r[2]) ^ gmul(MCc[r][3], z1r[3]));
    u8 Cd[256][16];
    for(int dv = 0; dv < 256; dv++){
        u8 w1[4] = { (u8)(w1r[0] ^ dv), w1r[1], w1r[2], w1r[3] };
        u8 z1[4];
        for(int r = 0; r < 4; r++) z1[r] = (u8)(gmul(iMCc[r][0], w1[0]) ^ gmul(iMCc[r][1], w1[1]) ^ gmul(iMCc[r][2], w1[2]) ^ gmul(iMCc[r][3], w1[3]));
        u8 PP[16]; memcpy(PP, K->P, 16);
        for(int i = 0; i < 4; i++) PP[DIAG0[i]] = (u8)(iSBOX[z1[i]] ^ K->km1[i]);
        aes_enc(&ctx, PP, Cd[dv]);
    }
    /* peel with true k7a */
    u8 x7ref[4]; for(int i = 0; i < 4; i++) x7ref[i] = iSBOX[Cd[0][ADIAG[i]] ^ K->k7a[i]];
    u8 eH[256] __attribute__((aligned(64))); eH[0] = 0;
    for(int dv = 1; dv < 256; dv++){
        u8 dy6 = 0;
        for(int i = 0; i < 4; i++) dy6 ^= gmul(iMCc[0][i], (u8)(iSBOX[Cd[dv][ADIAG[i]] ^ K->k7a[i]] ^ x7ref[i]));
        eH[dv] = Linv[dy6];
    }
    fp_ref(eH, &K->fH0, &K->fH1);
    /* bridge identity check h = a^2 d + a */
    u8 a = K->a_ref; int brok = 0;
    u8 x2ref0 = T.x[2][0];
    for(int dv = 1; dv < 256; dv++){
        u8 h  = GF_inv[eH[dv]];
        u8 om = (u8)(SBOX[x2ref0 ^ dv] ^ SBOX[x2ref0]);
        u8 e  = K->e_true[om], d = GF_inv[e];
        u8 hp = (u8)(gmul(gmul(a, a), d) ^ a);
        if(e == 0) hp = 0;
        if(e == a) hp = a;
        if(h == hp) brok++;
    }
    K->bridge_ok = brok;
}

#endif
