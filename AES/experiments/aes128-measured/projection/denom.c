// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* denom.c -- same-machine denominators: cycles per AES-128 encryption for
 *   (1) a software T-table implementation (classic 4x1KB Te tables, 10 rounds,
 *       last round via the S-box), and
 *   (2) AES-NI, pipelined throughput (8 independent streams) and 1-stream latency.
 * Both verified against the byte-wise reference aes_core.h (10 rounds) on
 * random blocks, then timed over 10^8 blocks with rdtscp + imul core-clock
 * calibration.  Counter-mode style independent blocks => throughput.
 *
 * Build: gcc -O3 -march=native -I../offline -o denom denom.c
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
#include <wmmintrin.h>
#include "aes_core.h"

static inline u64 rdt(void){ unsigned a; return __rdtscp(&a); }
static double calib_ratio(void){
    u64 best=~0ULL;
    for(int rep=0;rep<200;rep++){
        u64 y=3,t0=rdt();
        for(long i=0;i<20000;i++){
            __asm__ volatile("imulq %1,%0\n\timulq %1,%0\n\timulq %1,%0\n\timulq %1,%0\n\timulq %1,%0":"+r"(y):"r"(5UL));
        }
        u64 t1=rdt();
        __asm__ volatile(""::"r"(y));
        if(t1-t0<best)best=t1-t0;
    }
    return (3.0*100000.0)/(double)best;
}
static double tsc_hz(void){
    struct timespec a,b; clock_gettime(CLOCK_MONOTONIC,&a); u64 c0=rdt();
    struct timespec req={0,300*1000*1000}; nanosleep(&req,NULL);
    clock_gettime(CLOCK_MONOTONIC,&b); u64 c1=rdt();
    double dt=(b.tv_sec-a.tv_sec)+1e-9*(b.tv_nsec-a.tv_nsec);
    return (c1-c0)/dt;
}

/* ---------------- T-table AES-128 (encryption) ---------------- */
static u32 Te0[256], Te1[256], Te2[256], Te3[256];
static u32 RK[44];
static inline u32 rotr32(u32 x, int n){ return (x >> n) | (x << (32 - n)); }
static void ttab_init(void){
    for(int x = 0; x < 256; x++){
        u8 s = SBOX[x];
        u8 s2 = gmul(s, 2), s3 = gmul(s, 3);
        /* big-endian column word: (2s,1s,1s,3s) */
        u32 w = ((u32)s2 << 24) | ((u32)s << 16) | ((u32)s << 8) | (u32)s3;
        Te0[x] = w; Te1[x] = rotr32(w, 8); Te2[x] = rotr32(w, 16); Te3[x] = rotr32(w, 24);
    }
}
static void ttab_keysched(const u8 key[16]){
    /* standard FIPS-197 word key expansion, big-endian words */
    static const u8 Rcon[11]={0,1,2,4,8,16,32,64,128,0x1b,0x36};
    for(int i = 0; i < 4; i++)
        RK[i] = ((u32)key[4*i]<<24)|((u32)key[4*i+1]<<16)|((u32)key[4*i+2]<<8)|key[4*i+3];
    for(int i = 4; i < 44; i++){
        u32 t = RK[i-1];
        if(i % 4 == 0){
            t = (t << 8) | (t >> 24);
            t = ((u32)SBOX[(t>>24)&0xff]<<24)|((u32)SBOX[(t>>16)&0xff]<<16)|((u32)SBOX[(t>>8)&0xff]<<8)|SBOX[t&0xff];
            t ^= (u32)Rcon[i/4] << 24;
        }
        RK[i] = RK[i-4] ^ t;
    }
}
static inline void ttab_enc(const u8 in[16], u8 out[16]){
    u32 s0 = (((u32)in[0]<<24)|((u32)in[1]<<16)|((u32)in[2]<<8)|in[3]) ^ RK[0];
    u32 s1 = (((u32)in[4]<<24)|((u32)in[5]<<16)|((u32)in[6]<<8)|in[7]) ^ RK[1];
    u32 s2 = (((u32)in[8]<<24)|((u32)in[9]<<16)|((u32)in[10]<<8)|in[11]) ^ RK[2];
    u32 s3 = (((u32)in[12]<<24)|((u32)in[13]<<16)|((u32)in[14]<<8)|in[15]) ^ RK[3];
    u32 t0, t1, t2, t3;
    for(int r = 1; r < 10; r++){
        t0 = Te0[s0>>24] ^ Te1[(s1>>16)&0xff] ^ Te2[(s2>>8)&0xff] ^ Te3[s3&0xff] ^ RK[4*r];
        t1 = Te0[s1>>24] ^ Te1[(s2>>16)&0xff] ^ Te2[(s3>>8)&0xff] ^ Te3[s0&0xff] ^ RK[4*r+1];
        t2 = Te0[s2>>24] ^ Te1[(s3>>16)&0xff] ^ Te2[(s0>>8)&0xff] ^ Te3[s1&0xff] ^ RK[4*r+2];
        t3 = Te0[s3>>24] ^ Te1[(s0>>16)&0xff] ^ Te2[(s1>>8)&0xff] ^ Te3[s2&0xff] ^ RK[4*r+3];
        s0 = t0; s1 = t1; s2 = t2; s3 = t3;
    }
    /* final round: SubBytes + ShiftRows + AddRoundKey (no MixColumns) */
    t0 = ((u32)SBOX[s0>>24]<<24) ^ ((u32)SBOX[(s1>>16)&0xff]<<16) ^ ((u32)SBOX[(s2>>8)&0xff]<<8) ^ SBOX[s3&0xff] ^ RK[40];
    t1 = ((u32)SBOX[s1>>24]<<24) ^ ((u32)SBOX[(s2>>16)&0xff]<<16) ^ ((u32)SBOX[(s3>>8)&0xff]<<8) ^ SBOX[s0&0xff] ^ RK[41];
    t2 = ((u32)SBOX[s2>>24]<<24) ^ ((u32)SBOX[(s3>>16)&0xff]<<16) ^ ((u32)SBOX[(s0>>8)&0xff]<<8) ^ SBOX[s1&0xff] ^ RK[42];
    t3 = ((u32)SBOX[s3>>24]<<24) ^ ((u32)SBOX[(s0>>16)&0xff]<<16) ^ ((u32)SBOX[(s1>>8)&0xff]<<8) ^ SBOX[s2&0xff] ^ RK[43];
    out[0]=t0>>24; out[1]=t0>>16; out[2]=t0>>8; out[3]=(u8)t0;
    out[4]=t1>>24; out[5]=t1>>16; out[6]=t1>>8; out[7]=(u8)t1;
    out[8]=t2>>24; out[9]=t2>>16; out[10]=t2>>8; out[11]=(u8)t2;
    out[12]=t3>>24; out[13]=t3>>16; out[14]=t3>>8; out[15]=(u8)t3;
}

/* ---------------- AES-NI ---------------- */
static __m128i NIK[11];
static void aesni_keysched(const aesctx *ctx){
    for(int r = 0; r <= 10; r++) NIK[r] = _mm_loadu_si128((const __m128i*)ctx->rk[r]);
}
static inline __m128i aesni_enc1(__m128i s){
    s = _mm_xor_si128(s, NIK[0]);
    for(int r = 1; r < 10; r++) s = _mm_aesenc_si128(s, NIK[r]);
    return _mm_aesenclast_si128(s, NIK[10]);
}

static int cmpd(const void*a,const void*b){double x=*(const double*)a,y=*(const double*)b;return x<y?-1:x>y;}

int main(int argc, char **argv){
    u64 NB = argc > 1 ? strtoull(argv[1], 0, 0) : 100000000ULL;   /* blocks */
    int reps = argc > 2 ? atoi(argv[2]) : 5;
    aes_core_init(); ttab_init();
    u8 key[16]; for(int i = 0; i < 16; i++) key[i] = (u8)(37*i + 11);
    aesctx ctx; aes_set_key(&ctx, key, 10);
    ttab_keysched(key); aesni_keysched(&ctx);
    /* --- verification vs the byte-wise reference (10^4 random blocks) --- */
    {
        u64 s = 0x1234abcdULL; int bad_tt = 0, bad_ni = 0;
        for(int t = 0; t < 10000; t++){
            u8 P[16], Cref[16], Ctt[16], Cni[16];
            for(int i = 0; i < 16; i++){ s = s*6364136223846793005ULL + 1442695040888963407ULL; P[i] = (u8)(s >> 56); }
            aes_enc(&ctx, P, Cref);
            ttab_enc(P, Ctt);
            _mm_storeu_si128((__m128i*)Cni, aesni_enc1(_mm_loadu_si128((const __m128i*)P)));
            if(memcmp(Cref, Ctt, 16)) bad_tt++;
            if(memcmp(Cref, Cni, 16)) bad_ni++;
        }
        printf("[verify] 10000 random blocks: T-table == reference %s (mismatch %d); AES-NI == reference %s (mismatch %d)\n",
               bad_tt?"FAIL":"PASS", bad_tt, bad_ni?"FAIL":"PASS", bad_ni);
        if(bad_tt || bad_ni) return 1;
    }
    double hz = tsc_hz();
    printf("TSC %.4f GHz; blocks per rep = %llu; reps=%d\n", hz/1e9, (unsigned long long)NB, reps);

    double tt[16], tni8[16], tni1[16], ratios[16];
    u8 sinkb[16]; __m128i sinkv = _mm_setzero_si128();
    for(int rep = 0; rep < reps; rep++){
        double r0 = calib_ratio();
        /* T-table: independent counter blocks (throughput), chained input via counter */
        u8 in[16] = {0}; u8 out[16];
        u64 t0 = rdt();
        for(u64 n = 0; n < NB; n++){
            in[15] = (u8)n; in[14] = (u8)(n >> 8); in[13] = (u8)(n >> 16); in[12] = (u8)(n >> 24);
            ttab_enc(in, out);
            in[0] ^= out[0];   /* light dependency to defeat hoisting; counter dominates */
        }
        u64 t1 = rdt();
        memcpy(sinkb, out, 16);
        /* AES-NI 8-wide pipelined (independent streams) */
        __m128i p[8]; for(int i = 0; i < 8; i++) p[i] = _mm_set1_epi32(i + rep);
        u64 t2 = rdt();
        for(u64 n = 0; n < NB / 8; n++){
            for(int i = 0; i < 8; i++) p[i] = _mm_xor_si128(p[i], NIK[0]);
            for(int r = 1; r < 10; r++) for(int i = 0; i < 8; i++) p[i] = _mm_aesenc_si128(p[i], NIK[r]);
            for(int i = 0; i < 8; i++) p[i] = _mm_aesenclast_si128(p[i], NIK[10]);
        }
        u64 t3 = rdt();
        for(int i = 0; i < 8; i++) sinkv = _mm_xor_si128(sinkv, p[i]);
        /* AES-NI single stream (latency-bound), 10^7 blocks */
        __m128i q = _mm_set1_epi32(rep + 99);
        u64 t4 = rdt();
        for(u64 n = 0; n < NB / 10; n++) q = aesni_enc1(q);
        u64 t5 = rdt();
        sinkv = _mm_xor_si128(sinkv, q);
        double r1 = calib_ratio();
        ratios[rep] = (r0 + r1) * 0.5;
        tt[rep]   = (double)(t1 - t0) / NB;
        tni8[rep] = (double)(t3 - t2) / NB;
        tni1[rep] = (double)(t5 - t4) / (NB / 10);
        printf("  rep %d: T-table %.2f ticks/block | AES-NI 8-wide %.3f ticks/block | AES-NI 1-stream %.3f ticks/block | core/TSC %.3f->%.3f\n",
               rep, tt[rep], tni8[rep], tni1[rep], r0, r1);
        fflush(stdout);
    }
    double rmean = 0; for(int i = 0; i < reps; i++) rmean += ratios[i]; rmean /= reps;
    qsort(tt, reps, sizeof(double), cmpd); qsort(tni8, reps, sizeof(double), cmpd); qsort(tni1, reps, sizeof(double), cmpd);
    printf("DENOM T-table AES-128:  ticks/block median %.2f [min %.2f max %.2f]  => core-cyc/block %.1f (ratio %.3f)\n",
           tt[reps/2], tt[0], tt[reps-1], tt[reps/2] * rmean, rmean);
    printf("DENOM AES-NI 8-wide:    ticks/block median %.3f [min %.3f max %.3f]  => core-cyc/block %.2f\n",
           tni8[reps/2], tni8[0], tni8[reps-1], tni8[reps/2] * rmean);
    printf("DENOM AES-NI 1-stream:  ticks/block median %.3f [min %.3f max %.3f]  => core-cyc/block %.2f\n",
           tni1[reps/2], tni1[0], tni1[reps-1], tni1[reps/2] * rmean);
    printf("sink %02x%02x %llx\n", sinkb[0], sinkb[1], (unsigned long long)_mm_extract_epi64(sinkv, 0));
    return 0;
}
