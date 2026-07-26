// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
/* classify.c -- resolve what a duplicate / word-collision IS, from the
 * regenerated entries emitted by `table_build resolve` (TAB/KA lines with
 * full sequences).  For every collision group (records sharing fp0, or fp1,
 * or both), compares the underlying d-sequences:
 *   SEQ   : identical ordered sequences (distinct provenance => distinct
 *           parameter tuples giving the identical difference sequence)
 *   MSET  : same multiset, different order (same fingerprint by design)
 *   CHI   : same parity set chi(D) but different multisets (elements of
 *           even multiplicity cancel; same fingerprint by design)
 *   AGL   : different parity sets that are AGL(1,256)-equivalent
 *           (canonical sets coincide => same canonical hash by design)
 *   HASH  : none of the above: distinct canonical sets whose 64-bit FNV
 *           hashes collide (a genuine 64-bit hash collision)
 *
 * Usage: classify < resolve_output.txt   (lines starting TAB/KA)
 * Build: gcc -O3 -march=native -pthread -o classify classify.c -lm
 */
#include "common.h"

typedef struct {
    u64 fp0, fp1;
    u8 E[256];
    char tag[8];
    unsigned base, gidx;
    chi_t chiD;        /* parity set of d = 1/e */
    chi_t c0star;      /* canonical odd set (fp0) */
    chi_t c1star;      /* canonical even set (fp1) */
    u64 f0chk, f1chk;
} ent_t;

/* recompute fingerprints AND the canonical sets (mirrors kernel_full_v2) */
static void kernel_sets(const u8 *U0, u64 *fp0, u64 *fp1, chi_t *c0out, chi_t *c1out){
    u8 d[256]; u32 Sk = 0;
    for(int om = 1; om < 256; om++){ u8 dv = GF_inv[U0[om]]; d[om] = dv; Sk ^= POW1357[dv]; }
    u8 S1 = Sk & 0xff, S3 = (Sk >> 8) & 0xff, S5 = (Sk >> 16) & 0xff, S7 = (Sk >> 24) & 0xff;
    u8 a0 = alpha_from_S(S1, S3, S5, S7, 1, d + 1, 255);
    u8 b0 = gmul(a0, S1);
    chi_t c0; chi_clr(&c0);
    const u8 *am0 = MUL[a0];
    for(int om = 1; om < 256; om++) chi_flip(&c0, (u8)(am0[d[om]] ^ b0));
    *fp0 = chi_hash(&c0); *c0out = c0;
    u8 S2 = gmul(S1, S1), S4 = gmul(S2, S2), S6 = gmul(S3, S3);
    if(S1 == 0){
        chi_t ch; chi_clr(&ch);
        for(int om = 1; om < 256; om++) chi_flip(&ch, d[om]);
        chi_flip(&ch, 0);
        *fp1 = canon_even(&ch, c1out);
        return;
    }
    u8 P7_1 = (u8)(gmul(S1, S6) ^ gmul(S2, S5) ^ gmul(S3, S4));
    u8 a1;
    if(P7_1){ int l = GF_log[P7_1]; a1 = GF_exp[(255 - (l * 73) % 255) % 255]; }
    else {
        u32 Tk = 0; for(int om = 1; om < 256; om++) Tk ^= POW9_11_13[d[om]];
        u8 S9 = (u8)Tk, S8 = gmul(S4, S4), S10 = gmul(S5, S5);
        u8 P11 = (u8)(gmul(S1, S10) ^ gmul(S2, S9) ^ gmul(S3, S8));
        if(P11){ int l = GF_log[P11]; a1 = GF_exp[(255 - (l * 116) % 255) % 255]; }
        else {
            chi_t ch; chi_clr(&ch);
            for(int om = 1; om < 256; om++) chi_flip(&ch, d[om]);
            chi_flip(&ch, 0);
            a1 = kv2_alpha_deep(&ch);
            if(!a1){ *fp1 = canon_bruteAGL(&ch, c1out); return; }
        }
    }
    u8 A = gmul(S3, GF_inv[gmul(S1, S2)]);
    u8 u = CE_UROOT[A];
    u8 aS1 = gmul(a1, S1), ba = gmul(u, aS1);
    const u8 *am1 = MUL[a1];
    chi_t ca; chi_clr(&ca);
    for(int om = 1; om < 256; om++) chi_flip(&ca, (u8)(am1[d[om]] ^ ba));
    chi_flip(&ca, ba);
    chi_t cb = ca;
    for(int k = 0; k < 8; k++) if(aS1 & (1 << k)) chi_xorshift_bit(&cb, k);
    u64 ha = chi_hash(&ca), hb = chi_hash(&cb);
    if(ha < hb){ *fp1 = ha; *c1out = ca; } else { *fp1 = hb; *c1out = cb; }
}
static int chi_eq(const chi_t *a, const chi_t *b){
    return a->w[0]==b->w[0] && a->w[1]==b->w[1] && a->w[2]==b->w[2] && a->w[3]==b->w[3];
}
/* brute-force AGL(1,256) equivalence of two sets */
static int agl_equiv(const chi_t *A, const chi_t *B){
    if(ce_popcount(A) != ce_popcount(B)) return 0;
    for(int al = 1; al < 256; al++){
        chi_t S; ce_scale(A, (u8)al, &S);
        /* find beta mapping min element: try all beta via bit shifts */
        for(int bi = 0; bi < 256; bi++){
            if(bi){ int k = __builtin_ctz(bi); chi_xorshift_bit(&S, k); }
            if(chi_eq(&S, B)) return 1;
        }
    }
    return 0;
}

static ent_t *E; static int NE = 0, CE = 0;
static int parse_hex(const char *s, u8 *out, int n){
    for(int i = 0; i < n; i++){ unsigned v; if(sscanf(s + 2*i, "%2x", &v) != 1) return -1; out[i] = (u8)v; }
    return 0;
}
static int cmp0(const void *a, const void *b){ const ent_t *x = a, *y = b; return x->fp0 < y->fp0 ? -1 : x->fp0 > y->fp0; }
static int cmp1(const void *a, const void *b){ const ent_t *x = a, *y = b; return x->fp1 < y->fp1 ? -1 : x->fp1 > y->fp1; }

static const char *pair_class(const ent_t *a, const ent_t *b, int word){
    if(!memcmp(a->E + 1, b->E + 1, 255)) return "SEQ(identical sequence, distinct provenance)";
    /* multiset */
    u8 c1[256], c2[256]; memset(c1, 0, 256); memset(c2, 0, 256);
    for(int i = 1; i < 256; i++){ c1[GF_inv[a->E[i]]]++; c2[GF_inv[b->E[i]]]++; }
    if(!memcmp(c1, c2, 256)) return "MSET(same multiset, permuted sequence)";
    if(chi_eq(&a->chiD, &b->chiD)) return "CHI(same parity set, different multisets)";
    /* canonical set comparison for the relevant word */
    const chi_t *x = word ? &a->c1star : &a->c0star, *y = word ? &b->c1star : &b->c0star;
    if(chi_eq(x, y)){
        /* same canonical set => same AGL class of the underlying (odd/even) set */
        return word ? "AGL(even sets chi(D)u{0} affinely equivalent => fp1 equal by design)"
                    : "AGL(odd sets chi(D) affinely equivalent => fp0 equal by design)";
    }
    /* confirm with brute force (should agree) */
    chi_t W1 = a->chiD, W2 = b->chiD;
    if(word){ chi_flip(&W1, 0); chi_flip(&W2, 0); }
    if(agl_equiv(&W1, &W2)) return "AGL*(brute-force equivalent; canonical sets differ -- should not happen)";
    return "HASH(distinct canonical sets, genuine 64-bit FNV collision)";
}

int main(int argc, char **argv){
    (void)argc; (void)argv;
    aes_core_init(); ddt_init(); chifast_init(); kernel_v2_init();
    CE = 1 << 14; E = (ent_t*)malloc(sizeof(ent_t) * CE);
    char line[8192];
    while(fgets(line, sizeof line, stdin)){
        if(strncmp(line, "TAB", 3) && strncmp(line, "KA", 2)) continue;
        if(NE == CE){ CE *= 2; E = (ent_t*)realloc(E, sizeof(ent_t) * CE); }
        ent_t *e = &E[NE];
        memset(e, 0, sizeof *e);
        strncpy(e->tag, line, 3); e->tag[3] = 0;
        char *q;
        if((q = strstr(line, "base=")) != NULL) e->base = (unsigned)strtoul(q + 5, 0, 10);
        if((q = strstr(line, "gidx=")) != NULL) e->gidx = (unsigned)strtoul(q + 5, 0, 10);
        if((q = strstr(line, "fp0=")) == NULL) continue;  e->fp0 = strtoull(q + 4, 0, 16);
        if((q = strstr(line, "fp1=")) == NULL) continue;  e->fp1 = strtoull(q + 4, 0, 16);
        if((q = strstr(line, "seq=")) == NULL) continue;
        q += 4; e->E[0] = 0;
        if(parse_hex(q, e->E + 1, 255)) continue;
        /* recompute canonical data */
        kernel_sets(e->E, &e->f0chk, &e->f1chk, &e->c0star, &e->c1star);
        chi_of_seq(e->E, &e->chiD);
        NE++;
    }
    printf("[classify] %d resolved entries read\n", NE);
    int bad = 0;
    for(int i = 0; i < NE; i++) if(E[i].f0chk != E[i].fp0 || E[i].f1chk != E[i].fp1) bad++;
    printf("[classify] recomputed (fp0,fp1) from emitted sequences agree with the stored records for %d/%d\n", NE - bad, NE);

    /* group by fp0 */
    qsort(E, NE, sizeof(ent_t), cmp0);
    int g0 = 0, g0pairs = 0;
    for(int i = 0; i < NE; ){
        int j = i + 1; while(j < NE && E[j].fp0 == E[i].fp0) j++;
        if(j - i >= 2){
            g0++;
            for(int a = i; a < j; a++) for(int b = a + 1; b < j; b++){
                g0pairs++;
                int dup = (E[a].fp1 == E[b].fp1);
                printf("FP0-GROUP fp0=%016llx %s {%s b=%u i=%u} {%s b=%u i=%u}%s : %s\n",
                       (unsigned long long)E[i].fp0, dup ? "FULL-PAIR-DUP" : "fp0-word-collision",
                       E[a].tag, E[a].base, E[a].gidx, E[b].tag, E[b].base, E[b].gidx,
                       (E[a].base == E[b].base && E[a].gidx == E[b].gidx) ? " [SAME PROVENANCE]" : "",
                       pair_class(&E[a], &E[b], 0));
            }
        }
        i = j;
    }
    printf("[classify] fp0 groups: %d (pairs %d)\n", g0, g0pairs);
    /* group by fp1 */
    qsort(E, NE, sizeof(ent_t), cmp1);
    int g1 = 0, g1pairs = 0;
    for(int i = 0; i < NE; ){
        int j = i + 1; while(j < NE && E[j].fp1 == E[i].fp1) j++;
        if(j - i >= 2){
            g1++;
            for(int a = i; a < j; a++) for(int b = a + 1; b < j; b++){
                g1pairs++;
                int dup = (E[a].fp0 == E[b].fp0);
                printf("FP1-GROUP fp1=%016llx %s {%s b=%u i=%u} {%s b=%u i=%u}%s : %s\n",
                       (unsigned long long)E[i].fp1, dup ? "FULL-PAIR-DUP" : "fp1-word-collision",
                       E[a].tag, E[a].base, E[a].gidx, E[b].tag, E[b].base, E[b].gidx,
                       (E[a].base == E[b].base && E[a].gidx == E[b].gidx) ? " [SAME PROVENANCE]" : "",
                       pair_class(&E[a], &E[b], 1));
            }
        }
        i = j;
    }
    printf("[classify] fp1 groups: %d (pairs %d)\n", g1, g1pairs);
    return 0;
}
