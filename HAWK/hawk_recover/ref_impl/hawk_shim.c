// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0

// Minimal shim: sign with decoded sk, verify with pk.
// gcc -o hawk_shim hawk_shim.c hawk_sign.o hawk_vrfy.o hawk_kgen.o ng_*.o sha3.o -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hawk.h"

// Read exactly a whole file into buf. Returns the byte count, -1 if the file
// cannot be opened, -2 if it is longer than len (a wrong-size input must not
// be silently truncated to a plausible length).
static int read_file(const char *fn, void *buf, size_t len) {
  FILE *f = fopen(fn, "rb"); if (!f) return -1;
  size_t r = fread(buf, 1, len, f);
  int extra = fgetc(f); fclose(f);
  if (extra != EOF) return -2;
  return (int)r;
}

// RNG: use /dev/urandom
static void rng_urandom(void *ctx, void *dst, size_t len) {
  (void)ctx;
  FILE *f = fopen("/dev/urandom", "rb");
  if (!f || fread(dst, 1, len, f) != len) { perror("/dev/urandom"); exit(1); }
  fclose(f);
}

int main(int argc, char **argv) {
  unsigned logn = 8;
  unsigned n = 1u << logn;
  if (argc < 4) { fprintf(stderr, "usage: %s <cmd> <pk.bin> <skd592.bin> [msg]\n", argv[0]); return 1; }
  const char *cmd = argv[1], *pkf = argv[2], *skf = argv[3];
  const char *msg = argc > 4 ? argv[4] : "HAWK-256 key recovery test message";

  if (!strcmp(cmd, "hpub")) {
    // compute hpub = SHAKE256(pk)[:16]
    unsigned char pk[HAWK_PUBKEY_SIZE(8)];
    int pklen = read_file(pkf, pk, sizeof pk);
    if (pklen != (int)sizeof pk) { fprintf(stderr, "pk read %d != %zu\n", pklen, sizeof pk); return 1; }
    shake_context sc; shake_init(&sc, 256);
    shake_inject(&sc, pk, pklen); shake_flip(&sc);
    unsigned char hpub[16]; shake_extract(&sc, hpub, 16);
    fwrite(hpub, 1, 16, stdout);
    fprintf(stderr, "hpub written (16 bytes)\n");
    return 0;
  }

  if (!strcmp(cmd, "signverify")) {
    unsigned char pk[HAWK_PUBKEY_SIZE(8)];
    unsigned char skd[HAWK_PRIVKEY_DECODED_SIZE(8)];
    unsigned char sig[HAWK_SIG_SIZE(8)];
    unsigned char tmp[HAWK_TMPSIZE_VERIFY(8) + HAWK_TMPSIZE_SIGN(8)];
    int pklen = read_file(pkf, pk, sizeof pk);
    int sklen = read_file(skf, skd, sizeof skd);
    fprintf(stderr, "pk=%d (expect %zu), skd=%d (expect %zu)\n", pklen, sizeof pk, sklen, sizeof skd);
    if (pklen != (int)sizeof pk || sklen != (int)sizeof skd) return 1;

    // Sign
    shake_context sc; hawk_sign_start(&sc);
    shake_inject(&sc, msg, strlen(msg));
    int s = hawk_sign_finish_alt(logn, rng_urandom, NULL, sig, &sc, skd, sizeof skd, tmp, sizeof tmp);
    fprintf(stderr, "sign_finish_alt = %d\n", s);
    if (!s) { fprintf(stderr, "SIGN FAILED\n"); return 2; }

    // Verify
    shake_context scv; hawk_verify_start(&scv);
    shake_inject(&scv, msg, strlen(msg));
    int v = hawk_verify_finish(logn, sig, sizeof sig, &scv, pk, sizeof pk, tmp, sizeof tmp);
    fprintf(stderr, "verify_finish = %d\n", v);
    if (!v) { fprintf(stderr, "VERIFY FAILED\n"); return 3; }

    fprintf(stderr, "*** SIGN(recovered_sk)/VERIFY(original_pk) = PASS ***\n");
    // dump sig hex
    fprintf(stderr, "sig[:32]: ");
    for (int i=0;i<32;i++) fprintf(stderr, "%02x", sig[i]);
    fprintf(stderr, "\n");
    return 0;
  }
  return 1;
}
