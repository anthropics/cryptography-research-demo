# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

CC      ?= cc
# -O3: the scoring core is a vectorizable word-recompute loop; -O2 costs ~4x.
# Optional extra ~2x on modern x86: make CFLAGS="-O3 -march=native"
CFLAGS  ?= -O3
LDLIBS   = -lm

attack: attack.c lea.h oracle.h
	$(CC) $(CFLAGS) -pthread -o attack attack.c $(LDLIBS)

check: kat.c lea.h
	$(CC) $(CFLAGS) -o kat kat.c $(LDLIBS)
	./kat


clean:
	rm -f attack kat

.PHONY: check clean
