# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
CC=gcc
CFLAGS=-O3 -march=native -I../common
online_bench: online_bench.c $(wildcard ../common/*.h)
	$(CC) $(CFLAGS) -o $@ online_bench.c -lm
clean:
	rm -f online_bench
