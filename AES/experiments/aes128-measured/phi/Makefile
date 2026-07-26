# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
CC=gcc
CFLAGS=-O3 -march=native -I../common -pthread -Wall -Wno-unused-function -Wno-unused-variable
HDRS=$(wildcard ../common/*.h)
phi: phi.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ phi.c -lm
clean:
	rm -f phi
