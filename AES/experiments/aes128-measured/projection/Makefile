# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
CC=gcc
CFLAGS=-O3 -march=native -I../common -Wall -Wno-unused-function -Wno-unused-variable
HDRS=$(wildcard ../common/*.h)
all: offline3 online3 denom
%: %.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f offline3 online3 denom
