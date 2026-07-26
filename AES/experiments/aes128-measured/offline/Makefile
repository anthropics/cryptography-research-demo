# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
CC=gcc
CFLAGS=-O3 -march=native -I../common -pthread -Wall -Wno-unused-function
HDRS=$(wildcard ../common/*.h)
all: table_build analyze classify
table_build: table_build.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ table_build.c -lm
analyze: analyze.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ analyze.c -lm
classify: classify.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ classify.c -lm
clean:
	rm -f table_build analyze classify
