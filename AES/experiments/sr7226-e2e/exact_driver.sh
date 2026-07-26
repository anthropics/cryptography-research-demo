#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
set -x
cd "${WORKDIR:-$(pwd)}"
export OMP_NUM_THREADS=190
gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_exact2 sr7226_exact2.c -lm -lpthread
gcc -O3 -fopenmp -o sort_buckets sort_buckets.c
# ---- OURS: exact chi table ----
mkdir -p tbl_chi
taskset -c 0-207 ./sr7226_exact2 buildexact tbl_chi 0 4032 190 > build_chi.log 2>&1
taskset -c 0-207 ./sort_buckets tbl_chi 64 > sort_chi.log 2>&1
du -sh tbl_chi
# warm page cache
cat tbl_chi/*.sorted > /dev/null
mkdir -p sweep_exact_ours
printf '%s\n' $(seq 301 350) | xargs -P 50 -I{} bash -c 'taskset -c 0-207 ./sr7226_exact2 attack_exact tbl_chi seed={} > sweep_exact_ours/key_{}.log 2>&1; echo "{} rc=$?" >> sweep_exact_ours/rc.txt'
echo OURS_DONE >> sweep_exact_ours/rc.txt
# ---- DFJ: exact multiset-hash table ----
mkdir -p tbl_dfj
taskset -c 0-207 ./sr7226_exact2 buildexact tbl_dfj 0 4032 190 dfj > build_dfj.log 2>&1
taskset -c 0-207 ./sort_buckets tbl_dfj 64 > sort_dfj.log 2>&1
du -sh tbl_dfj
cat tbl_dfj/*.sorted > /dev/null
mkdir -p sweep_exact_dfj
printf '%s\n' $(seq 301 350) | xargs -P 50 -I{} bash -c 'taskset -c 0-207 ./sr7226_exact2 attack_exact_dfj tbl_dfj seed={} > sweep_exact_dfj/key_{}.log 2>&1; echo "{} rc=$?" >> sweep_exact_dfj/rc.txt'
echo DFJ_DONE >> sweep_exact_dfj/rc.txt
echo ALL_DONE > ALL_DONE
