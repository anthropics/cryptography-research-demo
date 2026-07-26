#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
# Full offline-table pipeline: build 2^36 -> sort/dedup/stats -> fp1 stats -> probe 2^34 -> KA -> dup resolution+classify
# Output goes to ROOT (default: this dir; ROOT=/big/disk for large runs).
set -x
SRC=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-$SRC}
B=$ROOT/buckets
LOG=$ROOT/logs
mkdir -p $B $LOG
cd $ROOT
echo "=== START $(date -u) ===" >> $LOG/pipeline.log
nproc --all; free -g | head -2; df -h /root | tail -1

# 1. BUILD 2^36 entries = 65536 bases x 2^20 gray-walked entries + 128 planted known answers
$SRC/table_build build nb=65536 threads=190 nka=128 dir=$B > $LOG/build.out 2> $LOG/build.log
echo "=== BUILD DONE $(date -u) rc=$? ===" >> $LOG/pipeline.log
du -sh $B >> $LOG/pipeline.log
df -h /root | tail -1 >> $LOG/pipeline.log

# 2. ANALYZE: sort+dedup+stats, fp1 stats, probe (2^18 inst x 2^16 wrong k7a = 2^34), known-answer
$SRC/analyze all threads=190 sortpar=48 pinst=262144 pg=65536 dir=$B > $LOG/analyze.out 2> $LOG/analyze.err
echo "=== ANALYZE DONE $(date -u) rc=$? ===" >> $LOG/pipeline.log

# 3. DUPLICATE / COLLISION RESOLUTION: 100-sample of each list -> regenerate provenance+sequences -> classify
cd $ROOT
: > t0.txt; : > t1.txt
[ -s dups.txt ]  && shuf -n 100 --random-source=<(yes) dups.txt  | awk "{print \$1}" >> t0.txt
[ -s coll0.txt ] && shuf -n 100 --random-source=<(yes) coll0.txt | awk "{print \$1}" >> t0.txt
[ -s coll1.txt ] && shuf -n 100 --random-source=<(yes) coll1.txt | awk "{print \$1}" >> t1.txt
wc -l t0.txt t1.txt >> $LOG/pipeline.log
if [ -s t0.txt ] || [ -s t1.txt ]; then
  ARGS=""
  [ -s t0.txt ] && ARGS="$ARGS targets0=$ROOT/t0.txt"
  [ -s t1.txt ] && ARGS="$ARGS targets1=$ROOT/t1.txt"
  $SRC/table_build resolve nb=65536 threads=190 nka=128 $ARGS > $ROOT/resolve.txt 2> $LOG/resolve.log
  echo "=== RESOLVE DONE $(date -u) rc=$? ===" >> $LOG/pipeline.log
  $SRC/classify < $ROOT/resolve.txt > $ROOT/classify.txt 2> $LOG/classify.err
  echo "=== CLASSIFY DONE $(date -u) rc=$? ===" >> $LOG/pipeline.log
else
  echo "=== NO DUPS/COLLISIONS TO RESOLVE ===" >> $LOG/pipeline.log
fi
echo "=== ALL DONE $(date -u) ===" >> $LOG/pipeline.log
