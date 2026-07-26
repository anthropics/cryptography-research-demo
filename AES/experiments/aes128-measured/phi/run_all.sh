#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
# run_all.sh -- full phi measurement campaign on the benchmark machine.
# Output goes to ROOT (default: this dir; ROOT=/big/disk for the ~950 GB).
# Usage: mkdir -p logs && setsid nohup ./run_all.sh > logs/run_all.log 2>&1 &
set -u
SRC=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-$SRC}
PY=${PY:-python3}
mkdir -p $ROOT/logs $ROOT/wrong $ROOT/ideal $ROOT/table
cd $ROOT
echo "=== start $(date) host $(hostname) nproc $(nproc) ==="
echo "--- gate1 ---"
$SRC/phi gate1 n=10000000 seed=4242
echo "--- gate2 (2e6 instances) ---"
$SRC/phi gate2 n=2000000 threads=96 seed=0xb21d6e
echo "--- table: 2^30 entries = 2^14 bases x 2^16 Gray steps ---"
$SRC/phi table logbases=14 walkbits=16 dir=$ROOT/table threads=96 seed=0x7ab1e5
echo "--- ideal baseline: 2^34 Phi(iid) ---"
$SRC/phi ideal logn=34 dir=$ROOT/ideal threads=96 seed=0x1dea1
echo "--- wrong-candidate side: 2^34 genuine AES-derived ---"
$SRC/phi wrong logn=34 dir=$ROOT/wrong threads=128 seed=0x57a7e1
echo "--- analysis: table ---"
$SRC/phi anat dir=$ROOT/table threads=48
echo "--- analysis: ideal ---"
$SRC/phi anaw dir=$ROOT/ideal threads=32
echo "--- analysis: wrong ---"
$SRC/phi anaw dir=$ROOT/wrong threads=32
echo "--- cross: wrong vs table ---"
$SRC/phi cross wdir=$ROOT/wrong tdir=$ROOT/table threads=32
echo "--- cross: ideal vs table ---"
$SRC/phi cross wdir=$ROOT/ideal tdir=$ROOT/table threads=32
echo "--- cross: wrong vs ideal (extra baseline) ---"
$SRC/phi cross wdir=$ROOT/wrong tdir=$ROOT/ideal threads=32
echo "--- chi-square reports ---"
$PY $SRC/chistats.py $ROOT/wrong/stats.txt $ROOT/ideal/stats.txt $ROOT/table/stats.txt
echo "--- two-sample real-vs-ideal ---"
$PY $SRC/compare2.py $ROOT/wrong/stats.txt $ROOT/ideal/stats.txt
$PY $SRC/compare2.py $ROOT/table/stats.txt $ROOT/ideal/stats.txt
echo "=== done $(date) ==="
du -sh $ROOT/wrong $ROOT/ideal $ROOT/table
