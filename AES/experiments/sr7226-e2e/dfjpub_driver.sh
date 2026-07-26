#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
# DFJ-as-published end-to-end run.  Workspace: the current directory.
# Phases: compile -> slice verify (graycheck) -> cold DFJ build (timed, 80 thr)
#  -> gates (spot x3, key recovery x3) -> hoisted DFJ rebuild (timed, 80 thr)
#  -> full-table cmp cold vs hoisted -> published-online 50-key sweep (25-way)
#  -> same-machine optimized-DFJ sweep -> ours table build (timed) + ours sweep.
set -u
cd "${WORKDIR:-$(pwd)}"
D="${WORKDIR:-$(pwd)}"
COLD=$D/bloom_dfj_cold_k7.bin
HOIST=$D/bloom_dfj_hoist_k7.bin
OURS=$D/bloom_ours_k7.bin
THR=80
PAR=25
export OMP_NUM_THREADS=$THR
log(){ echo "[$(date '+%F %T')] $*" | tee -a $D/driver.log; }

log "phase 0: compile"
gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_dfjpub sr7226_dfjpub.c -lm || exit 1
gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_gray   sr7226_gray.c   -lm || exit 1
nproc | tee -a $D/driver.log
lscpu | grep -E "Model name|^CPU\(s\)|Thread" | tee -a $D/driver.log

log "phase 1: slice verification (graycheck, 32 full work units x 6 builders, zero tolerance)"
./sr7226_gray graycheck 32 1 32 > $D/graycheck32.log 2> $D/graycheck32.err
grep -q "ALL MULTISETS IDENTICAL" $D/graycheck32.log || { log "SLICE CHECK FAILED - abort"; exit 2; }
log "slice check: $(grep 'GRAYCHECK RESULT' $D/graycheck32.log)"

log "phase 2: COLD (published, naive) DFJ multiset table build, $THR threads"
./sr7226_dfjpub buildcold $COLD 41 0 4032 $THR > $D/build_cold.log 2> $D/build_cold.err
log "cold build done: $(grep '^build:' $D/build_cold.log)"

log "phase 3a: gate (b) spot checks (3 seeds, stride 8192)"
mkdir -p $D/gates
for s in 301 302 303; do
  ./sr7226_dfjpub attack_dfj_pub dummy 41 seed=$s spot=8192 > $D/gates/spot_$s.log 2>&1
done
log "phase 3b: gate (a) key recovery, attack_dfj_pub on 3 seeds (601 602 603), cold table"
cat $COLD > /dev/null
for s in 601 602 603; do
  ( ./sr7226_dfjpub attack_dfj_pub $COLD 41 seed=$s > $D/gates/pub_$s.log 2>&1; echo "$s rc=$?" >> $D/gates/rc.txt ) &
done
wait
GA=$(grep -l "key recovered   : YES" $D/gates/pub_60*.log | wc -l)
GB=$(grep -h "mismatches" $D/gates/spot_*.log | grep -c "mismatches: 0")
log "gate a recovered=$GA/3   gate b spot-pass=$GB/3"
if [ "$GA" != 3 ] || [ "$GB" != 3 ]; then log "GATES FAILED - abort before sweep"; exit 3; fi

log "phase 4: HOISTED (optimized) DFJ multiset table rebuild for same-machine timing, $THR threads"
./sr7226_dfjpub builddfj $HOIST 41 0 4032 $THR > $D/build_hoist.log 2> $D/build_hoist.err
log "hoist build done: $(grep '^build:' $D/build_hoist.log)"
log "phase 4b: full-table bit compare cold vs hoisted (2 x 256 GB)"
if cmp $COLD $HOIST > $D/cmp_cold_hoist.log 2>&1; then
  log "FULL TABLE CMP: IDENTICAL (cold == hoisted, 2^41 bits)"; echo IDENTICAL > $D/cmp_result.txt
else
  log "FULL TABLE CMP: DIFFER !!"; echo DIFFER > $D/cmp_result.txt
fi
rm -f $HOIST

log "phase 5: published-online 50-key sweep, seeds 301..350, $PAR concurrent (attack_dfj_pub)"
cat $COLD > /dev/null
mkdir -p $D/sweep_pub
printf '%s\n' $(seq 301 350) | xargs -P $PAR -I{} bash -c './sr7226_dfjpub attack_dfj_pub '"$COLD"' 41 seed={} > sweep_pub/key_{}.log 2>&1; echo "{} rc=$?" >> sweep_pub/rc.txt'
echo PUB_DONE >> $D/sweep_pub/rc.txt
log "pub sweep done"

log "phase 6: same-machine optimized-DFJ online sweep (attack_dfj), seeds 301..350, $PAR concurrent"
mkdir -p $D/sweep_dfjopt
printf '%s\n' $(seq 301 350) | xargs -P $PAR -I{} bash -c './sr7226_dfjpub attack_dfj '"$COLD"' 41 seed={} > sweep_dfjopt/key_{}.log 2>&1; echo "{} rc=$?" >> sweep_dfjopt/rc.txt'
echo DFJOPT_DONE >> $D/sweep_dfjopt/rc.txt
log "dfj-opt sweep done"

log "phase 7: OURS chi table build re-time on this machine, $THR threads"
./sr7226_dfjpub build $OURS 41 0 4032 $THR > $D/build_ours.log 2> $D/build_ours.err
log "ours build done: $(grep '^build:' $D/build_ours.log)"
log "phase 7b: same-machine ours online sweep (attack_full), seeds 301..350, $PAR concurrent"
cat $OURS > /dev/null
mkdir -p $D/sweep_ours
printf '%s\n' $(seq 301 350) | xargs -P $PAR -I{} bash -c './sr7226_dfjpub attack_full '"$OURS"' 41 seed={} > sweep_ours/key_{}.log 2>&1; echo "{} rc=$?" >> sweep_ours/rc.txt'
echo OURS_DONE >> $D/sweep_ours/rc.txt
log "ours sweep done"
rm -f $OURS
log "ALL DONE"
echo ALL_DONE > $D/ALL_DONE
