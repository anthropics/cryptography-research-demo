#!/bin/sh
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# Validation battery: run the full attack (run.sh, including its escalation
# rule) on fresh random keys and log one line per key.
#
# usage: ./validate.sh [num_keys]          (default 10)
# output: validate.log in the current directory (fresh per invocation)

set -u
COUNT=${1:-10}
: > validate.log
# Quick build/health check on a fixed known-good key at default settings
# BEFORE burning hours on random keys: this key recovers at the first rung.
# (A fresh RANDOM key legitimately fails the first rung on about 44% of
# draws - the ladder handles those; see README. A failure on THIS key at
# the first rung indicates a broken build or environment.)
FIXEDKEY=149fe614a07b0049e8187f52268908bb
VTHREADS=$(nproc 2>/dev/null || echo 8)
[ "$VTHREADS" -gt 1024 ] && VTHREADS=1024
START=$(date +%s)
if ./attack 1048576 134217728 64 "$VTHREADS" 24 1 "$FIXEDKEY" > validate.key0.log 2>&1; then
    echo "key=$FIXEDKEY result=PASS seconds=$(($(date +%s) - START)) (fixed-key health check)" | tee -a validate.log
else
    echo "key=$FIXEDKEY result=FAIL (fixed-key health check - broken build?)" | tee -a validate.log
    exit 1
fi
i=0
while [ "$i" -lt "$COUNT" ]; do
    i=$((i + 1))
    KEY=$(od -An -N16 -tx1 /dev/urandom | tr -d ' \n')
    START=$(date +%s)
    if ./run.sh "$KEY" > "validate.key$i.log" 2>&1; then RES=PASS; else RES=FAIL; fi
    END=$(date +%s)
    ESC=$(grep -c "escalating:" "validate.key$i.log")
    echo "key=$KEY result=$RES seconds=$((END - START)) escalated=$ESC" | tee -a validate.log
done
