#!/bin/sh
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0

# Run the 13-round LEA-128 key-recovery attack with a bounded escalation ladder.
#
# usage:
#   ./run.sh                  attack the HTTP oracle at ORACLE_HOST:ORACLE_PORT
#                             (default 127.0.0.1:8666; see oracle_server.py)
#   ./run.sh <32 hex chars>   self-test against a built-in local oracle
#                             simulating that key (no server needed)
#   ./run.sh random           self-test with a fresh random key
#
# Settings can be overridden via environment variables:
#   N1 (default 2^20), N2 (2^27), BEAM (64), THREADS (nproc, capped at 1024),
#   MAXBIT (24), SEED (1).
# Exit codes: 0 key recovered; 1 not found after the ladder; 2 bad parameters;
# 3 resource allocation failure (memory or threads); 4 oracle communication failure;
# 128+N killed by signal N (e.g. 137 = OOM killer).
#
# If a run reports not-found, up to two escalated attempts are made, each with
# FRESH data (new seed) and a much wider beam. Beam width costs time, not
# memory (the engine stores 32 bytes per ciphertext pair regardless of width),
# so the rungs differ mainly in runtime; only the final rung, which also
# doubles the phase-2 data, needs more RAM (about 8.0 GiB vs 4.0 GiB).
# About 44% of random-key draws need at least one escalation; see README.
#   attempt 2: 2x N1, 8x BEAM
#   attempt 3: 2x N1, 2x N2, 16x BEAM

set -u

N1=${N1:-1048576}
N2=${N2:-134217728}
BEAM=${BEAM:-64}
THREADS=${THREADS:-$(nproc 2>/dev/null || echo 8)}
[ "$THREADS" -gt 1024 ] && THREADS=1024
MAXBIT=${MAXBIT:-24}
SEED=${SEED:-1}

KEY=${1:-}
if [ "$KEY" = "random" ]; then
    KEY=$(od -An -N16 -tx1 /dev/urandom | tr -d ' \n')
    echo "generated random key: $KEY"
fi

# returns 0 = key found, 1 = attack reported failure (escalate),
# aborts the ladder if the attack was killed by a signal (e.g. the kernel's
# OOM killer at exit code 137) - escalating would only use more memory.
run_once() {
    echo "=== attack attempt: N1=$1 N2=$2 BEAM=$3 THREADS=$THREADS MAXBIT=$MAXBIT SEED=$4 ==="
    if [ -n "$KEY" ]; then
        ./attack "$1" "$2" "$3" "$THREADS" "$MAXBIT" "$4" "$KEY"
    else
        ./attack "$1" "$2" "$3" "$THREADS" "$MAXBIT" "$4"
    fi
    rc=$?
    if [ "$rc" -eq 2 ]; then
        echo ""
        echo "ATTACK ABORTED: invalid parameters (see the usage message above)."
        exit 2
    fi
    if [ "$rc" -eq 4 ]; then
        echo ""
        echo "ATTACK ABORTED: oracle communication failed (is the oracle server"
        echo "running and reachable at ORACLE_HOST:ORACLE_PORT?)."
        exit 4
    fi
    if [ "$rc" -eq 3 ]; then
        echo ""
        echo "ATTACK ABORTED: resource allocation failed (memory or threads)."
        echo "If memory: this machine cannot fit the attempted settings (estimate"
        echo "printed at phase 2 start); escalation is not attempted because it"
        echo "would need even more memory."
        exit 3
    fi
    if [ "$rc" -ge 128 ]; then
        echo ""
        echo "ATTACK ABORTED: the attack process was killed by signal $((rc - 128))."
        echo "If this was the out-of-memory killer, this machine cannot fit the"
        echo "attempted settings (the memory estimate is printed at phase 2 start);"
        echo "escalation is not attempted because it would need even more memory."
        exit "$rc"
    fi
    return "$rc"
}

run_once "$N1" "$N2" "$BEAM" "$SEED" && exit 0

echo ""
echo "=== attempt 1 failed; escalating: fresh data, 2x phase-1 data, 8x beam width ==="
run_once $((N1 * 2)) "$N2" $((BEAM * 8)) $((SEED + 1)) && exit 0

echo ""
echo "=== attempt 2 failed; escalating: fresh data, 2x phase-1 and phase-2 data, 16x beam width ==="
run_once $((N1 * 2)) $((N2 * 2)) $((BEAM * 16)) $((SEED + 2)) && exit 0

echo ""
echo "ATTACK FAILED after escalation. Attempts:"
echo "  (N1=$N1, N2=$N2, BEAM=$BEAM), (N1=$((N1*2)), N2=$N2, BEAM=$((BEAM*8))), (N1=$((N1*2)), N2=$((N2*2)), BEAM=$((BEAM*16)))"
echo "Keys that survive this ladder have long runs of zero-signal key-bit"
echo "positions (see README: dead zones); wider beams or more data may help"
echo "(BEAM/N2 env vars), but some such keys need a qualitatively different"
echo "statistic and are not recoverable by this package."
exit 1
