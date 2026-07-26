#!/bin/bash
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
# OUR 50-key sweep (seeds 301..350, <=25 parallel) with the corrected canonicalizer
# against the new table.  Same protocol as the earlier Bloom sweeps.
set -x
cd "${WORKDIR:-$(pwd)}"
mkdir -p sweep_k7_canon2
cat bloom_ours_canon2_k7.bin > /dev/null          # warm page cache
printf '%s\n' $(seq 301 350) | xargs -P 25 -I{} bash -c './sr7226_canon2 attack_full bloom_ours_canon2_k7.bin 41 seed={} > sweep_k7_canon2/key_{}.log 2>&1; echo "{} rc=$?" >> sweep_k7_canon2/rc.txt'
echo OURS_CANON2_DONE >> sweep_k7_canon2/rc.txt
