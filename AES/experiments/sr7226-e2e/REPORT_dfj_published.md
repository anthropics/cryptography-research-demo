# "DFJ as published" end-to-end on SR(7,2,2,6) — third column for the paper table

Sources: `sr7226_dfjpub.c`; all runs on **the benchmark node**
(2x Xeon Platinum 8481C, 208 threads, 1.8 TB RAM, idle). Per-run logs, sweep
directories and comparison outputs are produced by the driver and not bundled.

## 1. What the variant is

Source: `sr7226_dfjpub.c` = the DFJ table-attack driver + one new online
mode. `build_table()` is byte-identical to the production file (diff-verified); only the
online branch `attack_dfj_pub` (plus a `spot=` gate and a DFJ self-test entry) was added.

* **Offline (published):** DFJ multiset-hash table, 2^41-bit Bloom, BLOOM_K=7, built with
  `COLD_MODE=1` (`buildcold`): the naive per-entry enumeration of the original builder —
  the a-side sequence y4a[] is rebuilt for every branch tuple (no hoisting over (b1,b2),
  no Gray amortization).
* **Offline (optimized):** same table via `builddfj` (hoisted a-side), re-timed on the same node.
* **Online (published), `attack_dfj_pub`:** identical to `attack_dfj` except inside the
  64-value u_5[0]=s loop the delta-set ciphertext peel `v[w]=T0[ct0[w]^k60]^T1[ct3[w]^k63]`
  is RECOMPUTED for every s (no sharing across the s loop), no early exit on Bloom hits,
  plain scalar. Candidate enumeration, Bloom probes and key-schedule disposal unchanged.

## 2. Verification (zero tolerance) — all PASS

* **Slice multiset identity** (`sr7226_gray graycheck 32 1 32`, node): 32 full work
  units (564,347,136 entries per builder, 0.79% of the table) enumerated by 6 builders
  incl. orig-hoist and orig-cold (= the production `build_table` in both COLD_MODE
  settings): sorted arrays of strong 64-bit hashes of every raw 64-element d-sequence
  identical, entry counts identical, and the commutative accumulators of the emitted
  DFJ multiset-hash (ms_sum/ms_xor) and chi fingerprints identical for every wu and
  builder. `GRAYCHECK RESULT: ALL MULTISETS IDENTICAL`, 0 mismatches (80.6 s).
* **Full-table identity:** the two 256 GB Bloom files (cold build, hoisted rebuild)
  compared with `cmp`: **bit-identical** (cmp_result.txt: IDENTICAL). Both report
  71,282,781,952 inserted d-sequences, load 0.1998.
* **Dev small-slice file identity** (wu 0..3, 2^34-bit filters, 80 thr): cold and hoisted
  files bit-identical, 69,553,152 entries each.
* **Gate (a)** — `attack_dfj_pub` recovers the key on 3 seeds (601, 602, 603) against the
  cold table: 3/3 YES (820.96 s, 539.68 s, 1575.32 s; 3 concurrent). PASS.
* **Gate (b)** — per-candidate multiset-hash stream of `attack_dfj_pub` (recomputed peel)
  vs `attack_dfj` (shared peel): seeds 301/302/303, random-phase stride-8192 sampling of
  the candidate stream, all 64 s per sampled candidate: 1017 / 1001 / 998 candidates =
  65088 / 64064 / 63872 (candidate,s) hash comparisons, **0 mismatches** each. PASS.

## 3. Table-build wall clocks (2^41-bit Bloom, BLOOM_K=7, 4032 work units, 80 threads)

| table | machine | build time (build_table) | entries/s | entries |
|---|---|---:|---:|---:|
| **DFJ, COLD / published** (buildcold) | the benchmark node | **934.3 s** (15.6 min) | 76.3 M/s | 71,282,781,952 |
| **DFJ, HOISTED / optimized** (builddfj) | the benchmark node | **823.7 s** (13.7 min) | 86.5 M/s | 71,282,781,952 |
| DFJ, HOISTED (earlier run, bloom_dfj_k7.bin) | local box (loaded) | 1397.9 s | 51.0 M/s | same |
| ours (chi, 2 parities), `build` re-timed | the benchmark node | 1618.6 s (27.0 min) | 44.0 M/s | same |
| ours (earlier run, bloom_7226_k7.bin) | local box (loaded) | 2745.6 s | 26.0 M/s | same |

Build times are the program's own `build_table` wall clock (excludes the ~3-4 min
256 GB file save; process spans from driver.log: cold 11:38:23->11:58:22, hoisted
12:25:39->12:43:46). Cold/hoisted = 1.13x on the node at 80 threads: the random
2^41-bit Bloom inserts (7 per entry) dominate the multithreaded build and mask most of
the sequence-generation saving (single-thread seq-gen gap measured earlier was 1.43 M
vs 1.89 M entries/s).

## 4. Online 50-key sweeps, seeds 301..350, concurrency = 25 processes (xargs -P 25)

All three sweeps use the same 25-way parallelism. The published and optimized-DFJ
sweeps both ran against the cold table (bit-identical to the hoisted one). Per-key logs
and aggregate tables are written by the driver and `aggregate.sh`; not bundled.

| online variant | machine | recovered | lookups mean / median | wall s/key mean / median / min / max |
|---|---|---:|---:|---:|
| **DFJ as published** (attack_dfj_pub) | the benchmark node | **50/50** | 1.373e9 / 1.045e9 | **757.3 / 571.0 / 274.8 / 1927.6** |
| DFJ + our online optimizations (attack_dfj) | the benchmark node | 50/50 | 1.373e9 / 1.045e9 | 680.0 / 517.3 / 249.9 / 1772.0 |
| DFJ + our online optimizations (earlier) | local box | 50/50 | 1.373e9 / 1.045e9 | 868.0 / 683.1 / 287.2 / 2354.0 |
| ours (chi, no u_5 guess), attack_full | the benchmark node | 50/50 | 2.146e7 / 1.633e7 | 38.39 / 35.80 / 20.80 / 73.28 |
| ours (earlier) | local box | 50/50 | 2.146e7 / 1.633e7 | 53.64 / 50.50 / 28.70 / 113.80 |

Published-online aggregates (node): recovered 50/50; lookups mean 1.373e9, median
1.045e9 (total 6.86658e10); Bloom hits/key mean 17,503 (candidate-level, identical to
attack_dfj) -- raw probe hits mean 17,509 (the no-early-exit tail adds 359 probes in
total); key-schedule tries/key mean 7.169e7; wall s/key mean 757.3, median 571.0,
min 274.8 (seed 346), max 1927.6 (seed 325).

Same-machine ratio published / optimized online = 757.3 / 680.0 = **1.11x** (per-key
lookups, hits and key-schedule tries identical between the two by construction; the
only extra work is the re-peel of the 64-element column inside the s loop, small next to
the 64 SINV + multiset-hash lookups per s). Machines differ (node idle vs. local box at
load >190/192), so cross-machine wall-clock ratios should use the node columns; the
local-box figures (868 / 53.64) are from the earlier report.

## 5. Exact commands (`export OMP_NUM_THREADS=80`)

    gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_dfjpub sr7226_dfjpub.c -lm
    gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_gray   sr7226_gray.c   -lm

    ./sr7226_gray graycheck 32 1 32                                   # slice identity
    ./sr7226_dfjpub buildcold bloom_dfj_cold_k7.bin  41 0 4032 80     # published table (timed)
    ./sr7226_dfjpub builddfj  bloom_dfj_hoist_k7.bin 41 0 4032 80     # optimized table (timed)
    cmp bloom_dfj_cold_k7.bin bloom_dfj_hoist_k7.bin                  # full-table identity
    ./sr7226_dfjpub attack_dfj_pub dummy 41 seed=301 spot=8192        # gate (b), also 302,303
    ./sr7226_dfjpub attack_dfj_pub bloom_dfj_cold_k7.bin 41 seed=601  # gate (a), also 602,603
    cat bloom_dfj_cold_k7.bin > /dev/null
    printf '%s\n' $(seq 301 350) | xargs -P 25 -I{} bash -c \
      './sr7226_dfjpub attack_dfj_pub bloom_dfj_cold_k7.bin 41 seed={} > sweep_pub/key_{}.log 2>&1'
    # same-machine references:
    printf '%s\n' $(seq 301 350) | xargs -P 25 -I{} bash -c \
      './sr7226_dfjpub attack_dfj bloom_dfj_cold_k7.bin 41 seed={} > sweep_dfjopt/key_{}.log 2>&1'
    ./sr7226_dfjpub build bloom_ours_k7.bin 41 0 4032 80             # ours table (timed)
    printf '%s\n' $(seq 301 350) | xargs -P 25 -I{} bash -c \
      './sr7226_dfjpub attack_full bloom_ours_k7.bin 41 seed={} > sweep_ours/key_{}.log 2>&1'

Driver: `dfjpub_driver.sh` (setsid nohup); aggregation: `aggregate.sh <sweep_dir>`.
Node disk after run: cold table `bloom_dfj_cold_k7.bin` (256 GB) kept at
the run's working directory; hoisted and rebuild files deleted after use; the
separately built exact `tbl_dfj` table was not touched.
