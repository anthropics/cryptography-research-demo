# SR(7,2,2,6) exact-fingerprint-table experiment - re-run (2026-07-18)

Exact-table re-run. No Bloom filter: every table hit
is a true 64-bit fingerprint (ours: chi* fp) / multiset-hash (DFJ) collision.
Node: the benchmark node (208c / 1.8 TB RAM, shared; we used <=100 cores, peak ~1.25 TB disk).
Sources: sr7226_exact2.c and exact_driver.sh (this directory). The sweep
directories, build/sort logs and results summary are produced by the driver and
not bundled.

## Code
sr7226_exact2.c = the exact-table modes (buildexact / attack_exact / attack_exact_dfj)
merged with the CORRECTED canonicalizer of sr7226_canon2.c (only merge conflict: a
duplicate cmp_u64, of which the exact-branch copy was dropped).
Build: gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_exact2 sr7226_exact2.c -lm -lpthread
Checks: ./sr7226_exact2 canoncheck 2000000 50000 40000000 -> VERIFY: ALL PASS.
Slice sanity wu 2000..2001, OLD binary vs NEW: identical d-sequence count 34,748,288 and
69,496,576 staged fingerprints; 99.25% of unique fp values shared (the ~0.75% difference =
corrected even-parity S1=0 fp1 values, ~1/2 of the 1/64 S1=0 tier, as expected). Full build:
71,282,781,952 d-seq and cumulative alpha-brute count 1,877,282 byte-identical to the previous
run / canonfix build log.

## Pipeline (exact_driver.sh, detached via setsid nohup, sequential)
1. ours: buildexact tbl_chi 0 4032 100 -> 1.14 TB staged in 499 s (142,565,563,904 fp)
   -> sort_buckets tbl_chi 24 -> 1.12 TB sorted -> warm (cat into page cache)
   -> 50 keys seeds 301..350, xargs -P 50, taskset -c 0-99, attack_exact tbl_chi seed=N
   -> 50/50 rc=0 (07:15-07:27 UTC).
2. Disk: after the ours sweep the node had 1,360 GB free and tbl_chi was 1,120 GB; the DFJ
   stage (~570 GB staged + sort transient) did not fit within our ~1.3 TB disk cap alongside tbl_chi,
   so tbl_chi was deleted before the DFJ build (its bucket listing was recorded first). tbl_dfj
   (532 GB) remained on the node.
3. DFJ: buildexact tbl_dfj 0 4032 100 dfj (183 s, 71,282,781,952 fp) -> sort -> warm ->
   50 keys attack_exact_dfj seed=301..350, <=50 parallel -> 50/50 rc=0 (07:46-11:04; key 325,
   7 structures / 3.7e9 lookups, is the 11,864 s straggler, same straggler as the old run).
RAM: build staging ~13 GB; sort <= ~110 GB (24 buckets in RAM); warm = table size page cache.

## Results
| | ours (canon2 chi fp) | DFJ (multiset hash) | previous run |
|---|---|---|---|
| keys recovered | 50/50 | 50/50 | ours 150/150, DFJ 50/50 |
| lookups/key mean (median) | 2.146e7 (1.633e7) | 1.373e9 (1.045e9) | -- |
| true-FP hits/key mean (median; range) | 3,741 (2,874; 1,318-10,064) | 6.06 (6); hits incl. true hit 7.06 | ours ~3,396; DFJ ~7 |
| FP hits / lookup ("per query") | 1.7434e-4 (mean of ratios 1.7431e-4) | 4.41e-9 | ours 1.75e-4 |
| FP per probe | 8.72e-5 (2 probes/lookup) | 4.41e-9 (1 probe) | -- |
| duplicate table entries (sort logs) | 1.7481% | 1.7438% | ~1.75% both |
| wall s/key mean / median / min / max | 217.3 / 191.9 / 120.9 / 435.4 | 4,621 / 3,740 / 1,702 / 11,864 | old node, Bloom-era |
| build time (staged) | 499 s (2^37.05 fp) | 183 s | -- |

Reproduction: the per-query FP rate and duplicate aggregates reproduce (ours 1.743e-4 vs
1.75e-4 per query; dup 1.748%/1.744% vs ~1.75%; DFJ ~7 hits/key; 100% key recovery). The ours
mean true-FP hits/key (3,741) exceeds the earlier quoted ~3,396 because this is a different
50-key subset (seeds 301..350, mean 2.146e7 lookups/key - the same counts as the canonfix 50-key
sweep) vs the old 150-key aggregate (mean lookups ~1.9e7/key): FP count scales with lookups, the
per-query rate is the invariant. Also the corrected canonicalizer changes ~0.75% of fp1 values
(even-parity S1=0 entries), so per-key hit counts would not match the earlier run
digit-for-digit anyway. Wall-clock is NOT comparable to the old Bloom runs: lookups here are binary searches in
mmap'd sorted buckets (1.1 TB / 0.5 TB), 50 keys concurrently on a 100-core cpuset of a shared node.
