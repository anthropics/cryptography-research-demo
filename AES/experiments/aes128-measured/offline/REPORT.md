# Offline table build — partial REAL table for the 7-round AES-128 attack

Partial table: **2^36 genuine entries** (of the paper's 2^88), staged into 256 bucket
files, sorted, de-duplicated and analysed; **2^34 genuine online wrong-candidate
probes**; **128 planted known answers**.  Everything below was produced by the
programs in this directory (`table_build`, `analyze`, `classify`), no simulation.

## Machine

- 2 x Intel Xeon Platinum 8481C @ 2.70 GHz (Sapphire Rapids), 104 physical cores /
  208 vCPU (SMT-2), 1842 GiB RAM, fast local NVMe scratch (~6 TB), Linux 6.6.97,
  gcc 11.4, AVX-512 + GFNI.
- 190 worker threads used (two threads per physical core over 95 cores); TSC = 2.700 GHz
  (reference clock; all per-entry cycle counts below are TSC ticks).
- Wall clock for the whole pipeline: 2026-07-18 06:28:24 - 07:11:43 UTC (43 min).

## Exact commands

```
make                                     # gcc -O3 -march=native -pthread (table_build, analyze, classify)
# pre-run bit-exactness verification (dense check strides), 2^24 entries:
./table_build build nb=16 threads=32 nka=8 chkkern=253 chkcold=997 dir=.../test_buckets
# compute-only throughput:
./table_build bench threads=190
# the real run (run_all.sh):
./table_build build nb=65536 threads=190 nka=128 dir=./buckets
./analyze all threads=190 sortpar=48 pinst=262144 pg=65536 dir=./buckets
#   (= sort/dedup/stats, fp1-word stats, 2^18 x 2^16 = 2^34 probes, known-answer)
shuf -n 100 coll0.txt -> t0.txt ; shuf -n 100 coll1.txt -> t1.txt   (dups.txt was empty)
./table_build resolve nb=65536 threads=190 nka=128 targets0=t0.txt targets1=t1.txt > resolve.txt
./classify < resolve.txt > classify.txt
```

What an entry is: a random DFJ'13 Prop-2-style 10-byte parameter base
(Din=Delta z_1[0], x_2[col0], Dout=Delta w_4[0], z_4[col0]) is rejection-sampled until
all 20 derived DDT constraints are satisfiable; the 20 DDT branch bits (4 x_4-diagonal +
16 x_3) are then Gray-walked, giving 2^20 genuine 24-reference states per base; each
state yields the 255-element difference sequence d and its brute-free chi-star
fingerprint (fp0, fp1) from the kernel_v2 AVX-512 kernel (the same headers as the online
bench).  65536 bases x 2^20 = 2^36 entries, 16-byte records (fp0,fp1), staged by the
top byte of fp0 into 256 buckets (1.1 TB), each bucket then sorted in RAM by (fp0,fp1)
and rewritten.

## Sizes actually built

| item | value |
|---|---|
| table entries | 68,719,476,736 = 2^36.000, + 128 planted known answers = 68,719,476,864 records |
| parameter bases | 65,536 (11,233 contain at least one DDT n=4 byte; the walk takes 2 of the 4 branches there) |
| on disk | 256 buckets, 1.00 TiB, all sorted and scanned for duplicates |
| probes | 17,179,869,184 = 2^34 (262,144 real AES-128 instances x 65,536 wrong k7a guesses) |
| known answers | 128 real 7-round AES-128 instances, true (x_3[col0],x_4,x_5[diag]) tuples forced in and disclosed in ka.bin |

Bit-exactness checks (all passed, counts printed by the programs):
start-up self-test 248/248 (vector walk vs cold scalar walk vs independent prop2 path,
vector kernel vs scalar reference kernel); pre-run dense build: 66,368/66,368 vector
kernel == scalar reference kernel, 16,832/16,832 walked sequences == cold_E ==
slow reference `prop2_all`; inside the real 2^36 build: 1,048,576/1,048,576 kernel
vec==scalar and 8,454,144/8,454,144 sequences vs cold_E and vs prop2_all; classify
recomputed (fp0,fp1) from 400 emitted sequences: 400/400 agree with the stored records;
records written 68,719,476,864 = expected.

## Throughput

| measurement | rate | per entry |
|---|---|---|
| 1 thread, compute only (walk + kernel) | 5.66 M entries/s/core | 405.6 TSC ticks (walk 102.9 + kernel 302.7) |
| 190 threads, compute only | 656.5 M/s aggregate; 3.46 M/s/thread (SMT: ~6.9 M/s per physical core) | 647.9 TSC ticks/entry |
| 190 threads, full build incl. 1 TiB staging I/O + checks | 286.6 M/s aggregate; 2^36 in 239.8 s + 79.9 s fsync (I/O 58% of thread time, 4.27 GiB/s sustained) | 688.9 TSC ticks compute + I/O |
| sort + dedup + stats of 1 TiB (48 buckets in parallel) | 835 s | |
| probe generation (real delta-set + honest peel + kernel) | 124.2 M probes/s aggregate, 0.654 M/s/thread | ~4131 TSC ticks/probe |

At the reference clock (2.7 GHz) a single core spends ~406 cycles per table entry
(~303 of them in the fingerprint kernel); the 2^36-entry compute is ~105 s of this
208-vCPU node, the 1 TiB write is what dominates the 240 s build.

Extrapolation to the paper's 2^88 entries (compute only, generously assuming the I/O
and the 2^92-byte store were free): 2^88 / 5.66e6 entries-per-core-second = 5.5e19
core-seconds = 1.7e12 core-years, or 2^88 / 6.57e8 = 4.7e17 node-seconds = 1.5e10
node-years on this machine; storage would be 2^88 x 16 B = 5e27 bytes.  The partial
table is 2^52 times smaller than the full one; nothing here shortens the 2^88 build, it
only measures the per-entry constant (~406 cycles) honestly.

## Duplicate / collision analysis (N = 2^36 + 128 records)

- **Full-record (128-bit pair) duplicates: 0.**  All 68,719,476,864 (fp0,fp1) pairs are
  distinct.  In particular no two distinct parameter tuples (base, Gray index) produced
  the same d-sequence: such a coincidence necessarily shows up as a duplicate record and
  there are none at N = 2^36.
- **fp0 64-bit word collisions** (equal fp0, different fp1): **135**, ideal uniform
  expectation C(N,2)/2^64 = 128.0 (Poisson sigma ~11.3); empirical pairwise rate
  5.72e-20 vs 2^-64 = 5.42e-20.
- **fp1 64-bit word collisions: 177.**  fp1 is the minimum over candidate canonical
  hashes, so its top byte is non-uniform (q=0x00 holds 535,844,855 records, q=0xff only
  1,400,540); the expectation given that observed skew and uniformity inside each
  top-byte class is 170.6 (pure-uniform would give 128.0).  Observed 177 matches the
  skew-adjusted model.
- **Bucket sizes** (top byte of fp0): min 268,394,512, max 268,478,139, mean 2^28;
  chi-square = 275.0 on 255 df (z = 0.88) -> consistent with uniform.  Top-20-bit
  prefix histogram: chi-square 1,049,610 on 1,048,575 df (z = 0.71) -> uniform.
- **100+100 collision sample, provenance-resolved and classified** (`table_build
  resolve` regenerates all 2^36 entries and emits the full 255-byte sequence +
  parameter tuple for every entry whose fp0 is in the 100-value sample t0.txt or whose
  fp1 is in t1.txt; `classify` then compares the two sequences of every colliding pair):
  400 entries / 200 pairs.  Same provenance: 0.  Identical sequences from distinct
  tuples: 0.  Same multiset / same parity set chi(D) / AGL(1,256)-equivalent canonical
  sets (the structural classes the fingerprint identifies by design): 0.  **200/200 =
  genuine 64-bit FNV hash collisions of distinct canonical sets** — exactly the
  mechanism, and the count, an ideal 64-bit word predicts.
- What N = 2^36 can and cannot bound: the 64-bit word statistics are measured to
  within a few percent (hundreds of events, chi-square and prefix counts on the
  uniform line).  For the 96-/128-bit behaviour the sample is silent: C(N,2) = 2^71
  offline pairs give an ideal 96-bit fingerprint only 2^-25 expected collisions (128-bit:
  2^-57); observing zero bounds the pairwise collision probability only to
  p < 3/2^71 ~ 2^-69.4 (95% rule of three).  The extrapolation to 2^-96 / 2^-128 rests
  on (a) both 64-bit words behaving ideally and (b) the absence of any structural
  (non-hash) collision class in the resolved sample — both hold here.

## Wrong-candidate probe test (2^34 genuine online probes)

2^18 real AES-128/7-round instances (random key, random plaintext, genuine 256-text
delta-set from the true k_0 diagonal) x 2^16 random wrong k7[anti-diag] guesses each,
honest peel -> eH -> kernel_v2, probed against the sorted table (fp0 against the
bucket files, fp1 against the in-RAM fp1 partitions).

- **128-bit joint hits (a probe equals a table entry in both words): 0**
  (model |T|x|P|/2^128 = 3.5e-18; 96-bit model 1.5e-8).
- **Exact 64-bit word matches: fp0 side 63, fp1 side 65**; ideal expectation
  |T|x|P|/2^64 = 64.0 each (fp1-side skew model: 85.3).  All 63+65 resolved to provenance
  (probe instance + guess vs table record): never both words of the same entry, so under
  the attack's EITHER-word acceptance rule 2^34 wrong candidates produce 128 single-word
  hits = the chance-word-match floor, and 0 full fingerprint hits.
- Partial-match (common-prefix) statistics, observed/expected ratio: fp0 side
  >=32 bits 1.000, >=40 1.000, >=48 1.000, >=56 1.011, =64 0.984 (uniform model);
  fp1 side 1.000 at >=32/40/48/56 bits and 0.762 at =64 against the top-byte-skew
  model (1.333 against a naive uniform model — the skew is real and fully explained by
  the min-of-hashes construction).
- Within the probe set itself 131,114 fH0 (131,116 fH1) values repeat; this is entirely
  the expected number of duplicated random 32-bit k7a guesses, 2^18 x C(2^16,2)/2^32 =
  131,072 — not a fingerprint artefact.

## Known-answer coverage

128 planted instances: real AES-128 reduced to 7 rounds, random key and plaintext,
true reference bytes (x_3[col0], x_4, x_5[diag]) -> Prop-2 sequence -> (fD0,fD1) written
as an ordinary record (disclosed with keys/plaintexts/positions in ka.bin, logs/
analyze.out).  Independently, the attacker-side online computation with the TRUE k_0
diagonal (delta-set) and TRUE k7 anti-diagonal (peel) gives (fH0,fH1).

- present in the sorted table by binary search: **128/128**
- genuine online fingerprint matches its own entry: fp0-word 64, fp1-word 64
  (one word always coincides — the raw-parity vs Add-0-parity case of the bridge
  identity; bridge_ok >= 253/255 for all 128, 255/255 in every instance listed), EITHER-word **128/128**
- **recovered (present AND online match) = 128/128**; 64 found directly by the fp0 index
  lookup, the other 64 are the fp1-word case (the real attack indexes both words).

## What this does and does not establish

Does: the offline entries are real (Prop-2 rebound bases, full DDT-Gray expansion,
255-element delta-sequences verified against two independent reference implementations,
fingerprints bit-exact to the scalar reference kernel); the per-entry generation cost is
measured (~406 reference cycles per core, ~303 in the kernel) and scales across 190
threads; at 2^36 real entries the fingerprint distribution is statistically ideal
(uniform buckets and prefixes, 64-bit word collision counts on the 2^-64 line, zero
duplicate sequences, and every sampled word collision resolved to a genuine hash
collision of distinct canonical sets); 2^34 genuine online wrong candidates behave as
independent uniform draws against the table (64-bit word matches 63/65 vs 64 expected,
no fingerprint hit); and true candidates, planted honestly, are always found (128/128).

Does not: this is 2^36 of the paper's 2^88 entries (a 2^-52 fraction); it cannot
measure 96- or 128-bit collision rates directly (zero events is the trivial outcome),
cannot exhibit the real attack's wrong-candidate hit count (which would need the full
table and ~2^72 candidates), says nothing about the data/time/memory complexities of
the attack as a whole, and the 1.5e10 node-year / 5e27-byte extrapolation is exactly
the reason the full table is not built.  The claim supported is the constant factor
and the hash-ideality assumption, not feasibility.

## Sources (this directory)

table_build.c, analyze.c, classify.c, the shared headers under ../common/, Makefile,
run_all.sh.  The run produces logs/ (build.log, analyze.out with all statistics + the
KA list, resolve.txt, classify.txt, the collision/duplicate lists) and the ~1 TB bucket
data under ./buckets; none of that output is bundled in this repository (see the
top-level AES README) - the pipeline regenerates it.
