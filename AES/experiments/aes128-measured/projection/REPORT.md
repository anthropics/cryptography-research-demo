# Full measured three-way comparison (offline + online) and full-scale projection
## 7-round AES-128: (A) OURS vs (B) DFJ-OPTIMIZED vs (C) DFJ-AS-PUBLISHED

What is **measured** here: per-table-entry offline build cost (real AES-128,
real Prop-2 bases, real 255-element d-sequences), per-(k_{-1}[diag],
u_6[anti-diag])-candidate online cost (real 7-round ciphertexts), and two
encryption denominators, all on the same machine.  What is **counted**: the
exponents (2^88 / 2^80 table entries, 2^88 online candidates, 2^96 inner
s-iterations) and the 160-lookup ledger.  **Nothing here runs 2^80..2^96
work** — every total below is (measured constant) x (counted exponent).

## Machine, load, exact commands

- the benchmark node: 2x Intel Xeon Platinum 8481C
  (Sapphire Rapids, GFNI/AVX-512), 104 cores / 208 threads, 1.8 TB RAM,
  TSC 2.70 GHz; gcc 11.4 `-O3 -march=native`.  No PMU in the VM: core cycles
  = `rdtscp` ticks x imul-latency core/TSC calibration performed around every
  configuration (core clock ~2.9 GHz, ratio 1.07-1.11; TSC ticks also quoted).
- Concurrent load: other work occupied 56-110 cores during the
  window (loadavg 53 -> 97 -> 28).  Our runs were pinned with `taskset` to
  core 94 (offline + denominators, node1) and core 36 (online, node0); their
  HT siblings (198, 140) were verified idle the whole time by a 10-s
  `/proc/stat` sampler.  5 reps per configuration,
  medians + [min,max]; spreads within ~1-2% except where noted.  Shared
  L3/memory contention inflates absolute online cycles ~10-20% vs the earlier
  IDLE-node run of the same binary (see `../online/REPORT.md`: A-allvec 581,
  A-plain 1916, B-plain 165103, B-allvec 113813 cyc); all ratios are
  self-consistent in this run and re-confirmed by the clock-neutral
  interleaved head-to-head.
- Sources `offline3.c` (offline
  three-way), `online3.c` (= online_bench.c + published-row breakdown +
  per-phase summary), `denom.c` (T-table + AES-NI); launcher `run_all.sh`.
  ```
  gcc -O3 -march=native -o offline3 offline3.c -lm     # headers via -I../common (see Makefile)
  gcc -O3 -march=native -o online3  online3.c -lm
  gcc -O3 -march=native -o denom    denom.c
  taskset -c 94 ./offline3 gate
  taskset -c 94 ./offline3 bench {A|B|C} 16 5 20 {2|1}  # 16 bases x 2^20 = 2^24 entries, 5 reps; stage 2=full, 1=gen-only
  taskset -c 36 ./online3 0xA15BE2C4 20 16 5 18          # 2^20 A cands, 2^16 B cands (2^24 s-iters), 2^18 Bpub cands (2^26 s-iters), 5 reps
  taskset -c 94 ./denom 100000000 5                        # 10^8 blocks each, 5 reps
  ```

## Correctness gates (verbatim) — ALL PASS

Offline (`offline3 gate`, 16,384 entries = 2 bases x 2^13 Gray states, > 10^4 asked):
```
[gate o-1] delta-sequence A(walk)==C(cold) bytewise        : 16384/16384
[gate o-2] raw-sequence hash A==B==C (FNV64 of d_1..d_255): 16384/16384
[gate o-3] B multiset-hash(walk seq) == hash of C cold seq : 16384/16384
[gate o-4] A chi* fp: kernel_v2_vec == kernel_full_v2 (both words, bit-exact): 16384/16384
=== gate: ALL PASS ===
```
Plus a full-stream cross-check: over the 2^24-entry timing streams the B and
C accumulators (XOR of all multiset hashes) are **identical**
(`827953b46068221f`): B's incremental walk and C's independent cold trace
produced the same delta-sequence for every one of 16.8 M entries.
Online gates (i-a .. iii, including the published-kernel gates i-c/i-d):
ALL PASS (`logs/online3.log`).

## 1. OFFLINE per-entry cost, all three (2^24 real entries each, one pinned core)

Entry = one genuine 255-element sequence d_w = Delta x_5[0] of a real Prop-2
parameter base (16 random valid 10-byte bases x all 2^20 DDT branch states).
Counters are scalar-**equivalent** lookups (one GFNI 64-lane
`VGF2P8AFFINEINV` = 64 S-box-class evaluations), ddt_gray-style.

| builder | COUNTED S-box/element | COUNTED other 256-table reads/element | ticks/entry median [min,max] | **core-cyc/entry** | cyc/element | split | entries/s/core |
|---|--:|--:|--:|--:|--:|---|--:|
| **(A) OURS**: DDT-Gray walk -> chi* (fp0,fp1), `kernel_v2_vec4` | **1.0664** (= (15*1+1*2)/16 + cold start/2^20) | 0 (fingerprint is bit-sliced GFNI, not a lookup table) | **347.4** [345.0,349.4] | **375.7** | 1.473 | walk+stage 49.8 cyc (walk alone ~29); chi* kernel ~326 = 1.28/elt | 7.77 M |
| **(B) DFJ-OPT**: same walk -> multiset hash of the raw sequence | **1.0664** | **1.00** (R64 gather-sum) | **264.9** [262.3,267.0] | **286.5** | 1.124 | walk 29.3; hash ~257 = 1.01/elt | 10.19 M |
| **(C) DFJ-PUB**: cold per-entry trace (rounds 2-4) + multiset | **24.094** (24 SB/elt + 24/entry; paper's loose "~20") | **39.0** GF-mul reads + **1.0** hash | **10554.6** [10519.4,10642.4] | **11493.3** | 45.07 | cold trace ~11.2 k (~98%); hash ~0.2 k | 0.256 M |

Per entry: C/A = 30.6x, C/B = 40.1x, A/B = 1.31x (the chi* canonicalization
costs ~1.3 cyc/element more than DFJ's multiset hash on the same walk).  (A)
reproduces the earlier table-build measurement (345-400 cyc/entry).

## 2. ONLINE per-candidate cost, all three (unit = one (k_{-1}[diag], u_6[anti-diag]) candidate)

For (B)/(C) the unit contains the whole 256-value u_5[0]=s loop; for (A) it is
one chi-star fingerprint + 2 probes (the Mobius bridge absorbs u_5[0]).
Verbatim three-way table (core cycles, median of 5; [min,max] ticks):

| variant | core-cyc/cand | [min,max] ticks | vs ours-allvec | vs ours-plain |
|---|--:|---|--:|--:|
| (C) DFJ-published (plain scalar, peel recomputed inside every s-iteration) | **371,236** | [337920, 358833] | 569x | 173x |
| (B) DFJ-optimized plain (shared peel) | 199,817 | [182657, 192611] | 308x | 93.8x |
| (B) DFJ-optimized gray (shared peel + Gray) | 199,025 | [180663, 184106] | 304x | 92.5x |
| (B) DFJ-optimized allvec (GFNI peel + per-s GFNI iSBOX) | **117,141** | [107647, 109626] | 179x | 54.3x |
| (A) ours plain (scalar fused peel, vector chi* kernel) | 2,150 | [1998, 2024] | 3.29x | 1.00x |
| (A) ours gray | 1,819 | [1685, 1708] | 2.78x | 0.84x |
| (A) ours allvec (GFNI peel + GFNI L^-1 + vector chi*) | **649** | [593, 648] | 1.00x | 0.30x |

Per-phase breakdown (core cycles per candidate; cycles per sequence element
in brackets — A has 255 elements/candidate, B/C have 256x255):

- **(A) plain = 2,150**: peel 1,004 [3.94] | sequence build (XOR+L^-1) 419
  [1.64] | chi* fingerprint incl. GFNI field inverse 466 [1.83] | 2 probes
  262.  **(A) allvec = 649**: GFNI peel 58 [0.23] | GFNI L^-1 + vector kernel
  407 | 2 probes 184.
- **(B) plain = 199,817**: shared peel 1,016 once | s-loop 256 x [iSBOX pass
  + R64 multiset hash] = 142,605 = 557 cyc per s [2.18/elt; the hash rides the
  iSBOX dependent-load chain, ~0 marginal] | 256 probes = 56,196 = 220/s ->
  777 cyc per s-iteration.  **(B) allvec = 117,141** = 458 cyc/s (GFNI peel
  + GFNI per-s iSBOX pass; R64 sum stays a gather).
- **(C) published = 371,236**: 256 x [full 255-element 4-iSBOX peel redone]
  = 241,049 = 942 cyc per s [3.69/elt], plus s-step iSBOX + hash + probe
  remainder 130,187 = 509/s -> 1450 cyc per s-iteration [5.69/elt].  (The
  standalone peel+iSBOX stage-2 cut of this row compiles to a slower
  byte-checksum loop than the full row, so (C)'s split is reported as
  peel-only + remainder rather than four cumulative cuts.)

Clock-neutral interleaved head-to-head (21 rounds, median [min,max]):
optimized-vs-optimized **allvec 277x [224,341]**, plain 68.6x, gray 95.3x;
published-vs-ours: 149x [132,170] vs ours-plain, **689x [510,777]** vs
ours-allvec.

## 3. Denominators (same node, 10^8 blocks each, pinned core 94, 5 reps)

Both implementations verified against the byte-wise reference on 10^4 random
blocks (T-table PASS, AES-NI PASS).

| denominator | ticks/block median [min,max] | core-cyc per AES-128 encryption |
|---|--:|--:|
| software T-table AES-128 (4x1KB Te, independent counter blocks) | 182.16 [180.28, 188.72] | **195.7** |
| AES-NI, 8-wide pipelined throughput | 5.241 [5.232, 5.258] | **5.63** |
| (AES-NI single-stream latency, for reference) | 31.34 [31.29, 31.54] | 33.7 |

## 4. PROJECTION at D = 2^105 (measured constants x counted exponents)

Counted exponents: (A) N_off = 2^88 entries (both parity words stored, ~26
B/entry -> M ~ 2^92.7 B), N_on = 2^88 candidates; (B),(C) N_off = 2^80
entries (DFJ's 10-byte parameterization; B ~16 B/entry -> 2^84 B, C ~64 B
multiset -> 2^86 B), N_on = 2^88 (k_{-1},u_6) candidates each containing
the 256-value s-loop = 2^96 inner iterations.  Each variant at its natural
implementation level (A: allvec online, walk+chi* offline; B: allvec; C:
plain published):

| variant | offline cyc | online cyc | **total cyc** | in T-table enc (/195.7) | in AES-NI enc (/5.63) | M (bytes) | D |
|---|--:|--:|--:|--:|--:|--:|--:|
| **(A) OURS** | 2^88 x 375.7 = 2^96.55 | 2^88 x 649 = 2^97.34 | **2^98.00** | 2^90.39 | 2^95.51 | 2^92.7 | 2^105 |
| **(B) DFJ-OPTIMIZED** | 2^80 x 286.5 = 2^88.16 | 2^88 x 117141 = 2^104.84 | **2^104.84** | 2^97.23 | 2^102.34 | 2^84 | 2^105 |
| **(C) DFJ-AS-PUBLISHED** | 2^80 x 11493 = 2^93.49 | 2^88 x 371236 = 2^106.50 | **2^106.50** | 2^98.89 | 2^104.01 | 2^86 | 2^105 |

**Measured-cycle ratios: total(B)/total(A) = 2^6.84 = 114x;
total(C)/total(A) = 2^8.50 = 363x.**  (Like-for-like scalar/plain-peel
online rows: A 2^99.30, B 2^105.61, C 2^106.50; B/A = 2^6.31 = 79x,
C/A = 2^7.20 = 147x.)  Offline is negligible for B and C (2^-12..2^-13 of
their totals) but is **2^-1.45 of A's total** — A's table has 2^8 more
entries and a heavier per-entry reduction, so at full scale its build sits
within ~0.8 bit of its online phase.

The same totals in the paper's **160-lookup units** (counted lookups, C_AES = 160):

| ledger row | l_off | l_on | T_off | T_on | total |
|---|--:|--:|--:|--:|--:|
| (A) as implemented here (DDT-Gray only; online 4-iSBOX peel, no W-cache) | 535 | 1538 | 2^89.74 | 2^91.26 | 2^91.70 |
| (A) paper-final (UU-cache offline / W-cache + chi*' online) | 263 | 388 | 2^88.72 | 2^89.28 | 2^90.02 |
| (B) DFJ-optimized (shared peel; hash charged 1/element) | 527 | 131,836 | 2^81.72 | 2^97.69 | 2^97.69 |
| (C) DFJ-as-published (per-s peel; hash charged 1/element) | 6,399 | 391,936 | 2^85.32 | 2^99.26 | 2^99.26 |
| (C) DFJ'13's own accounting (1280/s-iter, hash free, ~20 SB cold) | 5,355 | 327,680 | 2^85.06 | **2^99.00** | 2^99.00 |

**Counted-ledger ratios:** published T / paper-final T = 2^99 / 2^89.28 =
**2^9.72** (T-only; 2^8.98 with both offline builds included); structural
iteration-count ratio 2^96 / 2^88 = **2^8** exactly.

## 5. Does the measurement match the ledger?

- The structural **2^8** is fully realized in cycles: same-implementation-
  level online ratios are 69-93x plain and 179-277x with both sides
  GFNI-vectorized, and the published kernel's 2^96-vs-2^88 count is literal
  (its cost is s-count x constant: 371,236 / 256 = 1,450 cyc per s-iteration).
- Total (C)/(A) in measured cycles is **2^8.50** (natural levels) vs the
  ledger's **2^9.72** (T-only) / 2^8.98 (with offline): ~0.5-1.2 bit below
  the ledger.  The gap is explained by two measured-vs-counted constants:
  (i) our online row as implemented has **no W-cache** (ledger l_on = 1538,
  not the paper-final 388; the GFNI allvec kernel nevertheless runs it in
  649 cyc = 0.42 cyc per counted lookup), and (ii) the published kernel
  realizes ~0.95 cyc per counted lookup (371,236 cyc vs 391,936 counted) —
  the ledger prices DFJ's scalar lookups almost exactly at 1 cycle while
  pricing our GFNI-vectorized fingerprint above its wall-clock cost.  Against
  the as-implemented-A ledger (2^91.70) the counted ratio is 2^7.30, which
  brackets the measured 2^7.2 (plain-vs-plain) .. 2^8.5 (natural-level).
- Offline: measured A = 375.7 cyc/entry vs ledger l_off = 535 -> 0.70
  cyc/lookup (the walk runs at 0.19 cyc/element via GFNI); measured C =
  11,493 cyc vs 6,399 counted -> 1.8 cyc/lookup (the ~39 GF-mul reads per
  element are real work the S-box-only convention does not charge).  The
  counted walk rate **1.0664 S-box/element** matches the paper's
  (15x1+1x2)/16 = 1.0625 (+ the 2^-20-amortized cold start); C's counted
  cold rate is **24.09 S-box/element** (the paper's "~20" is the loose figure).
- Convention check on the denominators: one T-table encryption = 195.7 cyc
  ~ 196 one-cycle lookups, so the 160-lookup unit (A: 2^98.00-7.32 = 2^90.68)
  and the T-table unit (2^90.39) agree to 0.3 bit; the AES-NI unit is 5.1
  bits more demanding (2^95.51) — the known ~2^5 convention gap.

## Assumptions

1. Per-unit constants come from 2^24 offline entries and 2^16-2^20 online
   candidates on real 7-round AES-128; projection to 2^80..2^96 units assumes
   constant per-unit cost.  No memory-hierarchy change is modeled for the
   astronomically larger real tables: a real 2^88/2^80-entry table makes
   probes far slower than the L2-resident dummy set used here, which
   penalizes B/C (256 probes per candidate) more than A (2 probes), i.e. the
   direction is conservative for A.
2. D = 2^105 for all three (DFJ'13's data requirement; the bridge does not
   change D); data-collection cost is excluded from all totals by convention.
3. Exponents: A N_off = N_on = 2^88 (the Mobius bridge absorbs u_5[0]);
   B/C N_off = 2^80 (DFJ's 10-byte table), N_on = 2^88 (k_{-1},u_6)
   candidates x 256 s = 2^96 inner iterations.
4. Table memory: A 26 B/entry (two ~12-byte fingerprints + alignment) ->
   2^92.7 B; B 16 B/entry (8-byte multiset hash + payload) -> 2^84 B;
   C 64 B/entry (DFJ multiset, ~512 bits) -> 2^86 B.
5. Offline builds are priced as single-core per-entry costs; parallel builds
   scale linearly in cores (the real build is write-bandwidth bound, see
   `../offline/REPORT.md`, which only strengthens the conclusion).
6. Core cycles derive from TSC ticks via imul-latency calibration (no
   hardware PMU available); denominators measured by the same method on the same node, so
   all unit conversions are same-machine.  Concurrent load inflated absolute
   online cycles ~10-20% vs the idle-node run; ratios unchanged (interleaved
   head-to-head confirms).
7. (C) is DFJ'13's algorithm implemented plainly (scalar, peel recomputed per
   s, no caches, no early exit); their paper's own lookup accounting (1280
   lookups per (k_{-1},u_6,u_5[0]) iteration, ~20-S-box cold build, multiset /
   hash free) is the 2^99 row.

## Sources (this directory)

`offline3.c`, `online3.c`, `denom.c`, `project.py` (projection arithmetic),
`run_all.sh`.  The run writes `projection_out.txt` and per-stage logs, which are
not bundled.
