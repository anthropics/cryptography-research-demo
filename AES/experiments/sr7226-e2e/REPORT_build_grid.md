# SR(7,2,2,6) gray-stack: AES-128 table-enumeration optimizations ported to the small-scale demo — measured numbers

Sources for this experiment: `sr7226_gray.c` (this directory).
Source: `sr7226_gray.c` = the DFJ table-attack driver + the additions below.
Original production builders and both 2^41-bit tables untouched (tables only read).

## 0. What was built

`sr7226_gray.c` adds, next to the untouched original `build_table()` (hoisted default / `COLD_MODE` naive):

* four instrumented build kernels over the same work-unit enumeration
  (`blk_cold` = naive, `blk_hoist` = the existing a-side hoist, `blk_gray` =
  **mixed-radix reflected Gray walk over the 4 DDT branch digits (b0,b3,b1,b2)**
  with incremental update of the 64-element d-sequence (one flipped digit ->
  only the x3->x4->x5 terms through that cell are recomputed: 2 S-box/element),
  `blk_gray2` = gray **plus** the two x_4-layer DDT-partner digits
  (z4_0<->z4_0^dz4_0, z4_1<->z4_1^dz4_1, 1 S-box/element, placed innermost) - the
  exact analogue of the AES 16 x_3-bits + 4 x_4-bits walk),
  via a loopless Knuth-7.2.1.1 mixed-radix Gray generator (radices n_j in {2,4} = DDT multiplicities);
* op counters (S-box / GF-constant-mult / XOR / inverse lookups), split into
  *inner* (inside the 64-element loops = c_omega accounting) and *overhead*
  (per entry/block), compiled in only with `-DINSTR`;
* zero-tolerance verification mode `graycheck`: runs the ORIGINAL builder
  (hoisted and `COLD_MODE`) plus the four new kernels on the same work units
  and proves emitted-multiset identity;
* online `gray` option in `attack_full`/`attack_dfj` (k_6[0]-column cache,
  only the DDT-constrained k_6[3] partner column recomputed per candidate),
  with per-candidate op counters;
* timing modes `build3`, `buildorig` (original code), and C_post-isolation
  emit modes (`chionly`, `msonly`, `none`).

Build:
    gcc -O3 -march=native -fopenmp -DBLOOM_K=7          -o sr7226_gray     sr7226_gray.c -lm   # timing binary
    gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -DINSTR  -o sr7226_gray_cnt sr7226_gray.c -lm   # counting binary

## 1. Hard gate - emitted multisets identical across all builders

`./sr7226_gray graycheck 12 1 12` (12 diverse work units (dy1,x2_0) spread
over dy1=1..63, x2_0 varied; every work unit enumerated in FULL - all x2_1,
dw4, z4_0, z4_1, branch tuples; 210,074,240 table entries per builder).
Six builders each: **orig-hoist, orig-cold (original code), new-cold,
new-hoist, new-gray, new-gray2**. Per work unit and builder we collect
(a) the 64-bit strong order-sensitive hash of every raw 64-element
d-sequence into an array, sort, and memcmp against the orig-hoist array,
(b) the entry count, (c) commutative sum/xor accumulators of the hashes AND
of the *emitted* fingerprints of both table variants (chi canonical
fingerprint, both parities = ours; DFJ multiset hash) computed from each
sequence.

Result (`gc12.log`): **ALL MULTISETS IDENTICAL - 12 work units x 6 builders,
sorted-array equal, counts equal, all accumulators equal, 0 discrepancies.**
(A 2-wu pre-run, `gc_small.err`, same result.) Per-wu entry counts
16.80M-18.29M. The chi / multiset-hash fingerprints are deterministic
functions of the d-sequence, so identity of the raw-sequence multisets
implies identity of both emitted Bloom-entry multisets; the fingerprint
accumulators confirm it directly (chi_sum/chi_xor/ms_sum/ms_xor equal for all
6 builders on every wu). Zero tolerance met; no discrepancy to report.

## 2. Offline build cost - measured grid

Single-threaded, `taskset` to one core, one configuration per core, <=4
concurrent; **box load 110-150 on 192 cores during all runs** (shared
machine) - spread across 3 reps <=7%. Long pass: wu 1000..1005 (~106M
entries/run, 45-130 s/run); the 1-wu cross-check pass agrees within 3%.
Bloom for timing = anonymous 2^30-bit filter (not the 2^41-bit production
file), so absolute rates exceed the production 80-thread build (0.32M
entries/s/thread ours, 0.64M DFJ); the cross-mode comparison is the point.
INV = the d_w -> 1/d_w lookup feeding the chi fingerprint (ours only; POWINV analogue).
MULT = GF(64) constant-multiplier (x2,x3 MixColumns coefficients) table lookups.

### 2a. entries/sec (mean of 3, long pass) and counted ops per table entry (64 elements/entry)

| build mode | ours (chi, 2 parities) entries/s | DFJ (multiset hash) entries/s | seq-gen only entries/s | SBOX/entry | MULT/entry | XOR/entry | INV/entry |
|---|---:|---:|---:|---:|---:|---:|---:|
| naive (cold)            | 0.921 M | 1.779 M | 3.30 M | 388.5 | 648.5 | 712.2 | 64 |
| hoisted (existing)      | 1.002 M | 2.058 M | 4.46 M | 239.4 | 448.0 | 511.7 | 64 |
| **DDT-Gray (4 x_3 digits)** | **1.043 M** | **2.211 M** | **5.47 M** | **143.6** | **156.9** | **404.0** | 64 |
| DDT-Gray2 (+2 x_4-layer digits) | 1.101 M | 2.446 M | 7.66 M | 83.9 | 87.2 | 293.8 | 64 |
| original code, hoisted  | 0.989 M | 1.893 M | - | (== hoisted) | | | |
| original code, COLD_MODE| 0.844 M | 1.426 M | - | (== naive) | | | |

(new-hoist reproduces the production code's numbers: 1.00M vs 0.99M chi, 2.06M vs 1.89M DFJ.)

### 2b. per-(entry, element) op counts = measured c_omega^off (sequence generation)

inner-loop counters only (per-entry overhead above adds < 0.07/elem):

| build mode | SBOX/elem | MULT/elem | XOR/elem | INV/elem (ours) | analytic |
|---|---:|---:|---:|---:|---|
| naive   | **6.000** | 10.000 | 11.000 | 1 | 4 S-boxes (x_3 cells) + 2 (x_4 cells) |
| hoisted | **3.694** | 6.925 | 7.925 | 1 | 3 + 3/(n1 n2), a-side shared over (b1,b2) |
| gray    | **2.218** | 2.437 | 6.273 | 1 | 2 + 4/Nblk, Nblk = E[n0 n1 n2 n3] = 18.3 |
| gray2   | **1.305** | 1.359 | 4.568 | 1 | 1.25 + 4/N, N = 4 x 18.3 = 73 : (1/2 x 1 + 1/4 x 1 + 1/4 x 2) |

DDT of the 6-bit inversion S-box: per nonzero input difference, 33 zeros /
30 twos / 1 four -> E[n | n>0] = 64/31 = 2.065, block size E[prod n] = 2.065^4 =
18.2 (measured 18.3). Gray step cost: x_3-cell flip = 2 S-boxes/element
(round-3 cell then the round-4 cell it feeds), x_4-layer partner flip = 1
S-box/element - exactly the AES 2/1 structure; cold start (6 S-box/elem)
once per block. The ideal AES figure (15x1+1x2)/16 = 1.0625 has the same
form as our gray2 analytic 1/2 x 1 + 1/4 x 1 + 1/4 x 2 = 1.25 (4 cheap bits of 20 in AES
vs 2 cheap digits of 6 here).

Entries-per-second gains are smaller than op-count gains because at this
scale the per-entry combine dominates: see 2c.

### 2c. Per-entry cost split (hoisted kernel, long runs): C_post at small scale

| component | ns/entry |
|---|---:|
| sequence generation (hoisted, 64 elements) | 224-239 |
| + chi canonical fingerprint, both parities (ours' C_post) | **+303** |
| + DFJ multiset hash (64 lookups+adds) | + ~1 (== cost of the checksum it replaced; negligible) |
| + Bloom insert into 2^30-bit anon filter: 14 bit-sets (ours) / 7 (DFJ) | + ~455 / + ~246 |
| full ours build (hoist) 998 ns/entry; full DFJ 486 ns/entry | |

So C_post(chi-canonical, 2 parities) ~ 1.3x the whole hoisted 64-element
sequence generation here; at AES scale an entry has 255 elements, so
per-entry fixed costs weigh ~4x heavier at small scale (see section 4).

## 3. Online side

### 3a. Is a Gray/incremental update applicable? - partially, in a weaker form

AES online (`ddt_gray_online.c`): for a fixed (pair, k_-1-guess, gamma) the
k_6 candidates form a 2^4 DDT branch cube over the **4** active anti-diagonal
key bytes; walking that cube in Gray order recomputes one cell column per
step -> iSBOX/(cand,elem) 4 -> (4+15x1)/16 = 1.19.

SR(7,2,2,6): the last round has only **2** active key cells, k_6[0] and
k_6[3]. The attack enumerates k_6[0] freely (64 values - it plays gamma's
role: it fixes dx_6[0] and hence the required output difference for cell 3)
and only **k_6[3] is DDT-constrained**, with 2 or 4 partner solutions. So
there is no multi-cell branch cube to Gray-walk: the only reusable structure
is "for fixed k_6[0], its DDT-partner k_6[3] candidates share the cell-0 peel
column". Implemented as the `gray` online option: per valid k_6[0], build the
cell-0 column A0[w]=T0[ct0[w]^k60] once (64 iSBOX-class lookups), then each
of its n(k63) in {2,4} partners recomputes only the cell-3 column.
k_-1 guesses cannot be amortized (changing k_-1 changes the delta-set itself)
- same as AES - and DFJ's u_5[0] guess is a free additive shift of all 64
elements (1 SINV/element, irreducible).

Measured partner multiplicity: 8,467,520 candidates / 4,102,400 column
builds = **2.064 partners per column** (seed 336) -> iSBOX-class per
(cand,elem) = 1 + 1/2.064 = 1.4845.

### 3b. Measured per-candidate element costs (instrumented binary, 5 seeds each)

Counted over the whole attack; "cand" = one Bloom lookup; 64 elements each.

| variant | mode | iSBOX-class (T0/T1) /elem | GLUT (L^-1 o inv) /elem | SINV /elem | XOR /elem | mset add /elem |
|---|---|---:|---:|---:|---:|---:|
| ours (chi, no u_5) | plain | **2.0000** | 0.9844 | 0.031 | 4.015 | - |
| ours (chi, no u_5) | gray  | **1.4845** | 0.9844 | 0.031 | 3.500 | - |
| DFJ (mset, u_5 guessed) | plain | 0.0313 | - | **1.016** | 2.063 | 1.000 |
| DFJ (mset, u_5 guessed) | gray  | 0.0232 | - | **1.016** | 2.055 | 1.000 |

(plus the chi psum accumulate inside chi_fp for ours: one POW-row wide-XOR per
element per parity = the POWINV analogue.) Identical numbers on all 5 seeds
to 3 decimals. Ours gray: -25.8% iSBOX-class lookups (2.00 -> 1.48), total
per-element lookups 2.98 -> 2.47 (-17%); DFJ gray: -0.008/elem - negligible,
as at AES scale the u_5 peel dominates DFJ's online cost.

### 3c. Key recovery still works (existing 2^41-bit, k=7 tables)

20 runs: 5 seeds {329,336,341,346,348} x {ours via bloom_7226_k7.bin,
DFJ via bloom_dfj_k7.bin} x {plain, gray}: **20/20 correct key
recovered**, and for every seed the candidate count, lookup count, Bloom-hit
count and key-schedule-trial count are **byte-identical between plain and
gray** (the gray path computes exactly the same sequences). Candidates per
structure: ours 7.66M-8.47M (~2^23), DFJ 490M-542M (~2^29) = 64x (the u_5
guess). Wall-clock with the instrumented binary: ours ~58 s/seed, DFJ
~295-386 s/seed (per-run logs not bundled).

## 4. Mapping measured small-scale quantities -> AES-128 ledger (paper App. sec:ddtgray, sec:total)

| AES ledger quantity | AES value (paper) | SR(7,2,2,6) measured analogue |
|---|---|---|
| c_omega^off naive: delta-difference traced through rounds 2-4, "~20 S-boxes/omega" + 1 POWINV | 20 (+1) | **6.000 SBOX/elem** (4 x_3 cells + 2 x_4 cells) + 1 INV; +10 const-mult, +11 XOR (xtime-class / free in the AES ledger) |
| c_omega^off DDT-Gray amortized: (15x1 + 1x2)/16 | 1.0625 ideal, 1.0664 measured | **gray2 1.305 SBOX/elem** measured [ideal 1/2x1 + 1/4x1 + 1/4x2 = 1.25 + 4/Nblk]; x_3-digits-only walk (2 S-box flips only): 2.218 |
| c_omega^off = S-box term + 1 accumulate | 2.06 | gray2: 1.31 + 1 INV ~ 2.3; gray: 2.22 + 1 ~ 3.2; naive 6 + 1 = 7 |
| c_omega^on cold: 4 iSBOX + 1 L^-1 + 1 POWINV | 6 | 2.000 T-lookups (iSBOX folded with x2/x3) + 0.984 GLUT (L^-1 o inverse folded) + 1 psum-row XOR/parity ~ **4.0** |
| c_omega^on DDT-Gray: (4+15x1)/16 + 1 + 1 | 3.19 | 1.4845 + 0.984 + 1 ~ **3.47**; only a 1-cell partner walk exists here (2 active key cells, only k_6[3] DDT-constrained; 2.06 partners/column) vs AES's 2^4 cube over 4 cells |
| C_post ours (chi-star canonical combine; I_{m,n} 160 before chi-star) | ~15 (160) | two-parity chi canonicalization **303 ns/entry** ~ 1.3x a 64-element hoisted sequence generation; dominates the small-scale build because an entry has only 64 elements (AES: 255) |
| C_post DFJ (multiset) | sort/hash over 256 | 64 R64 lookups + 64 adds ~ 20 ns/entry (negligible) |
| N_off, entries/wu | 2^88 | ~17.5M entries per work unit (dy1,x2_0); 4032 work units; 7.13e10 entries total |

Structural notes for the writeup: (i) the offline Gray walk ports 1-for-1
(DDT-partner digits at x_3, x_4-layer partner digits as the cheap innermost
digits, mixed radix for the DDT-4 cells - the mixed-radix option of
sec:ddtgray - and bit-identical tables); (ii) the online Gray collapses to a
single-cell partner cache at this scale because the 2x2 state exposes only
one DDT-constrained cell per candidate; (iii) C_post weighs ~4x more per
entry than at AES scale, so the measured entries/s speedups (chi: 1.00 -> 1.04 -> 1.10
M/s; DFJ: 2.06 -> 2.21 -> 2.45 M/s) understate the c_omega reduction seen in the
seq-only throughput (4.46 -> 5.47 -> 7.66 M entries/s) and in the op counts
(3.69 -> 2.22 -> 1.31 SBOX/elem).

## 5. Caveats (reported, not hidden)

* Box heavily loaded (load 110-150/192 cores) during all timing; >=3 reps, spread <=7%.
* Timing builds insert into a 2^30-bit anonymous Bloom, not the 2^41-bit file; the
  fingerprint/Bloom split in 2c measured with dedicated chionly/msonly/none emit modes.
* gray2 enumerates the same blocks in a different order (z4 partner pairs
  merged into one walk) - verified bit-identical multiset like gray (6-builder check).
* The counting binary (-DINSTR) is separate from the timing binary; counts are
  deterministic (full-wu count runs in counts/).
* No change to what the attacks compute: plain vs gray have identical candidate
  sets, lookups, hits and key-schedule trials on every seed (verified), and
  all 20 runs recovered the true key.

## 6. Reproduce

    cd AES/experiments/sr7226-e2e
    gcc -O3 -march=native -fopenmp -DBLOOM_K=7         -o sr7226_gray     sr7226_gray.c -lm
    gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -DINSTR -o sr7226_gray_cnt sr7226_gray.c -lm

    # 1. multiset identity (6 builders x 12 full work units, with fingerprint accumulators)
    ./sr7226_gray graycheck 12 1 12                  # -> gc12.log  (ALL MULTISETS IDENTICAL)

    # 2. offline op counts (full work unit 1000), and throughput grid
    ./sr7226_gray_cnt build3 {cold|hoist|gray|gray2} {chi|dfj} 1000 1001 1 64 30   # counts/
    # timing grid: 3 reps x {cold,hoist,gray,gray2} x {chi,dfj,none} + original builder,
    # then the chionly / msonly / none split, via repeated invocations of the
    # build3 mode above with the corresponding arguments

    # 3. online: plain vs gray against the existing 2^41-bit k=7 tables (read-only)
    ./sr7226_gray_cnt attack_full bloom_7226_k7.bin 41 seed=NNN [gray]
    ./sr7226_gray_cnt attack_dfj  bloom_dfj_k7.bin  41 seed=NNN [gray]
    # 5 seeds x {ours,dfj} x {plain,gray} via the two attack modes above

The verification, op-count, throughput and attack outputs are produced by the
commands above; none of that output is bundled.
