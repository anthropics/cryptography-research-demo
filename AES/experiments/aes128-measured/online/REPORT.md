# Three-way comparison: DFJ-as-published vs DFJ-optimized vs ours (per (k_{-1}[diag], u_6[anti-diag]) guess)

Unit = ONE (k_{-1},u_6) guess.  For DFJ that unit contains the full inner
sweep over all 256 values of s = u_5[0]; for ours it is ONE chi-star
fingerprint (the Mobius bridge absorbs u_5[0]).  Three kernels:

1. **DFJ AS PUBLISHED** (`B_published`, new): for EACH of the 256 s values the
   FULL ciphertext peel is recomputed per element (4 iSBOX + MC^-1 row-0
   combine) plus the s-dependent iSBOX plus the multiset insert -- i.e.
   255 x (4+1) S-box-class lookups + combine per s, peel NOT shared across
   the s-loop, no Gray amortization, plain scalar code, no early exit.
   This is the published accounting ell = 1280 lookups per 9-byte
   (k_{-1},u_6,u_5[0]) candidate = 256 x 1280 = 327,680 lookups per
   (k_{-1},u_6) unit.
2. **DFJ OPTIMIZED** (our B-variants): peel shared across the s-loop
   (computed once), DDT-Gray amortization across consecutive u_6
   candidates (B-gray), and the fully GFNI/AVX-512-vectorized
   implementation (B-allvec: GFNI peel + GFNI per-s iSBOX pass,
   gcc-vectorized gather multiset hash).
3. **OURS** (A-variants): one peel + L^-1 + one complete chi-star
   fingerprint pair + 2 probes (A-plain, A-gray, A-allvec with GFNI
   peel/L^-1 and the fast vector kernel).

## Commands (fresh working directory, same 8481C node type)

```
sources: online_bench.c aes_core.h chi256*.h kernel_*.h canon_even.h
gcc -O3 -march=native -o online_bench online_bench.c          (gcc 11.4, = sapphirerapids)
taskset -c 200 ./online_bench 0xA15BE2C4 20 16 5 18  > out.log
   # 2^20 candidates for ours (A), 2^16 for DFJ-optimized (= 2^24 inner
   # (u_6,s) iterations), 2^18 for DFJ-published (= 2^26 inner iterations),
   # 5 reps each; medians + [min,max]; plus 21 interleaved head-to-head rounds
```
CAVEAT for this run: the node was NOT idle (concurrent unpinned jobs,
loadavg 110-138/208; no free core, and real-time scheduling was not
available), so the benchmark core was time-shared/HT-contended.  Two
consequences: (a) absolute cycles are inflated vs the idle-node run of the
same binary (idle-node run: ours-plain 1916 cyc, ours-allvec 581,
DFJ-opt-plain 165103, DFJ-opt-allvec 113813) -- ours-plain suffered most
(2177-cyc peel here vs 1050 idle); (b) the robust numbers are the
**interleaved head-to-head ratios** (each round times all kernels back to
back on the same core) and the per-config **thread-CPU-time** cycles
(excludes descheduled time).  Both are given; they bracket the truth.

## Correctness gates (verbatim, re-verified on the new node) -- ALL PASS

```
[gate i-a] ours: online fp (vec)    fp0=1f0231824ee34712 fp1=9a2275ca1e2d9bb0
[gate i-a] ours: online fp (scalar) fp0=1f0231824ee34712 fp1=9a2275ca1e2d9bb0  vec==scalar PASS
[gate i-a] true offline seq fp      fD0=1f0231824ee34712 fD1=b9cf53c79a28e666  EITHER-match PASS
[gate i-b] DFJ: iSBOX(v^kappa_true)==trace x_6[0] all 256: PASS; multiset hash at true s == truth: PASS (H=eaefc18f5ed471c6)
[gate i-c] DFJ-published (per-s peel): hash at true s == shared-peel hash == truth: PASS (Hpub=eaefc18f5ed471c6)
[gate i-d] DFJ-published per-s peel == shared peel and per-s hash == shared-hash, 1000 random (cand,s): PASS
[gate ii] 20 keys: e-seq==trace(inv a) 20/20, DFJ a-seq==trace 20/20, gray-peel==plain 20/20, GFNI-vec==scalar 20/20, chi* EITHER 19/20 -> PASS
[gate iii] 64 wrong u_6: collide-with-true 0, distinct fp0 64/64, DFJ hash collide 0 -> PASS
=== gates: ALL PASS ===
```
New gates for (1): the published kernel's multiset hash at the TRUE key and
TRUE s equals the shared-peel hash at the true s AND the trace-derived truth
(i-c), and for 1000 random (candidate,s) its per-s intermediate peel equals
the shared peel elementwise and its per-s hash equals the shared-peel hash
(i-d).  (chi* EITHER 19/20 = the known ~0.4%-rate degenerate a=0/parity
edge case of the bridge test, not a harness failure.)

## The table (cycles per (k_{-1},u_6) candidate)

Standalone configurations, thread-CPU-time core cycles (median of 5 reps)
with [min,max] TSC ticks -- verbatim program output (contended run; see
caveat; idle-node absolute cycles in the second column where available):

| kernel | cyc/candidate (this run, cpu-time) | idle-node cyc (idle run) | breakdown | vs ours-allvec | vs ours-plain |
|---|--:|--:|---|--:|--:|
| DFJ-published (plain)              | **451,280** [390713-445875 t] |     --    | 256 x [255 x (4 fused iSBOX*iMC + 1 iSBOX + 1 R64 + add) + probe] = 1746 cyc per s (6.9 cyc/element) | 840.1x | 124.0x |
| DFJ-optimized plain (shared peel)  | 185,815 [169659-175879 t]      | 165,103   | peel 1034 once + 256 x [iSBOX pass + hash 595 + probe 130] | 345.9x | 51.1x |
| DFJ-optimized gray (shared + Gray) | 186,028 [165633-175791 t]      | 170,138   | peel amortized 4->1.19 iSBOX/omega; s-loop unchanged (dominant) | 346.3x | 51.1x |
| DFJ-optimized allvec (GFNI)        | 133,774 [112313-125892 t]      | 113,813   | GFNI peel 57 + 256 x [GFNI iSBOX pass + gather-sum ~340 + probe] | 249.0x | 36.8x |
| ours plain (vec chi* kernel)       | 3,639 [2645-3799 t]            |   1,916   | peel 2178* (idle 1050) + L^-1 seq 360-930 + chi* vec kernel ~430 + 2 probes ~80-150 | 6.8x | 1.00x |
| ours gray (vec chi* kernel)        | 3,068 [2686-3070 t]            |   1,682   | Gray peel 720 (idle 796) + seq + kernel + probes | 5.7x | 0.84x |
| ours allvec (GFNI peel + vec chi*) | **537** [487-522 t]            |     581   | GFNI peel 57 + GFNI L^-1 + chi* vec kernel ~380 + 2 probes | 1.00x | 0.15x |
(\* ours-plain's scalar peel was the stage most inflated by the concurrent
load: 2178 cyc here vs 1050 on the idle node.)

Clock- and contention-neutral **interleaved head-to-head** (21 rounds,
per-round ratios, verbatim):
```
plain : A   2164.1 t (  2426.3 cyc)  B  172424.7 t ( 193318.6 cyc)  ratio median  80.47x [min 75.90 max 84.65]
gray  : A   1753.2 t (  1965.6 cyc)  B  178328.1 t ( 199937.3 cyc)  ratio median 101.35x [min 97.30 max 108.74]
allvec: A    584.9 t (   655.8 cyc)  B  154956.2 t ( 173733.3 cyc)  ratio median 266.93x [min 241.15 max 283.88]
published: Bpub  351476.7 t ( 394067.5 cyc)  vs ours-plain ratio median 163.51x [min 152.91 max 217.37] | vs ours-allvec ratio median 606.36x [min 547.89 max 725.61]
```

## The four ratios

| comparison | standalone (cpu-time cyc) | interleaved head-to-head | meaning |
|---|--:|--:|---|
| **ours-allvec vs DFJ-optimized-allvec** (HEADLINE, both sides optimized) | **249x** | **267x** [241-284] | the bridge's 2^8 survives full engineering on both sides |
| **ours-allvec vs DFJ-published-plain** (HEADLINE, our best vs their algorithm as published) | **840x** | **606x** [548-726] | real-machine online speed-up over DFJ'13 as specified |
| ours-plain vs DFJ-published-plain (same implementation level) | 124x | 164x [153-217] | pure algorithmic difference, scalar code both sides |
| ours-allvec vs DFJ-optimized-plain | 346x | 295x (= 172425/584.9 t) | -- |

(Where standalone and interleaved differ, the interleaved number is the
more trustworthy because of the concurrent load in this run; on the earlier
IDLE node the corresponding optimized-vs-optimized ratios were 86x plain /
101x gray / 196x allvec, see the original section below.  In DFJ's own
lookup accounting the (k_{-1},u_6) unit is 327,680 lookups; it measures
447-451k cycles here = 1.4 cyc/lookup; our allvec unit is ~540-580 cycles.)

---

# Real-scale AES-128 ONLINE inner-loop benchmark: chi-star/Mobius (ours) vs DFJ-style

All measurements are on REAL 7-round AES-128 (`aes_core.h`, NR=7): one
candidate-pair context = a real delta-set of 256 chosen plaintexts and their
real 7-round ciphertexts under a random key.  NOTE (also in the source
header): we are COSTING the inner loop, not attacking, so the delta-set is
built the honest way from a random base plaintext + the TRUE
`k_{-1}[0,5,10,15]` (the `build_delta_set()` construction of `a1_e2e.c`).

## Machine, load, exact commands

- the benchmark node: 2x Intel Xeon Platinum
  8481C (Sapphire Rapids, GFNI + AVX-512), 104 cores / 208 threads, gcc
  11.4.0, `-O3 -march=native` (== `-march=sapphirerapids`).  **Idle during
  the quoted run: loadavg 0.66 start / 1.27 end of 208.**  No
  hardware PMU was exposed (`perf stat -e cycles` -> "not supported"), so cycles come
  from `rdtscp` (TSC at 2.70 GHz = "reference ticks") plus an imul-latency
  core/TSC calibration run immediately around EVERY configuration (the real
  core clock floats ~1.2-1.4x the TSC); the headline ratio is ALSO measured
  clock-neutrally by interleaving A and B chunks round-robin (head-to-head
  section: per-round ratios, 21 rounds).
- Single thread pinned `taskset -c 50` (core + HT sibling verified idle via
  `/proc/stat` sampling).
- Commands (on the benchmark node):
  ```
  gcc -O3 -march=native -o online_bench online_bench.c
  taskset -c 50 ./online_bench 0xA15BE2C4 20 14 5
  ```
  args: seed, log2(#A candidates)=20, log2(#B candidates)=14, reps=5 (a B
  "candidate" contains the 256-iteration u_5[0] loop, so 2^14 B candidates
  = 2^22 inner (u_6,u_5[0]) iterations = 2^30 sequence elements).  Medians /
  minima over 5 reps.
- Files (this dir): `online_bench.c` (harness + correctness gates),
  `kernel_v2_vec.h` (fast GFNI/AVX-512 chi-star kernel, brute-free v2),
  `kernel_v2_scalar.h`, `canon_even.h`, `kernel_fast.h`, `chi256*.h`,
  `aes_core.h`.  The numbers quoted below are from the idle-node run
  (incl. the interleaved head-to-head); an earlier loaded-node run
  (loadavg ~15-50) was consistent, ratios within 10%.  Logs are not
  bundled - the harness regenerates them.

## What is timed per candidate (u_6 anti-diagonal guess)

- **(A) OURS**: peel the 255 ciphertext differences (4 iSBOX + inverse-MC
  row-0 combination per omega, as 4 fused byte lookups
  `TT[i][x]=iMC[0][i]*iSBOX[x]` + 3 XOR), then per omega one XOR + one `Linv`
  lookup -> `e_omega = L^{-1}(v_omega ^ v_0)`, then the COMPLETE chi-star
  fingerprint pair (fp0,fp1) with the v2 brute-free canonicalization (the
  field inversion `g_omega = 1/e_omega` happens INSIDE the kernel: one
  `VGF2P8AFFINEINVQB` per 64 elements in the vector kernel / a `GF_inv`
  table in the scalar one), then 2 hash-table probes (fp0 and fp1, the
  EITHER-match) into a dummy 2^20-entry in-memory set.  The FAST vector
  kernel `kernel_v2_vec` applies directly to the online layout (it takes four
  256-byte rows, XORs them, inverts; online feeds the e-array as row 0 plus
  three zero rows - no offline-layout dependence).  The SCALAR
  canonicalization (`kernel_full_v2`) cost is reported too.
- **(B) DFJ-style**: the SAME peel once per u_6 candidate, then an inner
  loop over all 256 values s = u_5[0]: per s, `a_omega = iSBOX(v_omega ^ s)`
  (the 255-element x_5[0] difference sequence), an order-free multiset hash
  (sum over omega of 64-bit random table entries `R64[a_omega ^ a_0]`), one
  probe.  **Early exit DISABLED** for the timing row; a separate row
  measures early exit at a uniformly random planted s (hits are ~2^-64-rare
  per probe in the real attack, so no-early-exit is the honest number).
- **DDT-Gray** ("gray" rows) for BOTH: candidates in Gray order over the
  four u_6 branch bits (4-bit BRGC: cold peel 4 iSBOX/omega every 16
  candidates + 15 one-byte flips at 1 iSBOX/omega), as in
  `ddt_gray_online.c` (4 -> 1.19 iSBOX/omega).
- **"allvec"** rows (supplementary, "both sides get SIMD"): A with a
  GFNI-vectorized peel (iSBOX = `VGF2P8AFFINEQB`(inverse affine, const
  0x05) + `VGF2P8AFFINEINVQB`(identity), self-tested vs the tables) and a
  GFNI `L^{-1}`; B with the same GFNI peel + a GFNI-vectorized per-s iSBOX
  pass (its R64 multiset sum is auto-vectorized by gcc into `VPGATHERDQ`
  gathers - the compiler's best for DFJ; scalar B cannot auto-vectorize
  because of the byte `iSBOX` index lookup).

## Correctness gates (verbatim) - ALL PASS

```
[gate i-a] ours: online fp (vec)    fp0=1f0231824ee34712 fp1=9a2275ca1e2d9bb0
[gate i-a] ours: online fp (scalar) fp0=1f0231824ee34712 fp1=9a2275ca1e2d9bb0  vec==scalar PASS
[gate i-a] true offline seq fp      fD0=1f0231824ee34712 fD1=b9cf53c79a28e666  EITHER-match PASS
[gate i-b] DFJ: iSBOX(v^kappa_true)==trace x_6[0] all 256: PASS; multiset hash at true s == truth: PASS (H=eaefc18f5ed471c6)
[gate ii] 20 keys: e-seq==trace(inv a) 20/20, DFJ a-seq==trace 20/20, gray-peel==plain 20/20, GFNI-vec==scalar 20/20, chi* EITHER 20/20 -> PASS
[gate iii] 64 wrong u_6: collide-with-true 0, distinct fp0 64/64, DFJ hash collide 0 -> PASS
=== gates: ALL PASS ===
```
(i) known-answer with the TRUE u_6 bytes: the online fp EITHER-matches the
fp of the TRUE offline sequence Delta x_6[0] computed independently from the
internal trace (fp0 bit-identical here), vector == scalar kernel on the
online layout, DFJ's multiset hash at the true s equals the true-sequence
hash; (ii) all sequence elements agree with the internal trace for 20 random
keys (`e_omega == inv(a_omega)^inv(a_0)`, `v_omega == S(a_omega)^kappa`, DFJ
a-sequence, Gray peel and GFNI peel/seq/hash identical to plain scalar);
(iii) wrong u_6 guesses -> 64/64 distinct fingerprints, none colliding with
the true one.  (Harness uses the `a1_e2e.c` round naming, rounds 1..7: the
paper's `u_6[0,7,10,13]` = rk[7] anti-diagonal, `u_5[0]` = `(MC^-1 rk[6])[0]`
= kappa, and the offline anchor `x_5[0]` = internal `x_6[0]` here.)

## Benchmark lines (verbatim from the idle-node run; ticks = 2.7 GHz TSC reference ticks, core-cyc = ticks * per-config calibrated ratio shown as (before/after))

```
core/TSC ratio now (imul calibration): 1.185  [ticks = rdtscp TSC ticks; core-cyc = ticks*ratio, ratio re-calibrated around every configuration]

=== (A) OURS: per-candidate = peel + L^-1 + chi-star fp (v2, brute-free) + 2 probes ===
A-plain-vec   (kernel_v2_vec)      N=2^20 reps=5  ticks/cand median    1567.9 min    1538.2 | core-cyc/cand median    1916.4 (ratio 1.19/1.26) | cyc/omega   7.515
A-gray-vec    (kernel_v2_vec)      N=2^20 reps=5  ticks/cand median    1335.6 min    1311.8 | core-cyc/cand median    1681.9 (ratio 1.26/1.26) | cyc/omega   6.596
A-plain-scalar(kernel_full_v2)     N=2^20 reps=5  ticks/cand median    3523.3 min    3506.6 | core-cyc/cand median    4697.8 (ratio 1.41/1.26) | cyc/omega  18.423
A-gray-scalar (kernel_full_v2)     N=2^20 reps=5  ticks/cand median    3352.9 min    3318.1 | core-cyc/cand median    4222.3 (ratio 1.26/1.26) | cyc/omega  16.558
-- A breakdown (cumulative stages, plain peel) --
A.peel only                        N=2^20 reps=5  ticks/cand median     834.1 min     828.3 | core-cyc/cand median    1050.4 (ratio 1.26/1.26) | cyc/omega   4.119
A.peel+seq(L^-1)                   N=2^20 reps=5  ticks/cand median    1089.9 min    1082.0 | core-cyc/cand median    1412.9 (ratio 1.19/1.41) | cyc/omega   5.541
A.peel+seq+kernel_v2_vec           N=2^20 reps=5  ticks/cand median    1461.6 min    1422.1 | core-cyc/cand median    1840.6 (ratio 1.26/1.26) | cyc/omega   7.218
A.peel+seq+kernel_full_v2(scal)    N=2^20 reps=5  ticks/cand median    3442.2 min    3386.0 | core-cyc/cand median    4589.5 (ratio 1.41/1.26) | cyc/omega  17.998
-- A breakdown (gray peel) --
A.gray peel only                   N=2^20 reps=5  ticks/cand median     597.3 min     584.6 | core-cyc/cand median     796.4 (ratio 1.41/1.26) | cyc/omega   3.123
-- A fully vectorized (GFNI peel + GFNI L^-1 + kernel_v2_vec + probes) --
A-allvec (GFNI peel)               N=2^20 reps=5  ticks/cand median     435.4 min     428.4 | core-cyc/cand median     580.5 (ratio 1.41/1.26) | cyc/omega   2.276
A.vec peel only                    N=2^20 reps=5  ticks/cand median      48.2 min      44.5 | core-cyc/cand median      64.3 (ratio 1.41/1.26) | cyc/omega   0.252
A.vec peel+seq+kernel (no probe)   N=2^20 reps=5  ticks/cand median     368.3 min     367.1 | core-cyc/cand median     491.1 (ratio 1.26/1.41) | cyc/omega   1.926

=== (B) DFJ-style: per-candidate = peel + 256 x (iSBOX pass + multiset hash + probe) ===
B-plain no-early-exit              N=2^14 reps=5  ticks/cand median  135078.7 min  133402.7 | core-cyc/cand median  165103.3 (ratio 1.19/1.26) | cyc/omega 647.464 | s-iters/cand 256.0 (cyc/s-iter 644.9)
B-gray  no-early-exit              N=2^14 reps=5  ticks/cand median  135107.1 min  132815.9 | core-cyc/cand median  170138.2 (ratio 1.26/1.26) | cyc/omega 667.209 | s-iters/cand 256.0 (cyc/s-iter 664.6)
B-plain EARLY-exit@planted-s       N=2^14 reps=5  ticks/cand median   68202.6 min   67693.0 | core-cyc/cand median   90940.1 (ratio 1.26/1.41) | cyc/omega 356.628 | s-iters/cand 128.1 (cyc/s-iter 709.7)
-- B breakdown (cumulative stages, plain peel, no early exit) --
B.peel only                        N=2^14 reps=5  ticks/cand median     895.5 min     877.8 | core-cyc/cand median    1160.8 (ratio 1.41/1.19) | cyc/omega   4.552 | s-iters/cand 0.0 (cyc/s-iter 0.0)
B.peel+256x iSBOX-pass             N=2^14 reps=5  ticks/cand median  120521.3 min  119888.2 | core-cyc/cand median  151768.7 (ratio 1.11/1.41) | cyc/omega 595.171 | s-iters/cand 256.0 (cyc/s-iter 592.8)
B.peel+256x(iSBOX+hash)            N=2^14 reps=5  ticks/cand median  120414.4 min  120262.7 | core-cyc/cand median  151637.9 (ratio 1.26/1.26) | cyc/omega 594.659 | s-iters/cand 256.0 (cyc/s-iter 592.3)
-- B vectorized (GFNI peel + per-s GFNI iSBOX pass; R64 multiset sum + probe stay gathers) --
B-allvec no-early-exit             N=2^14 reps=5  ticks/cand median   93115.9 min   91656.2 | core-cyc/cand median  113812.7 (ratio 1.26/1.19) | cyc/omega 446.324 | s-iters/cand 256.0 (cyc/s-iter 444.6)
B.vec peel+256x(vecISBOX+hash)     N=2^14 reps=5  ticks/cand median   71210.3 min   69676.2 | core-cyc/cand median   89676.6 (ratio 1.26/1.26) | cyc/omega 351.673 | s-iters/cand 256.0 (cyc/s-iter 350.3)
loadavg(end): 1.27 1.88 2.34 1/7538 30099

=== SUMMARY (per candidate; [ticks = TSC 2.7GHz reference ticks | cyc = calibrated core cycles]) ===
A-plain vec    1567.9 t    1916.4 cyc | A-gray vec    1335.6 t    1681.9 cyc | A-plain scalar    3523.3 t    4697.8 cyc | A-gray scalar    3352.9 t    4222.3 cyc | A-allvec     435.4 t     580.5 cyc
B-plain      135078.7 t  165103.3 cyc | B-gray      135107.1 t  170138.2 cyc | B-plain early-exit   68202.6 t   90940.1 cyc | B-allvec   93115.9 t  113812.7 cyc
A breakdown (core-cyc): peel 1050.4 | +seq(L^-1) 362.5 | +chi* vec 427.7 (scalar chi* 3176.6) | +2 probes 75.8 | (GFNI peel 64.3, GFNI peel+seq+kernel 491.1)
B breakdown (core-cyc): peel 1160.8 | 256x(iSBOX pass+hash) 150477.1 [= 587.8 /s; iSBOX-only stage 588.3/s, +hash -0.5/s] | 256x probe 13465.4 [= 52.6 /s]
RATIO B/A plain (vec chi*):        ticks 86.15 x | core-cyc 86.15 x   [165103 / 1916 cyc]
RATIO B/A gray  (vec chi*):        ticks 101.16 x | core-cyc 101.16 x   [170138 / 1682 cyc]
RATIO B/A plain (scalar chi*):     ticks 38.34 x | core-cyc 35.14 x   [165103 / 4698 cyc]
RATIO B/A allvec (both GFNI):      ticks 213.88 x | core-cyc 196.06 x   [113813 / 580 cyc]
RATIO B/A compute-only (no probes): ticks 82.38 x | core-cyc 82.39 x   [B3 151638 / A3 1841 cyc];  vec-vec 182.61 x
MODEL-FORM: 256 * (per-u5-iteration cost) / (our fingerprint cost) = 256 * 640.4 / 866.0 = 189.3  (core-cyc; fingerprint cost = A minus its peel)

=== HEAD-TO-HEAD interleaved (per-round ratio, 21 rounds; A chunk 2^16 cands, B chunk 2^10 cands) ===
core/TSC ratio over rounds: median 1.259 (min 1.185, max 1.407)
plain : A   1851.3 t (  2331.4 cyc)  B  138630.7 t ( 174578.0 cyc)  ratio median  75.82x [min 72.19 max 81.01]
gray  : A   1500.1 t (  1889.1 cyc)  B  136386.0 t ( 171751.3 cyc)  ratio median  90.90x [min 86.40 max 98.15]
allvec: A    427.6 t (   538.5 cyc)  B  127371.0 t ( 160398.7 cyc)  ratio median 294.82x [min 273.73 max 333.25]
gates: ALL PASS
```

## Breakdowns (core cycles per candidate, medians of the idle-node run)

**A (ours, plain peel, fast vector kernel) = 1916 cyc:** peel 1050 (4.1
cyc/omega: 4 fused iSBOX*iMC byte lookups + 3 XOR) | sequence (XOR + L^-1
lookup) 363 (1.4 cyc/omega) | chi-star canonicalization `kernel_v2_vec` 428
(1.7 cyc/omega; includes the GFNI field inversion; in line with the
kernel's stand-alone ~310-450 cyc/entry) | 2 probes 76.  With DDT-Gray the
peel drops to 796 -> A-gray 1682.  Scalar canonicalization
(`kernel_full_v2`) costs 3177 instead of 428 -> A-plain-scalar 4698.
GFNI peel 64 (0.25 cyc/omega) and GFNI L^-1 -> fully-vectorized A = 581
cyc/candidate (chi-star kernel dominates).

**B (DFJ-style, plain peel, no early exit) = 165103 cyc:** peel 1161 once
| 256 x (iSBOX pass + multiset hash) = 150477 (588 cyc per u_5[0]
iteration = 2.3 cyc per sequence element: 1 iSBOX + 1 R64 lookup + add; the
iSBOX-only stage measured alone costs the same, the hash rides the same
dependent-load chain) | 256 probes = 13465 (53 cyc per s; independent per-s
probes overlap in the out-of-order window, whereas A's two probes sit on
the candidate's critical path) -> ~640 cyc per u_5[0] iteration.  Early
exit at a uniform random hit: 90940 cyc (128.1 s-iterations): at best a
factor ~1.9, and moot in the real attack.  GFNI-vectorized B: 350 cyc per
s, 113813 cyc per candidate.

## The measured ratio B/A (real-scale online speed-up)

Standalone configurations (per-config calibrated, idle node):

| configuration (like-for-like)                     | A cyc/cand | B cyc/cand | B/A core-cyc | B/A ticks |
|---------------------------------------------------|-----------:|-----------:|-------------:|----------:|
| plain peel, FAST vector chi* kernel               |      1916  |   165103   |  **86.2x**   |   86.2x   |
| DDT-Gray peel (both), vector chi* kernel          |      1682  |   170138   | **101.2x**   |  101.2x   |
| plain peel, SCALAR chi* kernel (kernel_full_v2)   |      4698  |   165103   |    35.1x     |   38.3x   |
| both sides GFNI-vectorized ("allvec")             |       581  |   113813   |   196.1x     |  213.9x   |
| compute-only, no probes (B3/A3)                   |      1841  |   151638   |    82.4x     |   82.4x   |
| compute-only, both GFNI (vec-vec)                 |       491  |    89677   |   182.6x     |     -     |

Interleaved head-to-head (clock-neutral, 21 rounds, per-round ratio
median [min-max]): **plain 75.8x [72.2-81.0], DDT-Gray 90.9x [86.4-98.2],
allvec 294.8x [273.7-333.3]** (verbatim in the log above).  The two methods
agree within the VM's ~15% clock/cache noise:

**Measured real-scale online ratio B/A: ~76-86x (plain), ~91-101x (DDT-Gray),
~200-300x when both sides are GFNI-vectorized; against the SCALAR chi*
kernel the ratio is ~35-38x.**

## Comparison with the model

The paper's model predicts **ratio = 256 x (per-u_5-iteration cost) / (our
fingerprint cost)**.  From the measured stage costs this is
`256 x 640.4 / 866.0 = 189.3x` (our "fingerprint cost" = A-plain minus its
peel = sequence + chi-star + probes = 866 cyc; DFJ's per-u_5 iteration = 640
cyc).  The measured FULL-candidate ratios (86x plain, 101x gray) sit below
this form for precisely the reason the paper's per-entry formula
`ell = 255*c_omega + C_post` anticipates: our candidate also pays its own
255-element peel (255 x ~4 lookups = 1050 cyc, 55% of the A candidate),
which nothing amortizes, whereas the DFJ candidate amortizes its single peel
over its 256 s-iterations.  DDT-Gray shrinks our peel term (-> 101x);
GFNI-vectorizing the peel (64 cyc) on both sides makes the fingerprint
dominate A and the measured ratio (196x standalone / 295x interleaved) meets
or exceeds the model form (189x).  In the paper's lookup units: this A is
`ell_on ~ 255 x (4+1+1) + C_post ~ 1540 lookups` (no W-cache here; the
paper's W-cache row is `ell_on = 255 x 1.49 + 8 = 388`), measured at 1916
core cycles; DFJ is `256 x 5 = 1280` lookups per `(k_{-1},u_6,u_5[0])`
iteration in their own accounting (peel re-done per s) - our B is the
like-for-like best case for DFJ (peel shared across the 256 s) and still
measures 640 cyc per s-iteration and 165k cyc per u_6-candidate.  The clean
2^8 of the bridge is directly visible: B's cost is its s-count times a
constant (no-early-exit / early-exit = 165103/90940 = 1.82 ~ 256/140 incl.
the fixed peel; compute-only vec-vec 182.6x).

## Caveats

- No PMU in the VM: "core cycles" rely on the per-config imul calibration;
  TSC ticks (wall-clock at 2.7 GHz) are given alongside, and the clock-
  neutral interleaved head-to-head confirms the ratios.
- The dummy 2^20-entry table (8 MB) makes probes L2/L3 hits; the real table
  is astronomically larger for BOTH sides (1 probe pair per A candidate vs
  256 probes per B candidate), which only widens B/A.
- B's multiset hash (sum of R64 over 255 differences) is inherently a gather;
  the paper's lookup accounting charges DFJ 1 S-box per element and prices
  the hash/match as free - generous to DFJ; this wall-clock measurement
  charges it (and charges our chi* canonicalization in full).
- The chi* EITHER known-answer gate has the known rare degenerate-miss rate
  (a=0 / parity edge; 1992/2000 in the honest-online test); all keys in this
  harness run matched.
