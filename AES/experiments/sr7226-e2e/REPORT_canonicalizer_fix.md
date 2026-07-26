# canonfix — brute-free even-parity χ★ canonicalization for SR(7,2,2,6), GF(64)

Source: `sr7226_canon2.c` (this directory; the DFJ table-attack driver + canon2 tiers + verification/bench/combine modes).
Build: `gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o sr7226_canon2 sr7226_canon2.c -lm`.

## 0. Key finding about the starting point
The GF(64) production canonicalizer ALREADY implements the τ-branch (u²+u = A resp. A+τ,
TAU1 = 0x20, A = S3/S1³) for even |V| with S1≠0 — the 50 %-of-even-entries brute that the
GF(256) AES code had does NOT exist here. The only even-parity brute left was the
S1 = 0 tier (exact rate 1/64 of even calls): min over all 64 translates. Plus the
pre-existing α*=0 full brute (all five P_m = 0, degenerate, both parities) which is kept
(counted, unchanged, see §3). canon2 replaces the S1=0 tier by a complete algebraic rule.

## 1. GF(64) identities (all derived/verified numerically, VERIFY.log test (0))
Lucas: S_k(V+c) = Σ_{j⊆k} c^{k−j} S_j, S_0 = n mod 2, S_2j = S_j² (200,000 (V,c) × all k<32,
both parities: hold). Even n, S1 = 0: S3,S5,S9 translation-invariant;
S7(V+c) = S7 + L7(c), L7(c) = S3 c⁴ + S5 c² + S3² c (GF(2)-linear);
P5 = S1⁵, P11, P13 vanish with S1 (P-chain {5,11,13,23,31} ⟂ 63), so at S1=0
α* comes from P23 = S3 S5⁴ + S5 S9² + S3² S17 (explicit form verified) or P31.
If also S3=S5=0: T11(c) = S11 + S9 c² (bijective in c), T9,T13 invariant.
GF(64)-specific fact (found numerically, then explained): 17 ≡ 5·16 (mod 63) ⇒ S17 = S5¹⁶,
so S1=S3=S5=S9=0 ⇒ S17=0 ⇒ P23=P31=0 ⇒ the S15-coset tier of the GF(256) design is
provably UNREACHABLE at q=64 (goes to the α-brute). Stabilizer of (S1,S3,S5) for S1≠0 is
exactly {0,S1} (20,000/20,000 exhaustive) — the two-candidate rule is complete.
L7 kernel census over the 4,095 (s3,s5)≠0: dim0 1,428 / dim1 2,016 / dim2 651,
E[#candidates | S1=0 tier] = 1.969.

## 2. Chosen rule (a port of the power-sum canonicalizer design)
Frame W = α*V (α* from first nonzero P_m, unchanged); β-candidates B(W) equivariant
(B(W+c)=B(W)+c); canonical = min over the ≤4 candidate χ-permutations:
- odd |V|: {s1}.
- even, S1≠0 (E1, unchanged): {s1 u, s1(u+1)}, u²+u = A + [Tr A]·TAU1.
- even, S1=0, (s3,s5)≠0 (E2A): {c : s7 + L7(c) = min(s7 + Im L7)}, via tables
  C0TBL[s3][s5][s7] (256 KB) + KER[s3][s5] (≤4) built at init.
- S1=S3=S5=0, s9≠0 (E2B): single candidate (s11/s9)^32 (√ in GF(64)).
- else: (s7,s11,s13)≠0 → L15 coset rule ≤8 (kept, unreachable); else 64-translate brute
  (CHI_BETA_BRUTE counter) — never fired; α*=0 full brute (CHI_BRUTE_FALLBACKS) kept.
Invariance/completeness proof is the same as the powersum Lemma 1/Thm 1,2 (branch predicates
and candidate sets are AGL-invariant/equivariant; canon(V) is an element of orbit(V)).

## 3. Verification (`./sr7226_canon2 canoncheck 12000000 200000 160000000`, VERIFY.log, exit 0)
(a) invariance canon(αV+β)=canon(V): 12,000,000 even sets (9 generator families incl.
uniform, tiny/huge sizes, forced S1=0, GF(4)/GF(8)/affine-subspace cosets, multiplicative
cosets, attack-like d⁻¹-sequences) + 12,000,000 odd: 0 failures.
(b) ground truth vs exhaustive min over all 64·63 group elements: 200,000 sets —
orbit membership of the fast output + fast/brute orbit-relation agreement on paired
same-popcount sets (55,674 same-orbit pairs, 44,401 nontrivial): 0 mismatches.
(c) edge battery: |V| ∈ {0,2,4,…,64} extremes, S1=0 sets, subgroup/subspace/multiplicative
cosets: 180,640 invariance checks + orbit checks: 0 failures.
(e) targeted tiers: 640 searched S1=S5=0 sets + 320 S1=S3=S5=0 sets × 8 maps + exhaustive
orbit check (E2A 5,544 / E2B 2,466 / α-brute 630 firings): 0 failures.
(d) rates, N = 1.6×10^8 uniform even sets: S1=0 (= old 64-brute rate) 1.5620e-2
(theory 1/64 = 1.5625e-2) → new: E2A 1.5591e-2, E2B 3.3e-6, E2C 0, β-brute 0,
α-brute 2.53e-5 (pre-existing, identical in old code: all P_m=0, ~2^-12·P(P31=0)).
Attack-like (5×10^6 real d⁻¹-sequences): S1=0 1.5588e-2, α-brute 1.56e-5/call.
Residual brute after fix: α*=0 only, 1.3e-5 per call — genuinely degenerate chain
(every P_m in the Frobenius-class chain vanishes).

## 4. Cost, old vs new (rdtscp, TSC 2.2 GHz, timing/canonbench.log, 2 sessions × 3·10^6 calls,
8,192 attack-like even-parity inputs; box load 100–130/192)
| | NEW | OLD |
|---|---|---|
| even-parity canon call, mean ticks | 162–208 | 197–251 |
| S1=0 even calls (1.39 %), mean ticks | 265–363 | 2,657–3,636 (−2,390…3,275 = 10.1×) |
| chi_fp both parities per entry (incl. psums) | 495–668 ticks (225–304 ns) | 523–739 ticks (238–336 ns) |
Amortised removal ≈ 33–46 ticks per even call ≈ 8–10 ns per table entry (≈ −5 to −7 %
of the two-parity fingerprint). Offline single-thread combine, wu 1000–1002, pinned, 3 reps
each (timing/combine_final.log): seq-gen only 228 ns/entry; +χ fp old 236 ns/entry;
+χ fp new 230 ns/entry → corrected χ-combine for the graystack grid ≈ 230 ns/entry at this
load (their 303 ns was at load 110–150; relative correction −2.5…−5.5 %). Online isolation
(`NO_DISPOSE=1`, pinned core 93, seed 336, 130,889,216 lookups, new table, env
`OLD_CANON=1` switches to the old canonicalizer): new {153.5,142.0,135.7}s mean 143.7 s =
1.098 µs/lookup; old {150.9,140.3,148.1}s mean 146.4 s = 1.118 µs/lookup → −1.9 %, inside the
±6 % run spread: the online per-lookup cost is dominated by the 256 GB Bloom probes and the
64-element sequence build, not by canonicalization.

## 5. New table (`bloom_ours_canon2_k7.bin`, 256 GB, BLOOM_K=7, 2^41 bits; build_canon2_k7.log)
`build bloom_ours_canon2_k7.bin 41 0 4032 64`: 71,282,781,952 d-sequences (== old
71,282,781,952 ✓) in 2,472.5 s (41.2 min) at **64 threads** = 28.83 M entries/s
= 0.450 M/s/thread. Comparators: old k7 build 2,745.6 s at 80 threads = 25.96 M/s
(0.325 M/s/thread, load unknown); the quoted 41.3 M/s / 28.8-min run was build.log
(2^40-bit /dev/shm) at 90 threads = 0.459 M/s/thread. Bloom load 0.3597 and est FP
7.80e-04 identical to the old table; the α-brute count 1,877,282 is byte-identical to the
old build log (enumeration + fp0 path unchanged; only the S1=0 even fp1 values differ).
Build tier census: E2A 1,111,663,025 = 1.560e-2 of even calls, E2B 234,638, E2C 0,
β-brute 0, α-brute 1,877,282 (1.32e-5/call). DFJ table/sweep reused unchanged.

## 6. 50-key sweep, seeds 301..350, ≤25 parallel, new code + new table
50/50 keys recovered; χ lookups mean 2.146e7 (identical to old, canon-independent);
Bloom hits mean 37,195.7 (old 37,201.9); key-sched tries unchanged order.
Canon tiers over the sweep: E2A 1.673e7 (1.56e-2/even call), E2B 3,583, α-brute 28,313.
| sweep (wall s/key) | mean | median | min | max | vs DFJ 868.03 |
|---|---|---|---|---|---|
| DFJ k7 (existing k7 sweep) | 868.03 | 682.7 | 287.3 | 2354.0 | 1× |
| ours OLD, logged | 53.64 | 50.5 | 28.7 | 113.8 | 16.2× |
| ours OLD code+table, rerun now | 55.89 | 48.1 | 25.0 | 123.9 | 15.5× |
| ours NEW, run 1 (right after build, table hot) | 42.97 | 36.9 | 25.1 | 85.3 | **20.2×** |
| ours NEW, run 2 (paired with the OLD rerun, same hour) | 54.16 | 45.3 | 24.1 | 122.5 | 16.0× |
Honest reading: the paired A/B (old rerun 55.9 vs new run2 54.2, −3 %) and the pinned
NO_DISPOSE test agree that the canonicalization fix changes the end-to-end online wall by
only ~1–3 % at q=64 — because the production GF(64) code already had the τ-branch, the
residual S1=0 brute was 1/64 of even calls, and wall time is dominated by memory-bound
Bloom probes + sequence construction (plus 25-way contention and box load, spread 13 %).
The published ratio (16.2×) was therefore NOT inflated by an un-ledgered even-parity
brute at this scale; with the corrected (provably brute-free, min-of-≤4) canonicalization
the genuine numbers are 16.0–20.2× (median-of-runs ≈ 18×) vs DFJ k7.

## 7. Sources
sr7226_canon2.c and sweep_canon2.sh (this directory). Verification logs, the
Bloom table, sweep directories and timing outputs are produced by the run and not
bundled. New binary modes: `canoncheck Na Nb Nd`, `canonbench N`, `combine {new|old|none} wu_lo wu_hi`,
`finddeg`; attack env `OLD_CANON=1`, `NO_DISPOSE=1`.
