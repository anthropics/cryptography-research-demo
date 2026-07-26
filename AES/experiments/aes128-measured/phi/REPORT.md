# Φ (ξ-ordered difference-ratio fingerprint): false-positive / uniformity measurement

Target assumption (the paper's "Assumptions" section): *"a wrong candidate matches a
stored fingerprint at the random rate of 2^-8 per fingerprint byte"* (12-byte Φ,
i.e. 2^-96 per (candidate, entry) probe).  This turns the assumption into a
measurement on genuine 7-round AES-128 data, as far as such data can reach.

## 1. Exact definition used

Objects from the paper's §"The ξ-ordered fingerprint", in the round numbering
of the code (`AES/experiments/aes128-measured/common/common.h`: rounds
r=1..7, x_1 = P⊕k_0 — the paper's numbering shifted by one: paper k_{-1} ≡ code
k_0, paper u_6 ≡ code k_7[anti-diag], paper u_5[0]/x_5[0] ≡ code x_6[0] side):

* δ-byte β = x_2[0] (δ-set built on plaintext diagonal {0,5,10,15} via the
  k_0[diag] guess); match byte ρ = x_6[0] (peeled from ciphertext anti-diagonal
  {0,13,10,7} via the k_7[anti-diag] guess).
* Offline entry: d_η = Δx_6[0] between the reference and the member with
  post-S-box difference η = Δz_2[0] (exactly the 255-element sequence of the real
  offline generator: prop2_all / cold_E / DDT-Gray walk).  D_j := d_j^{-1}
  (0^{-1}:=0), j = 1..14.
* Online: after the two-round peel, g_v := (L^{-1}(Δy_6[0]^{(v)}))^{-1}, v = Δx_2[0]
  the δ-set index; bridge: g_v = s^2·D_{η(v)} ⊕ s with s = x_6^{(0)}[0];
  η(v) = S(ξ⊕v)⊕S(ξ), ξ := x_2^{(0)}[0] (supplied by the 9th key byte k_1[0] ≡
  paper k_0[β]).  Online 14-tuple: D'_j := g_{v_j}, v_j := S^{-1}(S(ξ)⊕j)⊕ξ, j=1..14.
* **Φ(D_1..D_14)** (12 bytes): let **m := the smallest j ≥ 2 with D_j ≠ D_1**
  (normally m = 2).  Φ := ( (D_j⊕D_1)/(D_m⊕D_1) )_{j ∈ {2..14}∖{m}, ascending}.
  **D_2 = D_1 rule (explicit):** D_m is swapped into the denominator role; the other
  twelve ratios keep ascending label order.  All-equal window → Φ:=0^12 (never
  seen: 0 in 3.5·10^10 real inputs).  The rule is in-window (only the 14 stored
  labels) and Φ stays invariant under D ↦ αD⊕β, α≠0 (equalities, hence m, are
  preserved; α cancels).
* **Intrinsic consequence (any in-window rule):** output byte 0 is 0 whenever
  D_3 = D_1 *or* D_2 = D_1, so P(Φ[0]=0) = 511/65536 ≈ 2·2^-8 exactly in the
  ideal model (14 iid uniform bytes); byte k≥1 has P(0) = 2^-8 + (255/256)·2^-8(k+1).
  This is a property of Φ, not of AES; it is predicted analytically and verified below.

What is fed to Φ:
* **wrong-candidate side (online style)**, N = 2^34: per candidate a fresh random
  key, random base plaintext, random *wrong* k_0[diag0] guess (δ-set construction),
  random *wrong* k_7[anti-diag] guess (peel), random wrong k_1[0] guess giving
  ξ = w_1^{(0)}[0] ⊕ k_1[0]_guess; the 15 δ-set members actually read by Φ (reference
  v=0 and v_1..v_14) are encrypted with real 7-round AES-128 (AES-NI, verified
  bit-exact against the scalar reference cipher); D'_j = g_{v_j}; Φ → 12 bytes.
* **table-entry side (offline style)**, N = 2^30: genuine δ-difference sequences from
  the real offline generator (base_sample → rebound → 2^14 bases × 2^16-step
  DDT-Gray walk); window = labels η=1..14; D_j = d_j^{-1}; Φ → 12 bytes; each record
  keeps its 14 raw d's, a 64-bit hash of the full 255-sequence, and provenance.
* **ideal baseline**, N = 2^34: Φ of 14 iid uniform bytes (comparator that isolates
  Φ-structure from AES-data effects).

## 2. Machine / commands

the benchmark node (2× Xeon Platinum 8481C, 208 threads, ~1.8 TB RAM; ~55 cores
busy with other work during the run); on-disk footprint: wrong/ 449 GB,
ideal/ 449 GB, table/ 57 GB buckets.  Sources (this directory): `phi.c`,
`chistats.py`, `compare2.py`, `run_all.sh`; the run writes logs/ and results/
files, which are not bundled.
```
gcc -O3 -march=native -pthread -o phi phi.c -lm
./phi gate1 n=10000000 seed=4242
./phi gate2 n=2000000 threads=96 seed=0xb21d6e
./phi table logbases=14 walkbits=16 dir=table threads=96   # 2^30 entries
./phi ideal logn=34 dir=ideal threads=96                   # 2^34
./phi wrong logn=34 dir=wrong threads=128                  # 2^34
./phi anat dir=table ; ./phi anaw dir=ideal ; ./phi anaw dir=wrong
./phi cross wdir=wrong tdir=table ; cross wdir=ideal tdir=table ; cross wdir=wrong tdir=ideal
python3 chistats.py {wrong,ideal,table}/stats.txt ; python3 compare2.py wrong ideal ; compare2.py table ideal
```
Total wall time ≈ 16 min.

## 3. Correctness gates

**(i) Invariance** — 10^7 random D-vectors (25% with forced ties to exercise the
D_2=D_1 rule) × random (α≠0, β): Φ(αD⊕β) = Φ(D) in 10,000,000/10,000,000.
**0 failures. PASS.**

**(ii) Bridge composition on real data** — 2,000,000 real right-pair instances
(random key and plaintext; true key disclosed to the harness; offline Φ from the true
parameters' d-sequence vs online Φ from the peeled g_v through the true key and the
correct ξ):
```
clean instances (no bad label, s≠0):  1,785,002 / 2,000,000 = 0.8925
   online Φ == offline Φ on clean:    1,785,002 / 1,785,002   → 0 mismatches  PASS
exception instances:                    214,998 = 0.1075
   d_j = 0 only: 101,258   d_j = s(=a) only: 101,089   both: 5,328   s=0: 7,711
   Φ matches on exception instances: 0 / 214,998 (exceptions always miss)
per-label rates: P(d_j=0)=0.00390, P(d_j=s)=0.00390  (ideal 1/256 = 0.00391)
```
**Note / discrepancy to carry back to the appendix:** the bridge has *two* bad values
per label (a_ω = 0 ⇔ d_j = s, and a_ω = s ⇔ d_j = 0), plus the s = 0 case, so the
clean (covered) fraction of the 14-label window is (254/256)^14·(255/256) = 0.8927
(measured 0.8925), i.e. the uncovered fraction is 10.7%, **not** the 1−η_cov = 5.7%
implied by the appendix's η_cov ≈ 0.943 (which corresponds to one bad value per
label).  This does not affect the uniformity question below, but it is a measured
fact the coverage accounting should reconcile.

## 4. Wrong-candidate side (N = 2^34), with ideal-model comparator

Per-byte marginal chi-square (256 cells, df=255), N = 2^34:
```
pos   chi2(wrong)   p       chi2(ideal)   p      comment
 0    6.69e7        0       6.69e7        0      P(Φ[0]=0)=0.007798 both; = 511/65536 (analytic);
                                                  vs corrected model: chi2=248, p=0.61 (wrong), 252, p=0.55 (ideal)
 1    1316.5        ~0      1217.4        ~0     intrinsic 2nd-order (P(0)=2^-8+2^-16); same in ideal
 2..11  208–288     0.07–0.99  242–307   0.01–0.70   uniform
```
Two-sample homogeneity wrong-vs-ideal (is the AES-derived Φ byte distribution the
ideal-Φ distribution?): per byte p ∈ {0.73,0.68,0.11,0.15,0.011,0.10,0.98,0.74,0.13,
0.17,0.045,0.30}; joints (0,1),(0,11),(3,9),(5,6): p = 0.99, 0.63, 0.25, 0.89.
16 tests, min p = 0.011 — consistent with H0 (no AES-specific deviation detected at
the 2^34 scale; per-cell sensitivity ≈ 10^-6 of N for marginals).

Collisions (C(N,2) = 2^67 pairs), prefix of k bytes, observed / ideal 2^67·2^-8k:
```
k  :      4            5          6        7       8      9..12 (incl. full 96)
wrong: 3.4494e10(1.0039) 1.3474e8(1.0039) 527266(1.0057) 2113(1.032) 5(8) 0 (ideal 0.03..2e-9)
ideal: 3.4494e10(1.0039) 1.3473e8(1.0038) 527084(1.0053) 2039(0.996) 4(8) 0
```
The +0.39% at k≤6 is the analytic byte-0 effect (Σ_v p_v² over byte 0 = 1.0036·2^-8),
identically present in the ideal run; k=7,8 fluctuate within Poisson error (σ≈45, 2.8).
96-bit duplicates: **0** among 2^34 fingerprints on both sides.
Entropy: collision (Rényi-2) entropy of Φ measured exactly at 32/40/48-bit
projections (≈ −0.006 bit below 8k), 55.95/64.7 bits at 56/64-bit projections
(sampling-limited), and ≥ 65.4 bits for the full Φ (0 collisions in 2^67 pairs ⇒
pairwise collision probability < 3/2^67 = 2^-65.4 at 95%).  The analytic ideal-Φ
value is H_2 ≈ 95.994 bits (collision probability ≈ 1.004·2^-96 per pair); the
data is consistent with it but cannot distinguish 2^-96 from, say, 2^-70.

## 5. Table-entry side (N = 2^30 genuine offline entries)

Per-byte chi-square: byte 0 = intrinsic 511/65536 effect (corrected-model p=0.33);
byte 1 chi2=321 (same intrinsic 2nd-order effect); bytes 2..11 p ∈ [0.01,0.98].
Two-sample table-vs-ideal: per byte p ∈ {0.42,0.48,0.86,0.08,0.63,0.99,0.014,0.05,
0.37,0.019,0.096,0.058}; joints p = 0.79, 0.16, 0.17, 0.71 — consistent.
zero-in-window rate 0.0531 (≈ 1−(255/256)^14 = 0.0533 ✓), relabel rate 2^-8.0 ✓.

Collisions (C(N,2)=2^59 pairs), prefix k bytes, observed/ideal:
k=4: 1.3473e8/1.3422e8 (1.0038); k=5: 526378 (1.0040); k=6: 2066 (1.0088);
k=7: 9 (8 ideal); k=8..12: 0.
**Full 96-bit duplicates: 0** ⇒ **structural duplicates (distinct parameter tuples,
same Φ): 0 observed in 2^59 entry-pairs** (2^14 bases × 2^16 Gray-walk states
each; within-base neighbours differ in one DDT branch and still never collide).
There is therefore nothing to classify; the instrumented classifier (identical
sequence / identical window / affine-equivalent window / other) found no group.
Pairwise Φ-collision probability on the table side < 3/2^59 = 2^-57.4 (95%).

## 6. Cross-match (wrong ≥ 2^31 probes vs table ≥ 2^30 entries)

Actually run at 2^34 probes × 2^30 entries = 2^64 probe-pairs:
```
32-bit prefix matches: observed 4,311,733,198   ideal 2^32 = 4.295e9   ratio 1.0039 (intrinsic byte-0)
64-bit prefix matches: observed 0               ideal 1
96-bit full matches:   observed 0               ideal 2^-32
```
Baselines: ideal×table (2^64 pairs): 32-bit ratio 1.0039, 64-bit 0 (ideal 1), 96-bit 0;
wrong×ideal (2^68 pairs): 32-bit ratio 1.0039, 64-bit 12 (ideal 16), 96-bit 0.
The 32-bit sanity check reproduces the 2^32 ideal to the expected 0.39% offset; the
64-bit counts (0, 0, 12 over runs with ideal 1, 1, 16) are Poisson-consistent.  From
0 full matches in 2^64 wrong×table pairs, the cross collision probability is
< 3/2^64 = 2^-62.4 (95%).

## 7. Conclusion

**SUPPORTS, down to the scale the data can reach — with one quantified, intrinsic
correction.**

* The wrong-candidate Φ's from real AES data are statistically indistinguishable
  from Φ of 14 iid uniform bytes at N = 2^34 (all marginal, pairwise-joint and
  collision statistics match the ideal-model run; two-sample p-values unremarkable),
  and the genuine table entries (N = 2^30) match it too.  So the assumption reduces
  to the properties of Φ itself on iid inputs.
* Φ on iid inputs is *not exactly* "2^-8 per byte": byte 0 carries P(0) = 511/65536
  (the relabelling rule must record the D_2=D_1 event in-window), which inflates the
  pairwise match probability by a factor 1.0039 (≈ +0.0056 bits) — immaterial for
  T_fp ≈ 2^79.7, but it is the real, analytic correction to the per-byte rate.  No
  other deviation was detected.
* Reach: the measurement bounds the pairwise 96-bit collision/match probability
  only from above — < 2^-65.4 (wrong-side self-collisions, 2^67 pairs),
  < 2^-62.4 (wrong×table cross, 2^64 pairs), < 2^-57.4 (table-side self, 2^59 pairs)
  — and pins the 56–64-bit projections to their ideal values within Poisson error.
  **It cannot certify 96-bit uniformity or the 2^-96 rate itself**: a defect living
  below ~2^-64 (e.g. a true collision rate of 2^-80) would be invisible here.  What
  *can* be said at the 2^-96 level is analytic, conditional on the (well-supported)
  iid behaviour: Φ's exact collision probability is ≈ 1.004·2^-96.
* Side finding for the appendix: bridge exception coverage on real data is 0.8925
  (two bad d-values per label plus s=0), not the 0.943 assumed; the τ-variant
  coverage should be recomputed with P(label clean) = 254/256.

## 8. Coverage (plain and tau)

Purpose: measure, on real 7-round AES-128, the right-pair coverage η of the
ξ-ordered difference-ratio fingerprint (plain, 14 labels) and of its τ-variant
(16 labels, 15 probes), so the appendix can quote measured values.
Mode `phi cov` (added to `phi.c`; same real right-pair generator as gate2:
random key + plaintext, true key disclosed to the harness, δ-byte x_2[0], match
byte s = x_6[0], offline d-sequence from the true parameters via `prop2_all`,
online g_v through the true k_0[diag]/k_7[anti-diag]/ξ; labels η = 1..16).
Runs on the benchmark node (per-run logs and the coverage aggregate are
produced by the run, not bundled); **N = 3×10^8 fresh instances**
(3 independent runs, seeds 0xc0be7a6e, 0x7a5e2917, 0x51de3c; ≈ 75 s each at
96 threads).

**Definitions.**  Label j is *bad* (bridge exception) iff d_j ∈ {0, s}
(d_j = s kills a_ω = 0 of the Möbius map, d_j = 0 has no g-image); every variant
additionally requires s ≠ 0.
* *plain*: covered iff labels 1..14 all clean and s ≠ 0.
* *τ-variant, documented window family* ("clean key + 14 τ-substituted keys"
  = 15 probes over the 16 walked labels): the adjacent-pair-deletion family
  W_t = {1..16} ∖ {t, t+1}, t = 1..15 (W_15 = {1..14} is the clean window,
  probed first, then t = 1..14).  Covered iff s ≠ 0 and some W_t is exception-free
  (equivalently: the bad set among labels 1..16 is ∅, a singleton, or an adjacent
  pair).  As a cross-check the alternative "spare-label" family
  (W_0 = {1..14}, W_j = ({1..14}∖{j}) ∪ {15}, j = 1..14; 15 probes) was also
  measured; both families have the *same* analytic coverage
  (255/256)(254/256)^14(1 + 14·2/256) because (1−p)^2 + 16p(1−p) + 15p^2 = 1 + 14p.
* *idealized best case*: any exception-free 14-subset of the first 16 labels,
  i.e. ≤ 2 bad labels among 1..16 (and s ≠ 0).

**Counts (N = 300,000,000; ± = binomial standard error; z vs analytic):**
```
eta_plain                 267,747,163/N = 0.8924905 ± 0.0000179   analytic (254/256)^14 (255/256)            = 0.8925093  z=-1.05
eta^tau (adjacent family) 297,038,947/N = 0.9901298 ± 0.0000057   analytic (255/256)(254/256)^14 (1+14·2/256)  = 0.9901275  z=+0.41
eta^tau (spare family)    297,038,687/N = 0.9901290 ± 0.0000057   (same analytic)                              = 0.9901275  z=+0.26
eta^tau (idealized)       298,753,930/N = 0.9958464 ± 0.0000037   analytic (255/256) P[Bin(16,2/256) <= 2]     = 0.9958473  z=-0.23
s=0 rate                    1,172,276/N = 0.0039076 ± 0.0000036   ideal 1/256 = 0.0039062                                   z=+0.37
per-label exception rate (labels 1..16, 4.8e9 label-instances):
   d=0: 0.0039065 ± 0.9e-6 (1/256=0.0039062); d=s (s!=0): 0.0038910 ± 0.9e-6 (255/65536=0.0038910);
   total bad: 0.0077975 ± 0.0000013 = 2/256 − 1/65536 (z=+0.21); conditional on s != 0 the rate is 2/256 exactly.
mean #bad labels per instance: 0.109177 among labels 1..14 (model 14·(2/256−2^-16) = 0.109161),
   0.124760 among labels 1..16 (model 0.124756); histogram (1..16): 0:264,683,755 1:33,280,271 2:1,962,145 3:71,875 4:1,916 5:38 6+:0.
```
The 2·10^6-instance measurement of §3 (0.8925) is re-confirmed: 0.8924905 ± 0.000018,
and the two-bad-values-per-label accounting (2/256, not 1/256) is confirmed at the
10^-6 level.  Per-run values (A/B/C) agree within 1–2 σ (plain 0.892487 / 0.892530 /
0.892454; τ 0.990131 / 0.990141 / 0.990117).

**Bridge-composition gate on covered instances** (online Φ over the chosen
window == offline Φ over the same window), all three runs:
```
plain     : 267,747,163 / 267,747,163  PASS
tau (adj) : 297,038,947 / 297,038,947  PASS   (chosen probe: clean window 89.25%, each of the 14 substituted windows ≈ 0.70%, t=14 window 5.3e-5)
tau (spare): 297,038,687 / 297,038,687 PASS
idealized : 298,753,930 / 298,753,930  PASS
```
0 mismatches in 1.13×10^9 covered-instance checks.

**Recommendation for the appendix numbers.**
* η (plain) = 0.8925 (measured 0.89249 ± 0.00002 = analytic
  (254/256)^14·(255/256) = 0.89251) — the appendix's current value.
* η^τ: quote the *documented-family* value **0.9901** (measured 0.990130 ± 0.000006;
  analytic (255/256)(254/256)^14(1+28/256) = 0.990128), which is what a 15-probe
  window family (clean + 14 substituted windows over 16 labels) achieves; it matches
  the appendix's "≈ 0.99".  The idealized any-14-of-16 best case is 0.9958 (measured
  0.995846 ± 0.000004) but would require probing all C(16,2) = 120 two-deletion
  windows (≥ 121 probes, not 15), so it should not be quoted at W_off^τ ≈ 64; it is
  the ceiling, not the τ-variant's coverage.  log2(1/0.9901) = 0.0144 bit
  (vs 0.01 quoted); the uncovered fraction of the τ-variant is 1 − η^τ = 0.0099,
  of which 0.0039 is the irreparable s = 0 event.
