# Experiments behind the "Measured costs" appendix

Two experiments, source code only (no prebuilt tables or binaries).

## 1. `sr7226-e2e/` — complete key recovery on small-scale AES SR(7,2,2,6)

Black-box end-to-end key recovery on a 7-round AES with a 2x2 state of
6-bit words (24-bit key), full offline table, 50 random keys, three
attack variants on the same code base:

| file | what it is |
|---|---|
| `sr7226_variants.c` | cipher + ours (`attack_full`) + DFJ with our optimizations (`attack_dfj`); table builders `build` (ours), `builddfj` (DFJ multiset), `buildcold` (cold/naive enumeration) |
| `sr7226_canon2.c` | same, with the corrected (brute-free) even-parity chi-star canonicalizer; modes `canoncheck`, `canonbench` |
| `sr7226_dfjpub.c` | DFJ as published: cold table build + online loop with the ciphertext peel recomputed for every u_5[0] (`attack_dfj_pub`) |
| `sr7226_exact2.c` | exact-fingerprint table (no Bloom filter): `buildexact`, `attack_exact`, `attack_exact_dfj` |
| `sr7226_gray.c` | instrumented build-cost grid: cold / hoisted / DDT-Gray / Gray2 builders with S-box, mult, XOR counters; multiset-identity checker `graycheck` |
| `sort_buckets.c` | bucket sort/dedup for the exact tables |
| `*_driver.sh`, `sweep_canon2.sh`, `aggregate.sh` | build/sweep/aggregate scripts used for the runs |
| `REPORT_*.md` | measured results for each part |

Build: `gcc -O3 -march=native -fopenmp -DBLOOM_K=7 -o <bin> <file.c> -lm`.
Typical usage: `./<bin> build <bloomfile> 41 0 4032 <threads>`;
`./<bin> attack_full <bloomfile> 41 seed=NNN` (ours),
`attack_dfj ... seed=NNN` (DFJ optimized), `attack_dfj_pub ... seed=NNN`
(DFJ as published); `buildexact <dir> 0 4032 <threads>`,
`attack_exact <dir> seed=NNN`.

## 2. `aes128-measured/` — per-unit costs of the real 7-round AES-128 attack

| dir | what it measures |
|---|---|
| `online/online_bench.c` | cycles per online candidate: ours, DFJ with our optimizations, DFJ as published (plain / Gray / fully GFNI-AVX512), with correctness gates against the cipher's true internal sequence |
| `offline/` | genuine offline table entries (real parameters, DDT-Gray walk, chi-star fingerprint): `table_build.c` builds a partial table, `analyze.c`/`classify.c` do throughput, duplicate/collision and probe analysis |
| `projection/` | `offline3.c`, `online3.c`, `denom.c` measure the per-entry and per-candidate constants for all three variants plus AES denominators; `project.py` multiplies them by the counted exponents |
| `phi/phi.c` | false-positive and coverage measurements of the ordered difference-ratio fingerprint (uniformity, collisions, cross-probing, plain and tau coverage on real right pairs) |
| `common/` | shared headers: AES-128 core, chi-star kernels (scalar + GFNI/AVX-512), even-parity canonicalizer |

Build: `make` in each of `offline/`, `online/`, `phi/`, `projection/`
(shared headers in `common/`).  Each subdirectory has a `REPORT.md` with
the machine, exact commands and the measured numbers.
