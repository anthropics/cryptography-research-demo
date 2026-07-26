# Artifacts for "Further Single-Key Recovery Attacks on Reduced-Round AES"

Source code for every experiment quoted in the paper.  No measurement
output or prebuilt table is bundled: each program prints or writes its own
results, so every number can be regenerated from scratch.

| directory | what it is |
|---|---|
| `experiments/` | the two experiments behind the *Measured costs* appendix: the complete black-box key recovery on SR(7,2,2,6), and the per-entry / per-candidate measurements of the real 7-round AES-128 attack |
| `c/` | standalone C programs, one per component claim (bridge, fingerprints, DDT-Gray walk, S-box caches, wrong-key randomisation, lookup floor) |
| `python/` | Python single-claim checks and the pedagogical reference implementation of the attack layers (`python/a1/`) |
| `rust/` | a second, independent implementation of the component experiments (e1–e13) |
| `figures/` | plotting scripts that consume the TSVs written by the C programs |

The paper's key-recovery claims rest on `experiments/`; `c/`, `python/` and
`rust/` check the individual components (identities, invariances, DDT
facts, cache correctness, randomisation statistics) and are not attacks.

## `experiments/`

See `experiments/README.md`.

* `sr7226-e2e/` — complete, black-box key recovery on AES SR(7,2,2,6)
  (2x2 state, 6-bit words, 24-bit key) with the full offline table, for
  three variants on one code base: ours, Derbez–Fouque–Jean with our
  table-enumeration optimizations, and Derbez–Fouque–Jean as published.
* `aes128-measured/` — the real 7-round AES-128 attack one unit at a time:
  genuine offline table entries (partial table), cycles per online
  candidate for the three variants, the AES denominators, the ordered
  fingerprint's false-positive and coverage rates, and the projection to
  the full count.

## `c/`

Single-file or self-contained per directory; build with
`gcc -O3 -march=native -fopenmp -o <bin> <file.c> -lm` (add the directory
as `-I` where a local header is used).  Usage is documented in the header
comment of each file.

| path | checks | paper location |
|---|---|---|
| `clump.c`, `chi_clump.c` | collision entropy of the `I_{m,n}` and χ★ fingerprints | fingerprint properties; app. collision entropy |
| `wrongkey.c` | wrong-key randomisation, exhaustive `2^32` sweep of `k_6` | app. wrong-key randomisation |
| `fp_solve.c` | DDT solving and branch-bit enumeration used by the table construction | §DDT-Gray |
| `bridge-consistency/bc_impl_{a,b}.c`, `bc_truehit.c` | two independent implementations of the bridge-consistency false-accept statistic | Alg. BridgeConsistency |
| `ddt-gray/ddt_gray_*.c` | DDT-aware Gray walk: incremental equals cold recomputation, amortised S-box cost, timing | §DDT-Gray, app. DDT-Gray |
| `sb-cache/sbc_repro.c` | XOR-separability behind the W / UU S-box caches; `c_iSB` closed form | §S-box caches |
| `mobius-honest/honest_online.c` | the online side with `\|K_on\| = 8`: the δ-set rebuilt from `k_{-1}` only, and the wrong-key survivor rate | §component lemmas |
| `mobius-honest/honest_rebound.c`, `degen_analysis.c`, `saes_imn.c` | pipeline checks of the offline rebound and of `I_{m,n}` matching on small-scale AES; **known-key shortcuts are documented in each header — these are pipeline checks, not key recoveries** | §component lemmas, app. degenerate cases |
| `floor/verify_ls_fp.c`, `floor/verify_incompat.c`, `floor/a1_floor.py` | the lookup floor for `\|K_on\| = 8`: `dim LS(fp) ∈ {0,8}` and `T ≥ 2^88.67` | app. lower bound |

## `python/`

Python 3.  `python/a1/` is standard library only; the top-level test
scripts, `figures/`, and `experiments/aes128-measured/phi/` need
`numpy`, `scipy`, and `matplotlib` (`pip install -r requirements.txt`,
exact versions verified).  `python/a1/README.md` describes the
layered reference implementation (`gf256.py`, `aes7r.py`, `chi_canon.py`,
`prop2.py`, the golden vectors in `test_vectors.py`).  Top-level scripts are
single-claim checks: bridge invariants (`test_invariants.py`,
`verify_Imn.py`), the bridge on real AES traces (`test_bridge_aes.py`),
`0^{-1}` edge behaviour (`test_zeroinv.py`), wrong-key randomisation
(`test_wrongkey.py`), clumping / collision entropy (`test_clumping*.py`),
the key-schedule inversion behind CompleteKey (`keysched_solve.py`), and
the fingerprint-prefix false-positive demo (`fp_concat_demo.py`).

## `rust/`

`cargo build --release`, then `./target/release/repro <exp> [args]`; each
run appends a timestamped block under `rust/results/`.  The experiment
list and what each one checks is in `rust/src/main.rs`: e1–e13 cover the
bridge identity, group and coset facts, the honest bridge on real AES-128
traces, wrong-key randomisation, collision entropy, χ★ invariance, the DDT
spectrum, AES sanity, per-entry timing, and the SR(7,4,4,4) pipeline.
