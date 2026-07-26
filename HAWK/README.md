# hawk_recover

End-to-end HAWK key recovery via the τ-cocycle lattice attack, as described in
the accompanying paper (*HAWK-n Key Recovery Reduces to SVP in Dimension
n/2+1*). Takes an encoded HAWK public key,
outputs a working secret key, verified by a sign/verify round-trip against the
NIST reference implementation.

## Usage

```bash
./scripts/setup.sh   # once: build the external dependencies with Nix
nix-shell nix/shell.nix --run \
    'sage -python -m hawk_recover path/to/pk.hex --out recovered_sk.bin'
```

(`nix-shell nix/shell.nix` provides Sage from the pinned environment; any
Sage ≥ 10.x install also works — the inner command is just
`sage -python -m hawk_recover …`.)

Input: a hex string, a `.hex` text file, or a raw binary `.bin` file.

| flag | meaning | default |
|---|---|---|
| `-o / --out` | output path for the recovered secret key | `recovered_sk.bin` |
| `-w / --workdir` | scratch + checkpoint directory | `hawk_recover_work` |
| `--ell` | cyclotomic conductor exponent | inferred from pk length |
| `--svp` | sieve backend: `auto`, `bgj` | `auto` |
| `--skip-verify` | skip the final reference-implementation check | off |

The recovered key is written in the reference implementation's *decoded*
private-key layout — `f ‖ g` (int8) `‖ F mod 2 ‖ G mod 2` (bit-packed)
`‖ hpub`, 592 bytes for HAWK-256 — the format `hawk_sign_finish_alt()`
consumes directly. It is not the 96-byte NIST-encoded sk: that encoding stores
a seed from which (f, g) are regenerated, and no attack can recover the seed;
the two are functionally equivalent for signing.

Runs are **restartable**: the cocycle, SVP (including every warm-up step and
sieve pass), and intertwiner stages checkpoint into the workdir, and
re-running resumes from the last completed one; the remaining stages take
minutes and re-run from their checkpointed inputs. A workdir belongs to one
public key — use a fresh `--workdir` for each key you attack. Unless `--skip-verify` is
given, the output file only ever contains a key that passed the reference
implementation's sign/verify (exit 0); a key that fails verification is left
at `<out>.unverified` (exit 1); if the sign/verify oracle itself cannot run,
the untested key stays at `<out>.unverified` (exit 3). A non-HAWK-256 input is
rejected up front (exit 2). With `--skip-verify`, the key is written
unverified and marked as such in the log.

The SVP schedule is adaptive: a light BKZ warm-up (β=44), then multithreaded
sieve passes; if a sieve attempt yields fewer than 2 shortest vectors, the
warm-up deepens by one step and the sieve re-runs (never observed to be
needed on real HAWK-256 keys — see Expected runtime).

## Requirements & hardware

* **Nix** (https://nixos.org/download). `scripts/setup.sh` builds every
  external dependency from pinned, hash-verified Nix derivations — the sieve
  binary, the NIST reference sources, the sign/verify shim (see `nix/` and
  the Provenance section below) — plus the test deps from PyPI with pinned
  hashes. Nothing else is fetched from anywhere.
* SageMath ≥ 10.x with fpylll (stock Sage) and numpy/scipy.
  `nix-shell nix/shell.nix` provides Sage 10.7 from the pinned nixpkgs
  (`nix/nixpkgs.nix`) — the exact reference environment of the paper; any
  Sage ≥ 10.x install also works.
* x86-64 Linux. Only the **sieve stage** has a CPU requirement — AVX512-VNNI
  (Cascade Lake or newer); all other stages run anywhere. The sieve binary
  (`bgj`, built from pinned upstream source by setup.sh — `nix/svp-tool.nix`,
  GPL-2) targets the Cascade Lake ISA baseline with the AMX kernels compiled
  in but runtime-gated, so one reproducible binary serves both AMX and
  non-AMX machines.

Reference environment for all paper results: Sage 10.7 (nixpkgs), PARI 2.17.3,
fpylll 0.6.4, numpy 2.4.4 — all from the nixpkgs revision pinned in
`nix/nixpkgs.nix`. `hawk_shim` (sign/verify oracle) is built by
`nix/hawk-shim.nix` from `ref_impl/hawk_shim.c` + the NIST reference sources;
the equivalent manual compilation (against the sources setup.sh links at
`ref_impl/hawk256_ref/`) is:

```bash
cd hawk_recover/ref_impl/hawk256_ref
gcc -O2 -I. -o ../hawk_shim ../hawk_shim.c hawk_sign.c hawk_vrfy.c hawk_kgen.c \
    ng_fxp.c ng_hawk.c ng_mp31.c ng_ntru.c ng_poly.c ng_zint31.c sha3.c -lm
```

## Package layout — paper section map

| module | paper | role |
|---|---|---|
| `params.py` | §3.1, spec Tab. 3 | parameter sets; every constant cites its source |
| `types.py` | §§2–4 | stage-boundary dataclasses; `require()` invariant checks |
| `decode.py` | §3.1, `hawk_vrfy.c` | pk bytes → public Gram matrix Q |
| `ring.py` | §3.1 | negacyclic + cyclotomic arithmetic + integer kernel (see Known Sage issue) |
| `cocycle.py` | §4 (step 1) | Q → τ-cocycle lattice Λ_τ(Q) with Gram |
| `svp.py` | §6.2 (step 2) | heuristic BKZ warm-up + sieve harvest (backend-agnostic), in place of the provable block reduction of §5 |
| `intertwiner.py` | §6.2, Remark 6.2 (step 3) | shell vectors → a-side cocycle V_τ via the intertwiner sublattice Λ_a |
| `tower.py` | §6.2, App. D (step 3) | V_τ → B_rec by τ-tower descent (replaces the Z^n-LIP solve) |
| `encode.py` | spec §3.4 | B_rec → HAWK secret-key bytes |
| `__main__.py` | §2.1 | CLI chaining all stages |
| `infra/sieve_backends/` | §6.2 | the sieve boundary + the `bgj` (svp_tool) driver |
| `infra/shim.py` | — | hpub + sign/verify against the reference implementation |
| `infra/checkpoint.py` | — | per-stage persistence |
| `evidence/` | — | inert JSON export of a run + third-party re-verification |

The tower's per-level basis extraction is the paper's Lemma D.4: the
descended module M_{V_τ} is generated by the columns of (I+V_τ)/2 and
(I−V_τ)/(2ζ), the alternating form μ = ½Tr_{L/L'}∘det is unimodular on it, and
so any (x,y) with μ(x,y)=1 — found by a Bézout step (`idealaddtoone`) — is a
free basis (see `tower.py`).

## Tests

```bash
./scripts/setup.sh           # once: build deps with Nix, install test deps
./run_tests.sh               # full suite — seconds
```

Synthetic runs at ℓ=5–6 validate the algebraic pipeline (cocycle → SVP →
intertwiner → tower) on freshly generated keys, plus decode on the shipped
HAWK-256 public keys, the sieve-output harvest bridge on a synthetic log,
and the recovered basis's key equation.
`tests/test_dispatch.py` covers the HAWK-256-scale orchestration (warm-up
resume, warm-up early exit, deepening, thread clamping) with mocked solvers;
the real sieve backend and its multi-pass restart run only in the CLI.
The full-scale recovery is the CLI run below.

`scripts/verify_lemmas.py` re-derives the paper's §4 statements — Lemmas
4.1–4.4 and Proposition 4.5 (the isometry class and the shortest-vector
count) — exactly, on a random key at ℓ = 4, 5, 6:

```bash
nix-shell nix/shell.nix --run 'sage -python scripts/verify_lemmas.py 4 5 6'
```

`test_pk_1.hex` and `test_pk_2.hex` are two HAWK-256 public keys whose
secret keys this package has recovered; attack either one to reproduce the
result, or (preferably) generate a fresh key pair with the reference implementation's
keygen and attack that. **The proof of a successful recovery is the CLI's
final step**: the recovered key signs a message that verifies under the NIST
reference implementation against the input public key.

To share a finished run or double-check one:

```bash
sage -python -m hawk_recover.evidence export <workdir> <out_dir> --pk <pk> \
    --sieve-log <workdir>/sieve_b<beta>/final.log
sage -python -m hawk_recover.evidence verify <out_dir>
```

`export` converts a run's checkpoints to inert JSON (no code executes on
loading it); `verify` re-derives the mathematical claims present in that
export — including that the output comes from the stated public key
(constraint-matrix binding). Checks whose inputs are absent (the sieve log
without `--sieve-log`, the key file without `--sk`) are reported SKIPPED.

## Expected runtime (HAWK-256, ℓ=9, n=256)

Expected stage times for a full key recovery with the shipped CLI, on the
paper's reference machine (a single 96-core Sapphire Rapids server, §6.2). Rerunning
the CLI on either included public key, or on any key you generate with the
reference implementation's keygen, follows this profile:

| stage | expected time |
|---|---|
| decode + cocycle lattice | ~6 min |
| BKZ warm-up (stops at β=44, the adaptive launch point) | ~7 min |
| sieve: sliding passes | ~64 min |
| sieve: final pass (runs to its own termination) | ~129 min |
| intertwiner (Λ_a) | ~16 min |
| tower descent (μ-symplectic) | ~4 min |
| encode + reference sign/verify | < 1 min |
| **total** | **≈ 3 h 42 min** |

The sieve's final pass dominates: the shortest vectors are typically present
in its database well before its termination criterion is met. In-progress
harvesting in the packaged backend is the planned follow-up that would close
this gap.

**Smaller conductors** (synthetic keys, single core — BKZ alone reaches λ₁
for n ≤ 128, no sieve needed). Reproduce with
`nix-shell nix/shell.nix --run 'sage -python scripts/synthetic_attack.py <ell>'`:

| ℓ | n | total |
|---|---|---|
| 5 | 16 | 0.7 s |
| 6 | 32 | 1.4 s |
| 7 | 64 | 19.2 s |
| 8 | 128 | 8.6 min |

## Known Sage issue — why cyclotomic matrices use `implementation="generic"`

Sage's `Matrix_cyclo_dense` (the default over `CyclotomicField`) computes
products with a multimodular algorithm whose height bound is wrong for
moderately large entries: at conductor 64 with ~12-bit coefficients, `A*B`
silently returns an incorrect matrix. All cyclotomic matrices here are
therefore built via `ring.gmat`/`ring.gidentity` (generic implementation) and
Galois maps applied via `ring.mat_apply` (Sage's `.apply_map()` silently
reconstructs into the buggy class). Real HAWK-256 keys have tiny coefficients
and don't trigger the bug; synthetic test instances do. A hygiene test
enforces the discipline.

## Provenance

Every external dependency is expressed as a pinned Nix derivation under
`nix/`; `scripts/setup.sh` builds them and links the results into
`hawk_recover/ref_impl/`. Nothing outside `nix/` fetches anything.

* **BGJ sieve source** (`nix/sieve-src.nix`) — Zhao–Ding–Yang, *Sieving with
  Streaming Memory Access*, fetched from the paper's GitHub repository at
  commit `62f86a53…` (GPL-2, hash-pinned). One build-time substitution raises
  the compiled-in thread cap from 112 to 192; no source is otherwise modified.
  Upstream's dependency-install scripts (which download GMP/NTL without
  integrity verification) are never run — Nix provides both as pinned build
  inputs.
* **`svp_tool`** (`nix/svp-tool.nix`) — the sieve binary, built from that
  source with clang/libc++/OpenMP and Nix-provided GMP/NTL, targeting
  the Cascade Lake ISA baseline (reproducible across build machines; pass
  `marchNative = true` for a machine-specific build).
* **HAWK reference implementation** (`nix/hawk256-ref.nix`) — the HAWK
  team's development repository (github.com/hawk-sign/dev, MIT-licensed),
  pinned to an exact commit.
* **`hawk_shim`** (`nix/hawk-shim.nix`) — the sign/verify oracle, compiled
  from `ref_impl/hawk_shim.c` (in git) + the reference implementation above.
