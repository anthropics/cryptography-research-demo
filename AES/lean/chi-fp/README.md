# chi-fp: a Lean 4 + Mathlib formalization of the `χ★` false-positive bound

This directory contains a machine-checked proof (Lean 4, Mathlib) of the
counting bound stated in Section 5.5 of the paper.

## What is proved

Let `F` be a finite field with `|F| = 256`, let `Tuple F = Fin 255 → F` be the
space of 255-tuples, and for a tuple `g` let `χ(g) : F → ZMod 2` be its parity
vector, `χ(g)_d = #{ω | g_ω = d} mod 2`.  The affine group
`AGL = {v ↦ α·v + β : α ≠ 0}` acts on parity vectors by permuting
bit-positions, and `χ★(g)` is the `AGL`-orbit of `χ(g)`.  The theorem says:
for every finite set `T` of orbits,

```
#{ g ∈ F^255 | χ★(g) ∈ T } · 2^200  ≤  |T| · 256^255,
```

i.e. a uniformly random tuple lands in a prescribed set of `|T|` orbits with
probability at most `|T| · 2^-200`.

The Lean statement (file `ChiFp/Defs.lean` for the definitions and
`ChiFp/MainTheorem.lean` for the theorem; the paper typesets both in a single
listing, the only deviation here is the split into two files so that the proof
files can sit in between) is:

```lean
import Mathlib

open Finset

variable (F : Type) [Field F] [Fintype F] [DecidableEq F]

abbrev Tuple     := Fin 255 → F                    -- 255 field elements
abbrev ParityVec := F → ZMod 2                     -- one bit per field element

def chi (g : Tuple F) : ParityVec F :=     -- χ(g)_d = (#{ω | g_ω = d}) mod 2
  fun d => ((univ.filter (fun ω => g ω = d)).card : ZMod 2)

structure AGL where                                -- v ↦ α·v + β,  α ≠ 0
  α : F
  β : F
  hα : α ≠ 0

def aglAct (f : AGL F) (v : ParityVec F) : ParityVec F :=
  fun d => v (f.α⁻¹ * (d - f.β))           -- permute bit-positions by f⁻¹

def AGLOrbit (v w : ParityVec F) : Prop := ∃ f, w = aglAct F f v

def chiStar (g : Tuple F) : Quot (AGLOrbit F) :=
  Quot.mk _ (chi F g)                      -- χ⋆(g) = AGL-orbit of χ(g)

theorem chi_fp (hF : Fintype.card F = 256)
    (T : Finset (Quot (AGLOrbit F))) :
    Nat.card {g : Tuple F // chiStar F g ∈ T} * 2 ^ 200
      ≤ T.card * 256 ^ 255
```

## Scope caveat

This theorem bounds collisions of the *full* 256-bit orbit invariant `χ★`
under a uniformly random (wrong-key) tuple: it says that a uniform tuple hits
any fixed set of `|T|` orbit classes with probability at most `|T| · 2^-200`.
It is *not* a statement about the truncated-table false-positive rate of the
attack; that rate is estimated separately in the paper.  The formal statement
is exactly the combinatorial inequality above and nothing more.

## Proof route (mirrors the paper)

* `ChiFp/CharacterSum.lean` — Step (i): the exact character-sum identity
  `2^|F| · A(v) = Σ_{S ⊆ F} (−1)^{Σ_{d∈S} v_d} (|F| − 2|S|)^255` for the fiber
  size `A(v) = #{g | χ(g) = v}`, and Step (ii): the triangle inequality and
  grouping by `k = |S|`, giving `A(v) · 2^256 ≤ Bsum` with
  `Bsum = Σ_{k=0}^{256} C(256,k) · |256 − 2k|^255`.
* `ChiFp/Orbit.lean` — Step (iii): `AGL F` is finite with
  `|AGL F| = (|F|−1)·|F| = 255·256`, `AGLOrbit` is an equivalence relation, and
  every orbit class contains at most `255·256` parity vectors.
* `ChiFp/Numeric.lean` — Step (iv): the closed numeric inequality
  `256 · 255 · Bsum · 2^200 ≤ 256^255 · 2^256`, checked by kernel evaluation
  (`decide +kernel`).  A plain-Python sanity check of the same inequality is in
  `tools/check_bound.py` (slack ≈ 0.4 bits; `log2 Bsum ≈ 2079.6`).
* `ChiFp/MainTheorem.lean` — Step (v): assembling the integer inequalities and
  cancelling the common factor `2^256`.

## Building

Toolchain `leanprover/lean4:v4.32.1`, Mathlib `v4.32.1`
(commit `520045ab14e26149ee970e2e617ca04b09bde5d6`).

```
lake exe cache get && lake build
```

## Axiom audit

`ChiFp/Audit.lean` runs `#print axioms chi_fp`; its output is

```
'chi_fp' depends on axioms: [propext, Classical.choice, Quot.sound]
```

i.e. only the three standard axioms of Lean/Mathlib.  The development contains
no `sorry` and no additional axioms.

## Files

```
ChiFp.lean                 -- top-level import
ChiFp/Defs.lean            -- definitions (as typeset in the paper)
ChiFp/Numeric.lean         -- Bsum and the closed numeric inequality
ChiFp/CharacterSum.lean    -- character-sum identity and the per-fiber bound
ChiFp/Orbit.lean           -- AGL, its cardinality, orbit equivalence, class size
ChiFp/MainTheorem.lean     -- theorem chi_fp
ChiFp/Audit.lean           -- #print axioms chi_fp
tools/check_bound.py       -- numeric sanity check of Step (iv)
```
