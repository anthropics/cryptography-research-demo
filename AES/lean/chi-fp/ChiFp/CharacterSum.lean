-- Copyright 2026 Anthropic PBC
-- SPDX-License-Identifier: Apache-2.0
import ChiFp.Defs
import ChiFp.Numeric

/-!
# The character-sum identity and the per-fiber bound (Steps (i) and (ii))
-/

open Finset

namespace ChiFp

variable (F : Type) [Fintype F] [DecidableEq F]

/-- Parity comparison: `(m : ZMod 2) = a` iff `m + a.val` is even. -/
theorem natCast_eq_iff_even (m : ℕ) (a : ZMod 2) :
    (m : ZMod 2) = a ↔ Even (m + a.val) := by
  rw [← ZMod.natCast_eq_zero_iff_even, Nat.cast_add, ZMod.natCast_zmod_val, add_eq_zero_iff_eq_neg,
    ZMod.neg_eq_self_mod_two]

/-- `1 + (-1)^(m + a.val)` is the (scaled) indicator of `(m : ZMod 2) = a`. -/
theorem indicator_factor (m : ℕ) (a : ZMod 2) :
    (1 : ℤ) + (-1) ^ (m + a.val) = if (m : ZMod 2) = a then 2 else 0 := by
  split_ifs with h
  · rw [natCast_eq_iff_even] at h
    rw [h.neg_one_pow]; norm_num
  · rw [natCast_eq_iff_even] at h
    rw [(Nat.not_even_iff_odd.mp h).neg_one_pow]; norm_num

/-- The product over all field elements of the indicator factors detects `chi g = v`. -/
theorem prod_indicator_factor (g : Tuple F) (v : ParityVec F) :
    ∏ d : F, ((1 : ℤ) + (-1) ^ (((univ.filter (fun ω => g ω = d)).card) + (v d).val))
      = if chi F g = v then (2 : ℤ) ^ Fintype.card F else 0 := by
  simp_rw [indicator_factor]
  by_cases h : chi F g = v
  · subst h
    have hd : ∀ d : F, (((univ.filter (fun ω => g ω = d)).card : ℕ) : ZMod 2) = chi F g d :=
      fun d => rfl
    simp [hd]
  · rw [if_neg h]
    obtain ⟨d, hd⟩ :
        ∃ d : F, ¬ ((((univ.filter (fun ω => g ω = d)).card : ℕ) : ZMod 2) = v d) := by
      by_contra hall
      push Not at hall
      exact h (funext hall)
    exact Finset.prod_eq_zero (Finset.mem_univ d) (if_neg hd)

/-- Summing the indicator products over all tuples counts the fiber of `chi` over `v`. -/
theorem two_pow_mul_fiber (v : ParityVec F) :
    (2 : ℤ) ^ Fintype.card F * ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card
      = ∑ g : Tuple F,
          ∏ d : F, ((1 : ℤ) + (-1) ^ (((univ.filter (fun ω => g ω = d)).card) + (v d).val)) := by
  simp_rw [prod_indicator_factor]
  rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul, mul_comm]

/-- Expanding the product and swapping the two sums. -/
theorem expand (v : ParityVec F) :
    ∑ g : Tuple F,
        ∏ d : F, ((1 : ℤ) + (-1) ^ (((univ.filter (fun ω => g ω = d)).card) + (v d).val))
      = ∑ S ∈ (univ : Finset F).powerset, ∑ g : Tuple F,
          ∏ d ∈ S, (-1 : ℤ) ^ (((univ.filter (fun ω => g ω = d)).card) + (v d).val) := by
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun g _ => ?_
  exact Finset.prod_one_add _

/-- The sign sum over the field: `∑ₓ (±1) = |F| - 2|S|`. -/
theorem sum_sign (S : Finset F) :
    ∑ x : F, (if x ∈ S then (-1 : ℤ) else 1) = (Fintype.card F : ℤ) - 2 * S.card := by
  have h : ∀ x : F, (if x ∈ S then (-1 : ℤ) else 1) = 1 - 2 * (if x ∈ S then 1 else 0) := by
    intro x; split_ifs <;> norm_num
  simp_rw [h, Finset.sum_sub_distrib, ← Finset.mul_sum, Finset.sum_ite_mem, Finset.univ_inter,
    Finset.sum_const, Finset.card_univ]
  simp

/-- The inner sum over tuples for a fixed subset `S`. -/
theorem inner_term (v : ParityVec F) (S : Finset F) :
    ∑ g : Tuple F, ∏ d ∈ S, (-1 : ℤ) ^ (((univ.filter (fun ω => g ω = d)).card) + (v d).val)
      = (-1 : ℤ) ^ (∑ d ∈ S, (v d).val) * ((Fintype.card F : ℤ) - 2 * S.card) ^ 255 := by
  have key : ∀ g : Tuple F,
      ∏ d ∈ S, (-1 : ℤ) ^ (((univ.filter (fun ω => g ω = d)).card) + (v d).val)
        = (-1 : ℤ) ^ (∑ d ∈ S, (v d).val) *
            ∏ ω : Fin 255, (if g ω ∈ S then (-1 : ℤ) else 1) := by
    intro g
    rw [Finset.prod_pow_eq_pow_sum, Finset.sum_add_distrib, pow_add, mul_comm]
    congr 1
    have hc : ∑ d ∈ S, ((univ.filter (fun ω => g ω = d)).card)
        = ∑ ω : Fin 255, (if g ω ∈ S then 1 else 0) := by
      simp_rw [Finset.card_filter]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun ω _ => ?_
      rw [Finset.sum_ite_eq]
    rw [hc, ← Finset.prod_pow_eq_pow_sum]
    refine Finset.prod_congr rfl fun ω _ => ?_
    split_ifs <;> simp
  simp_rw [key, ← Finset.mul_sum]
  congr 1
  rw [← Fintype.prod_sum (fun (ω : Fin 255) (x : F) => if x ∈ S then (-1 : ℤ) else 1)]
  rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin, sum_sign]

/-- Step (i): the exact character-sum identity for the fiber size of `chi`. -/
theorem fiber_identity (v : ParityVec F) :
    (2 : ℤ) ^ Fintype.card F * ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card
      = ∑ S ∈ (univ : Finset F).powerset,
          (-1 : ℤ) ^ (∑ d ∈ S, (v d).val) * ((Fintype.card F : ℤ) - 2 * S.card) ^ 255 := by
  rw [two_pow_mul_fiber, expand]
  exact Finset.sum_congr rfl fun S _ => inner_term F v S

/-- Step (ii): triangle inequality and grouping the subsets by cardinality. -/
theorem fiber_bound (v : ParityVec F) :
    (2 : ℤ) ^ Fintype.card F * ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card
      ≤ ∑ k ∈ Finset.range (Fintype.card F + 1),
          ((Fintype.card F).choose k : ℤ) * |(Fintype.card F : ℤ) - 2 * k| ^ 255 := by
  rw [fiber_identity]
  calc ∑ S ∈ (univ : Finset F).powerset,
          (-1 : ℤ) ^ (∑ d ∈ S, (v d).val) * ((Fintype.card F : ℤ) - 2 * S.card) ^ 255
      ≤ ∑ S ∈ (univ : Finset F).powerset,
          |(-1 : ℤ) ^ (∑ d ∈ S, (v d).val) * ((Fintype.card F : ℤ) - 2 * S.card) ^ 255| :=
        Finset.sum_le_sum fun S _ => le_abs_self _
    _ = ∑ S ∈ (univ : Finset F).powerset, |(Fintype.card F : ℤ) - 2 * S.card| ^ 255 := by
        refine Finset.sum_congr rfl fun S _ => ?_
        rw [abs_mul, abs_pow, abs_pow, abs_neg, abs_one, one_pow, one_mul]
    _ = ∑ k ∈ Finset.range (Fintype.card F + 1),
          ((Fintype.card F).choose k : ℤ) * |(Fintype.card F : ℤ) - 2 * k| ^ 255 := by
        rw [Finset.sum_powerset, Finset.card_univ]
        refine Finset.sum_congr rfl fun j _ => ?_
        rw [Finset.sum_powersetCard j univ (fun m : ℕ => |(Fintype.card F : ℤ) - 2 * (m : ℤ)| ^ 255),
          nsmul_eq_mul, Finset.card_univ]

/-- Steps (i)+(ii) specialised to `|F| = 256`: every fiber of `chi` satisfies
`A(v) · 2²⁵⁶ ≤ Bsum`. -/
theorem fiber_card_mul_le (hF : Fintype.card F = 256) (v : ParityVec F) :
    ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card * 2 ^ 256 ≤ Bsum := by
  have h := fiber_bound F v
  rw [hF] at h
  have hB : (∑ k ∈ Finset.range (256 + 1), ((256).choose k : ℤ) * |((256 : ℕ) : ℤ) - 2 * k| ^ 255)
      = (Bsum : ℤ) := by
    unfold Bsum
    push_cast [babs_cast]
    rfl
  rw [hB] at h
  have h' : (((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card : ℤ) * 2 ^ 256
      ≤ (Bsum : ℤ) := by
    linarith
  exact_mod_cast h'

end ChiFp
