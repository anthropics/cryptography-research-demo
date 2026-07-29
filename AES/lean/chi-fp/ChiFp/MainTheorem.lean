-- Copyright 2026 Anthropic PBC
-- SPDX-License-Identifier: Apache-2.0
import ChiFp.Defs
import ChiFp.Numeric
import ChiFp.CharacterSum
import ChiFp.Orbit

/-!
# Main theorem: the false-positive bound for `χ★`
-/

open Finset

variable (F : Type) [Field F] [Fintype F] [DecidableEq F]

theorem chi_fp (hF : Fintype.card F = 256)
    (T : Finset (Quot (AGLOrbit F))) :
    Nat.card {g : Tuple F // chiStar F g ∈ T} * 2 ^ 200
      ≤ T.card * 256 ^ 255 := by
  classical
  have hN : Nat.card {g : Tuple F // chiStar F g ∈ T}
      = ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card := by
    rw [Nat.card_eq_fintype_card, Fintype.card_subtype]
  rw [hN]
  -- Step (iii-a): split the hitting tuples along the fibers of `chi`.
  have h1 : ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card
      ≤ ∑ v ∈ (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T),
          ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card := by
    rw [Finset.card_eq_sum_card_fiberwise (f := chi F)
      (s := (univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T))
      (t := (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T)) ?_]
    · refine Finset.sum_le_sum fun v _ => Finset.card_le_card ?_
      intro g hg
      rw [Finset.mem_filter] at hg ⊢
      exact ⟨Finset.mem_univ _, hg.2⟩
    · intro g hg
      simp only [Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_univ,
        true_and] at hg ⊢
      exact hg
  -- Step (i)+(ii): every fiber is small.
  have h2 : ∀ v ∈ (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T),
      ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card * 2 ^ 256 ≤ Bsum :=
    fun v _ => ChiFp.fiber_card_mul_le F hF v
  -- Step (iii-b): each class contains at most `255 * 256` parity vectors.
  have h3 : ((univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T)).card
      ≤ T.card * (255 * 256) := by
    rw [Finset.card_eq_sum_card_fiberwise
      (f := fun v : ParityVec F => Quot.mk (AGLOrbit F) v)
      (s := (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T))
      (t := T) ?_]
    · calc ∑ c ∈ T, (((univ : Finset (ParityVec F)).filter
            (fun v => Quot.mk (AGLOrbit F) v ∈ T)).filter
            (fun v => Quot.mk (AGLOrbit F) v = c)).card
          ≤ ∑ c ∈ T, (255 * 256) := by
            refine Finset.sum_le_sum fun c _ => ?_
            rcases (((univ : Finset (ParityVec F)).filter
                (fun v => Quot.mk (AGLOrbit F) v ∈ T)).filter
                (fun v => Quot.mk (AGLOrbit F) v = c)).eq_empty_or_nonempty
              with hne | ⟨v₀, hv₀⟩
            · rw [hne, Finset.card_empty]
              exact Nat.zero_le _
            · rw [Finset.mem_filter] at hv₀
              refine le_trans (card_class_le F v₀ _ ?_) (le_of_eq ?_)
              · intro v hv
                rw [Finset.mem_filter] at hv
                rw [hv.2, hv₀.2]
              · rw [hF]
        _ = T.card * (255 * 256) := by rw [Finset.sum_const, smul_eq_mul]
    · intro v hv
      simp only [Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_univ,
        true_and] at hv ⊢
      exact hv
  -- Assemble the integer inequalities.
  have h4 : ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card * 2 ^ 256
      ≤ T.card * (255 * 256) * Bsum := by
    calc ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card * 2 ^ 256
        ≤ (∑ v ∈ (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T),
            ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card) * 2 ^ 256 :=
          Nat.mul_le_mul_right _ h1
      _ = ∑ v ∈ (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T),
            ((univ : Finset (Tuple F)).filter (fun g => chi F g = v)).card * 2 ^ 256 :=
          Finset.sum_mul _ _ _
      _ ≤ ∑ v ∈ (univ : Finset (ParityVec F)).filter (fun v => Quot.mk (AGLOrbit F) v ∈ T),
            Bsum := Finset.sum_le_sum h2
      _ = ((univ : Finset (ParityVec F)).filter
            (fun v => Quot.mk (AGLOrbit F) v ∈ T)).card * Bsum := by
          rw [Finset.sum_const, smul_eq_mul]
      _ ≤ T.card * (255 * 256) * Bsum := Nat.mul_le_mul_right _ h3
  have h5 : ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card * 2 ^ 200 * 2 ^ 256
      ≤ T.card * 256 ^ 255 * 2 ^ 256 := by
    calc ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card * 2 ^ 200 * 2 ^ 256
        = ((univ : Finset (Tuple F)).filter (fun g => chiStar F g ∈ T)).card * 2 ^ 256
            * 2 ^ 200 := by ring
      _ ≤ T.card * (255 * 256) * Bsum * 2 ^ 200 := Nat.mul_le_mul_right _ h4
      _ = T.card * (256 * 255 * Bsum * 2 ^ 200) := by ring
      _ ≤ T.card * (256 ^ 255 * 2 ^ 256) := Nat.mul_le_mul_left _ numeric_ineq
      _ = T.card * 256 ^ 255 * 2 ^ 256 := by ring
  exact Nat.le_of_mul_le_mul_right h5 (by positivity)
