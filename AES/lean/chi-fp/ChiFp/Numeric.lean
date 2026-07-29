-- Copyright 2026 Anthropic PBC
-- SPDX-License-Identifier: Apache-2.0
import Mathlib

/-!
# The closed numeric inequality (Step (iv))
-/

open Finset

/-- `|256 - 2k|` as a natural number. -/
def babs (k : ℕ) : ℕ := if 2 * k ≤ 256 then 256 - 2 * k else 2 * k - 256

/-- `Bsum = ∑ₖ C(256, k) · |256 − 2k|²⁵⁵`. -/
def Bsum : ℕ := ∑ k ∈ Finset.range 257, Nat.choose 256 k * babs k ^ 255

theorem babs_cast (k : ℕ) : (babs k : ℤ) = |(256 : ℤ) - 2 * k| := by
  unfold babs
  split_ifs with h
  · rw [abs_of_nonneg (by omega)]; omega
  · rw [abs_of_neg (by omega)]; omega

theorem Bsum_desc :
    Bsum = ∑ k ∈ Finset.range 257, Nat.descFactorial 256 k / Nat.factorial k * babs k ^ 255 := by
  unfold Bsum
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [Nat.choose_eq_descFactorial_div_factorial]

theorem numeric_ineq : 256 * 255 * Bsum * 2 ^ 200 ≤ 256 ^ 255 * 2 ^ 256 := by
  rw [Bsum_desc]
  decide +kernel
