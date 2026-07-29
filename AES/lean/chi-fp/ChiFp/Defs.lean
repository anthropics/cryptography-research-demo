-- Copyright 2026 Anthropic PBC
-- SPDX-License-Identifier: Apache-2.0
import Mathlib

/-!
# Definitions for the `χ★` false-positive bound

These are the definitions exactly as typeset in the paper (Section 5.5).
-/

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
