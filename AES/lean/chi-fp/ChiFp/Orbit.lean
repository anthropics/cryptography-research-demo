-- Copyright 2026 Anthropic PBC
-- SPDX-License-Identifier: Apache-2.0
import ChiFp.Defs

/-!
# The affine group `AGL F` and its orbits on parity vectors (Step (iii))
-/

open Finset

section AGLGroup

variable (F : Type) [Field F]

namespace AGL

/-- `AGL F` is in bijection with `{α : F // α ≠ 0} × F`. -/
def equivProd : AGL F ≃ {a : F // a ≠ 0} × F where
  toFun f := (⟨f.α, f.hα⟩, f.β)
  invFun p := ⟨p.1.1, p.2, p.1.2⟩
  left_inv f := by cases f; rfl
  right_inv p := by rfl

/-- The identity affine map. -/
def one : AGL F := ⟨1, 0, one_ne_zero⟩

/-- The inverse affine map. -/
def inv (f : AGL F) : AGL F := ⟨f.α⁻¹, -(f.α⁻¹ * f.β), inv_ne_zero f.hα⟩

/-- Composition of affine maps: acting with `comp g f` is acting with `f` then with `g`. -/
def comp (g f : AGL F) : AGL F := ⟨f.α * g.α, g.β + g.α * f.β, mul_ne_zero f.hα g.hα⟩

theorem aglAct_one (v : ParityVec F) : aglAct F (one F) v = v := by
  funext d
  simp [aglAct, one]

theorem aglAct_inv (f : AGL F) (v : ParityVec F) :
    aglAct F (inv F f) (aglAct F f v) = v := by
  funext d
  simp only [aglAct, inv]
  congr 1
  field_simp [f.hα]
  ring

theorem aglAct_comp (g f : AGL F) (v : ParityVec F) :
    aglAct F g (aglAct F f v) = aglAct F (comp F g f) v := by
  funext d
  simp only [aglAct, comp]
  congr 1
  field_simp [f.hα, g.hα]
  ring

end AGL

theorem AGLOrbit.refl (v : ParityVec F) : AGLOrbit F v v :=
  ⟨AGL.one F, (AGL.aglAct_one F v).symm⟩

theorem AGLOrbit.symm {v w : ParityVec F} (h : AGLOrbit F v w) : AGLOrbit F w v := by
  obtain ⟨f, rfl⟩ := h
  exact ⟨AGL.inv F f, (AGL.aglAct_inv F f v).symm⟩

theorem AGLOrbit.trans {u v w : ParityVec F} (huv : AGLOrbit F u v) (hvw : AGLOrbit F v w) :
    AGLOrbit F u w := by
  obtain ⟨f, rfl⟩ := huv
  obtain ⟨g, rfl⟩ := hvw
  exact ⟨AGL.comp F g f, AGL.aglAct_comp F g f u⟩

theorem AGLOrbit.equivalence : Equivalence (AGLOrbit F) :=
  ⟨AGLOrbit.refl F, AGLOrbit.symm F, AGLOrbit.trans F⟩

/-- Two parity vectors have the same class iff they are related by an affine map. -/
theorem quot_mk_eq_iff {v w : ParityVec F} :
    Quot.mk (AGLOrbit F) v = Quot.mk (AGLOrbit F) w ↔ AGLOrbit F v w := by
  constructor
  · intro h
    rw [← (AGLOrbit.equivalence F).eqvGen_eq]
    exact Quot.eqvGen_exact h
  · exact Quot.sound

end AGLGroup

section Counting

variable (F : Type) [Field F] [Fintype F] [DecidableEq F]

instance AGL.instFintype : Fintype (AGL F) := Fintype.ofEquiv _ (AGL.equivProd F).symm

theorem AGL.card_eq : Fintype.card (AGL F) = (Fintype.card F - 1) * Fintype.card F := by
  rw [Fintype.card_congr (AGL.equivProd F), Fintype.card_prod, Fintype.card_subtype_compl,
    Fintype.card_subtype_eq]

/-- Every class of `Quot (AGLOrbit F)` has at most `|AGL F| = (|F|-1)·|F|` representatives. -/
theorem card_class_le (v₀ : ParityVec F) (s : Finset (ParityVec F))
    (hs : ∀ v ∈ s, Quot.mk (AGLOrbit F) v = Quot.mk (AGLOrbit F) v₀) :
    s.card ≤ (Fintype.card F - 1) * Fintype.card F := by
  classical
  rw [← AGL.card_eq F, ← Finset.card_univ]
  have hsub : s ⊆ (univ : Finset (AGL F)).image (fun f => aglAct F f v₀) := by
    intro v hv
    obtain ⟨f, hf⟩ := (quot_mk_eq_iff F).mp (hs v hv).symm
    exact Finset.mem_image.mpr ⟨f, Finset.mem_univ _, hf.symm⟩
  exact le_trans (Finset.card_le_card hsub) Finset.card_image_le

end Counting
