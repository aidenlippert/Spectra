import Std

-- These are transport and accounting lemmas, not a proof of the Python program.
theorem car_delta_transport (f : Nat → Nat)
    (injective : ∀ a b, f a = f b ↔ a = b) (a b : Nat) :
    (if f a = f b then (1 : Int) else 0) = (if a = b then 1 else 0) := by
  simp only [injective]

theorem inversion_transport (f : Nat → Nat)
    (ordered : ∀ a b, f a < f b ↔ a < b) (a b : Nat) :
    (if f b < f a then (1 : Int) else 0) = (if b < a then 1 else 0) := by
  simp only [ordered]

theorem increasing_injective (f : Nat → Nat)
    (ordered : ∀ a b, f a < f b ↔ a < b) (a b : Nat) :
    f a = f b ↔ a = b := by
  constructor
  · intro equal
    have ab := ordered a b
    have ba := ordered b a
    omega
  · intro equal
    simp only [equal]

theorem all_cross_terms (x y t u : Int) :
    (x+t)*(y+u) + (y+u)*(x+t) - (x*y+y*x) =
      x*u + t*y + t*u + y*t + u*x + u*t := by
  simp only [Int.add_mul, Int.mul_add]
  omega

theorem error_composition (actual prediction error1 error2 lower upper : Int)
    (first : lower - error1 ≤ actual)
    (second : actual ≤ upper + error1)
    (round_lo : prediction - error2 ≤ lower)
    (round_hi : upper ≤ prediction + error2) :
    prediction - (error1 + error2) ≤ actual ∧ actual ≤ prediction + (error1 + error2) := by
  omega

#print axioms car_delta_transport
#print axioms inversion_transport
#print axioms increasing_injective
#print axioms all_cross_terms
#print axioms error_composition
