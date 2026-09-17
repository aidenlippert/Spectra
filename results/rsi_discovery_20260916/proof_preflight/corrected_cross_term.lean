import Std
theorem cross_term_bounds (x y : Int) :
    -(x*x + y*y) ≤ 2*(x*y) ∧ 2*(x*y) ≤ x*x + y*y := by
  have hs (z : Int) : 0 ≤ z*z := by
    by_cases hz : 0 ≤ z
    · exact Int.mul_nonneg hz hz
    · exact Int.mul_nonneg_of_nonpos_of_nonpos (by omega) (by omega)
  have hp := hs (x+y)
  have hm := hs (x-y)
  simp only [Int.add_mul, Int.mul_add] at hp
  simp only [Int.sub_mul, Int.mul_sub] at hm
  have hc : y*x = x*y := Int.mul_comm y x
  omega
#print axioms cross_term_bounds
