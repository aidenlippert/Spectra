import Std
theorem compose_lower_errors (value lower first second : Int)
    (h : lower ≤ value + first + second) :
    lower - first - second ≤ value := by omega
#print axioms compose_lower_errors
