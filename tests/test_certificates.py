from fractions import Fraction
import unittest
from experiments.certificates import PauliTerm, certify_hamiltonian, check_certificate

class CertificateTests(unittest.TestCase):
  def test_clique_and_rational_bounds(self):
    terms=[PauliTerm("X",3),PauliTerm("Z",4),PauliTerm("I",2)]
    c=certify_hamiltonian(terms)
    self.assertLessEqual(Fraction(c["lower"]), -3); self.assertEqual(Fraction(c["upper"]), 6)
    self.assertTrue(check_certificate(terms,c))

  def test_omitted_charge_and_tamper_rejection(self):
    terms=[PauliTerm("XX",1),PauliTerm("ZZ",1)]
    c=certify_hamiltonian(terms,omitted_coefficients=[Fraction(1,2)])
    self.assertLessEqual(Fraction(c["lower"]), -Fraction(5,2)); self.assertEqual(Fraction(c["upper"]), Fraction(3,2))
    c["upper"]="999"
    self.assertFalse(check_certificate(terms,c))

  def test_bad_partition_rejected(self):
    terms=[PauliTerm("X",1),PauliTerm("Z",1)]
    c=certify_hamiltonian(terms); c["cliques"][0]["terms"].append("X")
    self.assertFalse(check_certificate(terms,c))

  def test_nonanticommuting_terms_separate(self):
    c=certify_hamiltonian([PauliTerm("X",1),PauliTerm("X",2)])
    self.assertEqual(len(c["cliques"]), 2)

  def test_omitted_certificate_is_self_consistent(self):
    terms=[PauliTerm("X",1), PauliTerm("Z",1)]
    c=certify_hamiltonian(terms, omitted_coefficients=[Fraction(1,3)])
    self.assertTrue(check_certificate(terms, c))
    self.assertGreaterEqual(Fraction(c["upper"]), Fraction(1,3))

  def test_corrupt_sqrt_or_coefficient_is_rejected(self):
    terms=[PauliTerm("X",3), PauliTerm("Z",4)]
    c=certify_hamiltonian(terms, sqrt_digits=3)
    self.assertTrue(check_certificate(terms, c))
    c["cliques"][0]["sqrt_upper"]="4"
    self.assertFalse(check_certificate(terms, c))
    c=certify_hamiltonian(terms)
    c["cliques"][0]["terms"][0]["coeff"]="999"
    self.assertFalse(check_certificate(terms, c))

  def test_malformed_empty_certificate_is_rejected(self):
    self.assertFalse(check_certificate([], {}))

  def test_random_small_dense_reference_is_bracketed(self):
    """Dense eigvalsh is external validation, never the certificate."""
    import numpy as np
    paulis={"I":np.eye(2), "X":np.array([[0,1],[1,0.]]),
            "Y":np.array([[0,-1j],[1j,0]]), "Z":np.diag([1.,-1.])}
    rng=np.random.default_rng(7)
    for _ in range(8):
      terms=[PauliTerm("X", Fraction(int(rng.integers(-3,4)),2)),
             PauliTerm("Z", Fraction(int(rng.integers(-3,4)),2))]
      c=certify_hamiltonian(terms, omitted_coefficients=[Fraction(1,5)])
      h=sum(float(t.coeff)*paulis[t.pauli] for t in terms)
      e=float(np.linalg.eigvalsh(h)[0])
      self.assertLessEqual(float(Fraction(c["lower"])), e + 1e-10)
      self.assertGreaterEqual(float(Fraction(c["upper"])), e - 1e-10)

  def test_multi_qubit_bounds_with_actual_unknown_perturbation(self):
    import numpy as np
    from experiments.run import dense_hamiltonian
    rng = np.random.default_rng(29)
    labels = ['XXI', 'ZII', 'IZX', 'YYI', 'IZZ', 'IIX']
    for _ in range(4):
      h = {p: Fraction(int(rng.integers(-5, 6)), 3) for p in labels}
      terms = [PauliTerm(p, c) for p, c in h.items()]
      eta = Fraction(2, 7)
      cert = certify_hamiltonian(terms, omitted_coefficients=[eta], sqrt_digits=5)
      self.assertTrue(check_certificate(terms, cert))
      for sign in (-1, 1):
        perturbed = dict(h, XYZ=sign * eta)
        ground = float(np.linalg.eigvalsh(dense_hamiltonian(perturbed))[0])
        self.assertLessEqual(float(Fraction(cert['lower'])), ground + 1e-10)
        self.assertGreaterEqual(float(Fraction(cert['upper'])), ground - 1e-10)
