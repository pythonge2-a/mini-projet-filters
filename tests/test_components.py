import unittest
from filters import Filters

class TestComponents(unittest.TestCase):
    def setUp(self):
        self.filters = Filters()

    def test_tchebychev_order_2(self):
        # Définitions des entrées
        r1, r2, c1, c2 = 1000, 2000, 1e-6, 2e-6
        print(f"\nFonction : filters.snk.tchebychev.2.components")
        print(f"Entrée : r1={r1}, r2={r2}, c1={c1}, c2={c2}")

        # Appel de la méthode
        result = getattr(self.filters.snk.tchebychev, "2").components(r1, r2, c1, c2)
        print(f"Sortie : {result}\n")

        # Résultat attendu
        expected_result = {"R1": 2000, "R2": 4000, "C1": 2e-06, "C2": 4e-06}
        self.assertEqual(result, expected_result)

if __name__ == "__main__":
    unittest.main()
