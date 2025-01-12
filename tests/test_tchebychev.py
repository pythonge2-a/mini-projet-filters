import unittest
from filters.snk.tchebychev import FilterBase, Filters


class TestTchebychev(unittest.TestCase):
    def setUp(self):
        self.filters = Filters()
        self.f = 2500  # Fréquence en Hz
        self.tolerance = 0.1  # Tolérance pour les comparaisons (en ohms)

    def test_tchebychev_order_1(self):
        capacitors = [10e-9]
        order = 1
        params = FilterBase.get_params("Tchebychev (r=1dB)", order, "lp")

        # Vérifie les paramètres récupérés
        expected_params = [
            (1.9652, None, 1.9652),
        ]
        self.assertEqual(len(params), len(expected_params))
        for param, expected in zip(params, expected_params):
            self.assertAlmostEqual(param[0], expected[0], places=4)
            self.assertAlmostEqual(param[1] or 0, expected[1] or 0, places=4)
            self.assertAlmostEqual(param[2], expected[2], places=4)

        # Vérifie les résistances calculées
        result = self.filters.tchebychev["lp"].order(f=self.f, capacitors=capacitors, order=order)
        expected_resistances = {
            "R1": 3239.5,
        }
        for res_name, expected_value in expected_resistances.items():
            self.assertAlmostEqual(
                result[res_name], expected_value, delta=self.tolerance,
                msg=f"Échec pour {res_name}: {result[res_name]} != {expected_value}"
            )

    def test_tchebychev_order_2(self):
        capacitors = [10e-9, 10e-9]
        order = 2
        params = FilterBase.get_params("Tchebychev (r=1dB)", order, "lp")

        # Vérifie les paramètres récupérés
        expected_params = [
            (1.0500, 0.9565, 1.0500),
        ]
        self.assertEqual(len(params), len(expected_params))
        for param, expected in zip(params, expected_params):
            self.assertAlmostEqual(param[0], expected[0], places=4)
            self.assertAlmostEqual(param[1], expected[1], places=4)
            self.assertAlmostEqual(param[2], expected[2], places=4)

        # Vérifie les résistances calculées
        result = self.filters.tchebychev["lp"].order(f=self.f, capacitors=capacitors, order=order)
        expected_resistances = {
            "R1": 3169.4,  # Valeur corrigée
            "R2": 11598.6,  # Valeur corrigée
        }
        for res_name, expected_value in expected_resistances.items():
            self.assertAlmostEqual(
                result[res_name], expected_value, delta=self.tolerance,
                msg=f"Échec pour {res_name}: {result[res_name]} != {expected_value}"
            )



    def test_tchebychev_order_3_lp(self):
        capacitors = [10e-9, 10e-9, 5e-9]
        order = 3
        params = FilterBase.get_params("Tchebychev (r=1dB)", order, "lp")

        # Vérifie les paramètres récupérés
        expected_params = [
            (0.4942, None, 0.4942),
            (0.9971, 2.0177, 0.9971),
        ]
        self.assertEqual(len(params), len(expected_params))
        for param, expected in zip(params, expected_params):
            self.assertAlmostEqual(param[0], expected[0], places=4)
            self.assertAlmostEqual(param[1] or 0, expected[1] or 0, places=4)
            self.assertAlmostEqual(param[2], expected[2], places=4)

        # Vérifie les résistances calculées
        result = self.filters.tchebychev["lp"].order(f=self.f, capacitors=capacitors, order=order)
        expected_resistances = {
            "R1": 12881.8,  # Ajoutez ici la valeur corrigée si nécessaire
            "R3": 2109.6,
            "R4": 38647.3,
        }
        for res_name, expected_value in expected_resistances.items():
            self.assertAlmostEqual(
                result[res_name], expected_value, delta=self.tolerance,
                msg=f"Échec pour {res_name}: {result[res_name]} != {expected_value}"
            )

    def test_tchebychev_order_4_hp_ex_FA9(self):
        capacitors = [10e-9, 10e-9, 10e-9, 5e-9]
        order = 4
        params = FilterBase.get_params("Tchebychev (r=1dB)", order, "hp")

        # Vérifie les paramètres récupérés
        expected_params = [
            (1 / 0.5286, 0.7845, 0.5286),
            (1 / 0.9932, 3.5590, 0.9932),
        ]
        self.assertEqual(len(params), len(expected_params))
        for param, expected in zip(params, expected_params):
            self.assertAlmostEqual(param[0], expected[0], places=4)
            self.assertAlmostEqual(param[1], expected[1], places=4)
            self.assertAlmostEqual(param[2], expected[2], places=4)

        # Vérifie les résistances calculées
        result = self.filters.tchebychev["hp"].order(f=self.f, capacitors=capacitors, order=order)
        expected_resistances = {
            "R1": 2144.8,  # Ajoutez ici la valeur corrigée si nécessaire
            "R2": 5280.0,
            "R3": 1184.4,
            "R4": 67509.7,
        }
        for res_name, expected_value in expected_resistances.items():
            self.assertAlmostEqual(
                result[res_name], expected_value, delta=self.tolerance,
                msg=f"Échec pour {res_name}: {result[res_name]} != {expected_value}"
            )

    def test_tchebychev_order_5(self):
        capacitors = [10e-9, 10e-9, 10e-9, 10e-9, 5e-9]
        order = 5
        params = FilterBase.get_params("Tchebychev (r=1dB)", order, "lp")

        # Vérifie les paramètres récupérés
        expected_params = [
            (0.2895, None, 0.2895),
            (0.6552, 1.3988, 0.6552),
            (0.9941, 5.5564, 0.9941),
        ]
        self.assertEqual(len(params), len(expected_params))
        for param, expected in zip(params, expected_params):
            self.assertAlmostEqual(param[0], expected[0], places=4)
            self.assertAlmostEqual(param[1] or 0, expected[1] or 0, places=4)
            self.assertAlmostEqual(param[2], expected[2], places=4)

        # Vérifie les résistances calculées
        result = self.filters.tchebychev["lp"].order(f=self.f, capacitors=capacitors, order=order)
        expected_resistances = {
            "R1": 21990.3,
            "R3": 3473.1,
            "R4": 27182.7,
            "R5": 768.4,
            "R6": 106749.2,
        }
        for res_name, expected_value in expected_resistances.items():
            self.assertAlmostEqual(
                result[res_name], expected_value, delta=self.tolerance,
                msg=f"Échec pour {res_name}: {result[res_name]} != {expected_value}"
            )



if __name__ == "__main__":
    unittest.main()
