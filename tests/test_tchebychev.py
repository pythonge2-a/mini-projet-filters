import unittest
import numpy as np
from filters.snk.tchebychev import TchebychevFilter  # <-- adaptez l'import selon votre structure
# Par exemple :
# from filters.snk.tchebychev import TchebychevFilter


class TestTchebychev(unittest.TestCase):
    def setUp(self):
        self.filter_designer = TchebychevFilter()
        self.f = 2500.0     # Fréquence de coupure (Hz) pour les tests
        self.tolerance = 0.1  # Tolérance absolue pour les comparaisons (en ohms)

    # ----------------------------------------------------------------
    # Test de l'ordre 1 (passe-bas)
    # ----------------------------------------------------------------
    def test_tchebychev_order_1(self):
        order = 1
        capacitors = [10e-9]  # Un seul condensateur pour la cellule d'ordre 1

        # 1) Vérification des pôles attendus dans la table interne
        #    Pour order=1, le code renvoie [(1.9652, 0.0)] dans TCHEBYCHEV_TABLE
        poles = self.filter_designer.tchebychev_poles(order)
        expected_poles = [(1.9652, 0.0)]
        self.assertEqual(len(poles), len(expected_poles))
        for (omega0_norm, q0), (omega0_exp, q0_exp) in zip(poles, expected_poles):
            self.assertAlmostEqual(omega0_norm, omega0_exp, places=4)
            self.assertAlmostEqual(q0, q0_exp, places=4)

        # 2) Construction du filtre
        tf_lp, stages_lp = self.filter_designer.design_filter(
            order=order,
            cutoff_freq=self.f,
            filter_type="lowpass",
            c_vals=capacitors
        )
        # On s'attend à une seule cellule (ordre 1 => q=0 => 1er ordre)
        self.assertEqual(len(stages_lp), 1)

        # 3) Vérification de la résistance calculée
        #    La cellule 1er ordre retourne {'R': ..., 'C': ...}
        stage1_params = stages_lp[0]["params"]
        calculated_R = stage1_params["R"]
        # Valeur de référence
        expected_R = 3239.5
        self.assertAlmostEqual(
            calculated_R, expected_R, delta=self.tolerance,
            msg=f"R calculé={calculated_R} vs attendu={expected_R}"
        )

    # ----------------------------------------------------------------
    # Test de l'ordre 2 (passe-bas)
    # ----------------------------------------------------------------
    def test_tchebychev_order_2(self):
        order = 2
        # Ici on aura une unique cellule d'ordre 2 => besoin de 2 condos
        capacitors = [10e-9, 10e-9]

        # Vérification des pôles attendus
        # Dans la table: 2: [(1.0500, 0.9565)]
        poles = self.filter_designer.tchebychev_poles(order)
        expected_poles = [(1.0500, 0.9565)]
        self.assertEqual(len(poles), len(expected_poles))
        for (omega0_norm, q0), (omega0_exp, q0_exp) in zip(poles, expected_poles):
            self.assertAlmostEqual(omega0_norm, omega0_exp, places=4)
            self.assertAlmostEqual(q0, q0_exp, places=4)

        # Construction du filtre
        tf_lp, stages_lp = self.filter_designer.design_filter(
            order=order,
            cutoff_freq=self.f,
            filter_type="lowpass",
            c_vals=capacitors
        )
        # Il ne doit y avoir qu'une seule cellule (d'ordre 2)
        self.assertEqual(len(stages_lp), 1)

        # Vérif des résistances: cellule d'ordre 2 => params = {R1, R2, C1, C2}
        stage2_params = stages_lp[0]["params"]
        R1_calc = stage2_params["R1"]
        R2_calc = stage2_params["R2"]

        expected_resistances = {
            "R1": 3169.4,   # Valeurs de référence
            "R2": 11598.6,
        }
        self.assertAlmostEqual(R1_calc, expected_resistances["R1"], delta=self.tolerance,
                               msg=f"R1 calculé={R1_calc} vs attendu={expected_resistances['R1']}")
        self.assertAlmostEqual(R2_calc, expected_resistances["R2"], delta=self.tolerance,
                               msg=f"R2 calculé={R2_calc} vs attendu={expected_resistances['R2']}")

    # ----------------------------------------------------------------
    # Test de l'ordre 3 (passe-bas)
    # ----------------------------------------------------------------
    def test_tchebychev_order_3_lp(self):
        order = 3
        capacitors = [10e-9, 10e-9, 5e-9]

        # La table mentionne : 3: [(1.4942, 0.0), (0.9971, 2.0177)]
        poles = self.filter_designer.tchebychev_poles(order)
        expected_poles = [(1.4942, 0.0), (0.9971, 2.0177)]
        self.assertEqual(len(poles), len(expected_poles))
        for (om, q), (om_exp, q_exp) in zip(poles, expected_poles):
            self.assertAlmostEqual(om, om_exp, places=4)
            self.assertAlmostEqual(q, q_exp, places=4)

        # Construction du filtre
        tf_lp, stages_lp = self.filter_designer.design_filter(
            order=order,
            cutoff_freq=self.f,
            filter_type="lowpass",
            c_vals=capacitors
        )
        # Pour un ordre 3: on a une cellule d'ordre 1 (q=0) + une cellule d'ordre 2
        self.assertEqual(len(stages_lp), 2)

        # Stage 1 => 1er ordre => params = {"R", "C"}
        R1 = stages_lp[0]["params"]["R"]

        # Stage 2 => 2e ordre => params = {"R1", "R2", "C1", "C2"}
        R3 = stages_lp[1]["params"]["R1"]  # On l'appelle R3 pour coller aux "noms" attendus
        R4 = stages_lp[1]["params"]["R2"]  # On l'appelle R4

        expected_resistances = {
            "R1": 4260.6,
            "R3": 2109.6,
            "R4": 38647.3,
        }
        self.assertAlmostEqual(R1, expected_resistances["R1"], delta=self.tolerance,
                               msg=f"R1 calculé={R1} vs attendu={expected_resistances['R1']}")
        self.assertAlmostEqual(R3, expected_resistances["R3"], delta=self.tolerance,
                               msg=f"R3 calculé={R3} vs attendu={expected_resistances['R3']}")
        self.assertAlmostEqual(R4, expected_resistances["R4"], delta=self.tolerance,
                               msg=f"R4 calculé={R4} vs attendu={expected_resistances['R4']}")

    # ----------------------------------------------------------------
    # Test de l'ordre 4 (passe-haut)
    # ----------------------------------------------------------------
    def test_tchebychev_order_4_hp_ex_FA9(self):
        order = 4
        capacitors = [10e-9, 10e-9, 10e-9, 5e-9]

        # Dans la table: 4: [(0.5286, 0.7845), (0.9932, 3.5590)]
        # L'utilisateur semble vouloir un *passe-haut* => on utilise filter_type="highpass"
        poles = self.filter_designer.tchebychev_poles(order)
        expected_poles = [(0.5286, 0.7845), (0.9932, 3.5590)]
        self.assertEqual(len(poles), len(expected_poles))
        for (om, q), (om_exp, q_exp) in zip(poles, expected_poles):
            self.assertAlmostEqual(om, om_exp, places=4)
            self.assertAlmostEqual(q, q_exp, places=4)

        # Construction du filtre
        tf_hp, stages_hp = self.filter_designer.design_filter(
            order=order,
            cutoff_freq=self.f,
            filter_type="highpass",
            c_vals=capacitors
        )
        # Pour un ordre 4, on s'attend à 2 cellules d'ordre 2
        self.assertEqual(len(stages_hp), 2)

        # 1ʳᵉ cellule => "params" = {R1, R2, C1, C2}
        R1 = stages_hp[0]["params"]["R1"]
        R2 = stages_hp[0]["params"]["R2"]
        # 2ᵉ cellule => on va l'appeler R3 et R4
        R3 = stages_hp[1]["params"]["R1"]
        R4 = stages_hp[1]["params"]["R2"]

        expected_resistances = {
            "R1": 2144.8,
            "R2": 5280.0,
            "R3": 1184.4,
            "R4": 67509.7,
        }
        for (name, val) in zip(["R1", "R2", "R3", "R4"], [R1, R2, R3, R4]):
            self.assertAlmostEqual(val, expected_resistances[name], delta=self.tolerance,
                                   msg=f"{name} calculé={val} vs attendu={expected_resistances[name]}")

    # ----------------------------------------------------------------
    # Test de l'ordre 5 (passe-bas)
    # ----------------------------------------------------------------
    def test_tchebychev_order_5(self):
        order = 5
        capacitors = [10e-9, 10e-9, 10e-9, 10e-9, 5e-9]

        # 5: [(0.2895, 0.0), (0.6552, 1.3988), (0.9941, 5.5564)]
        poles = self.filter_designer.tchebychev_poles(order)
        expected_poles = [(0.2895, 0.0), (0.6552, 1.3988), (0.9941, 5.5564)]
        self.assertEqual(len(poles), len(expected_poles))
        for (om, q), (om_exp, q_exp) in zip(poles, expected_poles):
            self.assertAlmostEqual(om, om_exp, places=4)
            self.assertAlmostEqual(q, q_exp, places=4)

        # Construction du filtre (passe-bas)
        tf_lp, stages_lp = self.filter_designer.design_filter(
            order=order,
            cutoff_freq=self.f,
            filter_type="lowpass",
            c_vals=capacitors
        )
        # Ordre 5 => une cellule d'ordre 1 + deux cellules d'ordre 2
        self.assertEqual(len(stages_lp), 3)

        # Répartissons les résistances comme dans l'exemple attendu:
        #  - 1er bloc (q=0) => R1
        #  - 2ᵉ bloc (q!=0) => R3, R4
        #  - 3ᵉ bloc => R5, R6
        R1 = stages_lp[0]["params"]["R"]
        R3 = stages_lp[1]["params"]["R1"]
        R4 = stages_lp[1]["params"]["R2"]
        R5 = stages_lp[2]["params"]["R1"]
        R6 = stages_lp[2]["params"]["R2"]

        expected_resistances = {
            "R1": 21990.3,
            "R3": 3473.1,
            "R4": 27182.7,
            "R5": 768.4,
            "R6": 106749.2,
        }

        for (name, val) in zip(["R1", "R3", "R4", "R5", "R6"],
                               [R1, R3, R4, R5, R6]):
            self.assertAlmostEqual(val, expected_resistances[name], delta=self.tolerance,
                                   msg=f"{name} calculé={val} vs attendu={expected_resistances[name]}")

# -------------------------------------------------------------------
# Lance tous les tests
# -------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
