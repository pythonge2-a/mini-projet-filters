import unittest
from filters.snk.bessel import lowpass,highpass

class TestBesselFilters(unittest.TestCase):
    def setUp(self):
        self.cutoff_freq = 1000  # Fréquence de coupure en Hz
        self.tolerance = 0.1  # Tolérance pour les comparaisons
        self.lowpass_instance = lowpass()  # Créer une instance de `lowpass`
        self.highpass_instance = highpass()


    def test_first_order_lowpass(self):
        order = 1
        r_vals = [159.15]  # Valeur donnée pour la résistance (valeur attendue)
        c_vals = None  # Aucun condensateur spécifié

        # Appel de la méthode avec l'instance
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        # Imprimer la structure de stages pour comprendre sa composition
        print("Structure de stages :", stages)

        self.assertAlmostEqual(
            stages[0]["params"]["R"], r_vals[0], delta=self.tolerance,
            msg=f"Erreur pour R : {stages[0]['params']['R']} != {r_vals[0]}"
        )
    def test_second_order_lowpass(self):
        order = 2
        r_vals = [350.59,116955.2]  # Résistances spécifiées
        c_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        # Vérification des paramètres pour le deuxième ordre
        self.assertAlmostEqual(stages[0]["params"]["R1"], r_vals[0], delta=self.tolerance, msg=f"Erreur pour R1 : {stages[0]['params']['R1']} != {r_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["R2"], r_vals[1], delta=self.tolerance, msg=f"Erreur pour R2 : {stages[0]['params']['R2']} != {r_vals[1]}")
    def test_third_order_lowpass(self):
        order = 3
        c_vals = [1.2e-7,0,8.35e-8,1.44e-8]  # Résistances spécifiées
        r_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour C : {stages[0]['params']['C']} != {c_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[2], delta=self.tolerance, msg=f"Erreur pour C : {stages[0]['params']['C']} != {c_vals[2]}")
        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[3], delta=self.tolerance, msg=f"Erreur pour C : {stages[0]['params']['C']} != {c_vals[3]}")
    def test_fourth_order_lowpass(self):
        order = 4
        c_vals = [6.38e-8,1.93e-8,8.79e-8,1.11e-8]  # Résistances spécifiées
        r_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour C1 : {stages[0]['params']['C1']} != {c_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[1], delta=self.tolerance, msg=f"Erreur pour C2 : {stages[0]['params']['C2']} != {c_vals[1]}")
        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[2], delta=self.tolerance, msg=f"Erreur pour C1 : {stages[0]['params']['C1']} != {c_vals[2]}")
       
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[3], delta=self.tolerance, msg=f"Erreur pour C2 : {stages[0]['params']['C2']} != {c_vals[3]}")  
    
    def test_first_order_highpass(self):
        order = 1
        c_vals = [1.59e-7]  # Résistances spécifiées
        r_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.highpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour C : {stages[0]['params']['C']} != {c_vals[0]}")
    def test_second_order_highpass(self):
        order = 2
        r_vals = [350.59,116955.2]  # Résistances spécifiées
        c_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.highpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals,c_vals=c_vals)

   
        self.assertAlmostEqual(stages[0]["params"]["R1"], r_vals[0], delta=self.tolerance, msg=f"Erreur pour R1 : {stages[0]['params']['R1']} != {r_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["R2"], r_vals[1], delta=self.tolerance, msg=f"Erreur pour R2 : {stages[0]['params']['R2']} != {r_vals[1]}")

    def test_third_order_highpass(self):
        order = 3
        c_vals = [2.10e-7,0,1.67e-8,3.16e-7]  # Résistances spécifiées
        r_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.highpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals,c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour R : {stages[0]['params']['C']} != {c_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[2], delta=self.tolerance, msg=f"Erreur pour R1 : {stages[0]['params']['C']} != {c_vals[2]}")
        self.assertAlmostEqual(stages[0]["params"]["C"], c_vals[3], delta=self.tolerance, msg=f"Erreur pour R2 : {stages[0]['params']['C']} != {c_vals[3]}")

    def test_fourth_order_highpass(self):
        order = 4
        c_vals = [1.22e-8,4.24e-7,2.21e-8,2.94e-7]  # Résistances spécifiées
        r_vals = None # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.highpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals,c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour C1 : {stages[0]['params']['C1']} != {c_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[1], delta=self.tolerance, msg=f"Erreur pour C2 : {stages[0]['params']['C2']} != {c_vals[1]}")
        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[2], delta=self.tolerance, msg=f"Erreur pour C1 : {stages[0]['params']['C1']} != {c_vals[2]}")
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[3], delta=self.tolerance, msg=f"Erreur pour C2 : {stages[0]['params']['C2']} != {c_vals[3]}")

if __name__ == "__main__":
    unittest.main()
