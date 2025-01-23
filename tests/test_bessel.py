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
        r_vals = [150]  # Valeur donnée pour la résistance
        c_val = 1.061032953945969e-06  # (valeur attendue)

        # Appel de la méthode avec l'instance
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals)

        self.assertAlmostEqual(
            stages[0]["params"]["C"], c_val, delta=self.tolerance,
            msg=f"Erreur pour C : {stages[0]['params']['C']} != {c_val}"
        )
    def test_second_order_lowpass(self):
        order = 2
        r_vals = [1000,10000]  # Résistances spécifiées
        c_vals = [7.93e-8,1.97e-8] # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals)

        # Vérification des paramètres pour le deuxième ordre
        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour R1 : {stages[0]['params']['C1']} != {c_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[1], delta=self.tolerance, msg=f"Erreur pour R2 : {stages[0]['params']['C2']} != {c_vals[1]}")
    
    def test_fourth_order_lowpass(self):
        order = 4
        c_vals = [1e-6,1e-9,1e-6,1e-9] # Résistances spécifiées
        r_vals = [58.06,213046.72,79.96,123079.84]# Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.lowpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[0], delta=self.tolerance, msg=f"Erreur pour C1 : {stages[0]['params']['C1']} != {c_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[1], delta=self.tolerance, msg=f"Erreur pour C2 : {stages[0]['params']['C2']} != {c_vals[1]}")
        self.assertAlmostEqual(stages[0]["params"]["C1"], c_vals[2], delta=self.tolerance, msg=f"Erreur pour C1 : {stages[0]['params']['C1']} != {c_vals[2]}")
       
        self.assertAlmostEqual(stages[0]["params"]["C2"], c_vals[3], delta=self.tolerance, msg=f"Erreur pour C2 : {stages[0]['params']['C2']} != {c_vals[3]}")  
    
    def test_first_order_highpass(self):
        order = 1
        c_vals = [1.59e-7]  # Résistances spécifiées
        r_vals = [159.15] # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.highpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals, c_vals=c_vals)

        self.assertAlmostEqual(stages[0]["params"]["R"], r_vals[0], delta=self.tolerance, msg=f"Erreur pour C : {stages[0]['params']['R']} != {r_vals[0]}")
    def test_second_order_highpass(self):
        order = 2
        c_vals = [1e-6,1e-9]  # Résistances spécifiées
        r_vals = [72.2,216725.55] # Capacités spécifiées

        # Résultats obtenus via la fonction
        _, stages = self.highpass_instance.components(order=order, cutoff_freq=self.cutoff_freq, r_vals=r_vals,c_vals=c_vals)

   
        self.assertAlmostEqual(stages[0]["params"]["R1"], r_vals[0], delta=self.tolerance, msg=f"Erreur pour R1 : {stages[0]['params']['R1']} != {r_vals[0]}")
        self.assertAlmostEqual(stages[0]["params"]["R2"], r_vals[1], delta=self.tolerance, msg=f"Erreur pour R2 : {stages[0]['params']['R2']} != {r_vals[1]}")

if __name__ == "__main__":
    unittest.main()
