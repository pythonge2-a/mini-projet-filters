import unittest
from filters.snk.butterworth import Butterworth_HighPass,Butterworth_LowPass

class TestButterWorthFilters(unittest.TestCase):
    def setUp(self):
        self.cutoff_frequency = 1000  # Fréquence de coupure en Hz
        self.tolerance_resistance = 0.1  # Tolérance pour les comparaisons des resistances
        self.tolerance_condensateur = 1e-6  # Tolérance pour les comparaisons des condensateurs
        self.lowpass_instance = Butterworth_LowPass()  # Créer une instance de `lowpass`
        self.highpass_instance = Butterworth_HighPass()

    def test_first_order_lowpass_C(self):   # Test pour calcul des condensateurs
        order = 1   # Ordre du filtre
        res_values = [1500]  # Valeur donnée pour la résistance
        condo_wanted = 1.061032953945969e-07 # Valeur attendue à la sortie de la fonction

        # Appel de la méthode avec l'instance
        values = self.lowpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        value_test = values["C"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            value_test, condo_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condo_wanted} != {value_test}"
        )

    def test_first_order_lowpass_R(self):   # Test pour calcul des resistances
        order = 1   # Ordre du filtre
        condo_values = [2.2e-6]  # Valeur donnée pour le condensateur
        res_wanted = 72.34315595086152 # Valeur attendue à la sortie de la fonction

        # Appel de la méthode avec l'instance
        values = self.lowpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, condo_values=condo_values)
        value_test = values["R"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            value_test, res_wanted, delta=self.tolerance_resistance,
            msg=f"Erreur pour C : {res_wanted} != {value_test}"
        )
        
    def test_second_order_lowpass(self):    # Test pour calcul des condensateurs
        order = 2
        res_values = [1000, 5000]  # Résistances spécifiées
        condos_wanted = [3.751353959644919e-08, 1.3504615231233503e-07] # Valeur attendue à la sortie de la fonction

        # Résultats obtenus via la fonction
        values = self.lowpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        values_test = values["C"]
        
        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            values_test, condos_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condos_wanted} != {values_test}"
        )

    def test_third_order_lowpass(self): # Test pour calcul des condensateurs
        order = 3
        res_values = [1000, 5000, 12000]  # Résistances spécifiées
        condos_wanted = [1.5915494309189535e-07, 9.362055475993845e-09, 4.5093900542703685e-08] # Valeur attendue à la sortie de la fonction

        # Résultats obtenus via la fonction
        values = self.lowpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        values_test = values["C"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            values_test, condos_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condos_wanted} != {values_test}"
        )

    def test_fourth_order_lowpass(self):  # Test pour calcul des condensateurs
        order = 4
        res_values = [1000, 5000, 12000, 6000]  # Résistances spécifiées
        condos_wanted = [4.901297828649154e-08, 1.0336158624160052e-07, 6.767137060219711e-09, 5.1987962160967617e-08] # Valeur attendue à la sortie de la fonction

        # Résultats obtenus via la fonction
        values = self.lowpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        values_test = values["C"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            values_test, condos_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condos_wanted} != {values_test}"
        )
            
    def test_first_order_highpass_C(self):   # Test pour calcul des condensateurs
        order = 1   # Ordre du filtre
        res_values = [1500]  # Valeur donnée pour la résistance
        condo_wanted = 1.061032953945969e-07 # Valeur attendue à la sortie de la fonction

        # Appel de la méthode avec l'instance
        values = self.highpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        value_test = values["C"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            value_test, condo_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condo_wanted} != {value_test}"
        )

    def test_first_order_highpass_R(self):   # Test pour calcul des resistances
        order = 1   # Ordre du filtre
        condo_values = [2.2e-6]  # Valeur donnée pour le condensateur
        res_wanted = 72.34315595086152 # Valeur attendue à la sortie de la fonction

        # Appel de la méthode avec l'instance
        values = self.highpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, condo_values=condo_values)
        value_test = values["R"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            value_test, res_wanted, delta=self.tolerance_resistance,
            msg=f"Erreur pour C : {res_wanted} != {value_test}"
        )
        
    def test_second_order_highpass(self):    # Test pour calcul des condensateurs
        order = 2
        res_values = [1000, 5000]  # Résistances spécifiées
        condos_wanted = [2.5366472992716085e-08, 1.997147645859791e-07] # Valeur attendue à la sortie de la fonction

        # Résultats obtenus via la fonction
        values = self.highpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        values_test = values["C"]
        
        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            values_test, condos_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condos_wanted} != {values_test}"
        )

    def test_third_order_highpass_invalid_resistances(self): # Test pour une erreur au cas ou la condition n'est pas respecter
        order = 3
        res_values = [1000, 5000, 12000]  # Résistances spécifiées
        condos_wanted = ValueError  # Valeur attendue à la sortie de la fonction

        # Vérification qu'un ValueError est levé
        with self.assertRaises(ValueError) as context:
            self.highpass_instance.components(
                order=order,
                cutoff_frequency=self.cutoff_frequency,
                res_values=res_values
            )

        # Vérification du message d'erreur (optionnel)
        self.assertIn("Condition non respectée au stage", str(context.exception))

    def test_fourth_order_highpass(self):  # Test pour calcul des condensateurs
        order = 4
        res_values = [1000, 5000, 1000, 12000]  # Résistances spécifiées
        condos_wanted = [1.837507364548829e-08, 2.75702796073461e-07, 2.0923392245338568e-08, 1.0088507483861624e-07] # Valeur attendue à la sortie de la fonction

        # Résultats obtenus via la fonction
        values = self.highpass_instance.components(order=order, cutoff_frequency=self.cutoff_frequency, res_values=res_values)
        values_test = values["C"]

        # Comparation des valeurs attendues avec obtenues avec une certaine tolérance
        self.assertAlmostEqual(
            values_test, condos_wanted, delta=self.tolerance_condensateur,
            msg=f"Erreur pour C : {condos_wanted} != {values_test}"
        )

if __name__ == "__main__":
    unittest.main()
