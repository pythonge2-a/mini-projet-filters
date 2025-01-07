import numpy as np
from scipy.signal import TransferFunction

class Butterworth_LowPass:
    def __init__(self):    
        # Tableau des pulsations et facteurs de qualité des filtres passe-bas normalisés
        self.BUTTERWORTH_TABLE = {
            1: [(0.0)],
            2: [(0.7071)],
            3: [(0.0), (1.0)],
            4: [(0.5412), (1.3066)],
            5: [(0.0), (0.6180), (1.6180)],
            6: [(0.5176), (0.7071), (1.9319)],
            7: [(0.0000), (0.5550), (0.8019), (2.2470)],
            8: [(0.5098), (0.6013), (0.8999), (2.5629)],
            9: [(0.0000), (0.5321), (0.6527), (1.0), (2.8794)],
           10: [(0.5062), (0.5612), (0.7071), (1.1013), (3.1962)],
        }

    def components(self, order, cutoff_frequency, res_values=None, condo_values=None):
        r'''
        Pour les composants inséré veuillez mettre des liste car a partir du 2ème ordre plusieurs éléments devront etre inséré
        '''
        # Verifie si l'ordre du filtre est entre 1 et 10
        if order not in self.BUTTERWORTH_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")    
        quality_Q0 = self.BUTTERWORTH_TABLE[order][0]
        if cutoff_frequency is None:
            raise ValueError("Veuillez fournir une fréquence de coupure.")
        pulsation_W0 = 2 * np.pi * cutoff_frequency
        if order is 1:
            if res_values is not None:
                if isinstance(res_values[0], (int, float)):
                    nbr_elements = len(res_values)
                    if nbr_elements > 2:
                        raise ValueError("Pour le 1er ordre veuillez mettre qu'une seule résistance.")
                    condo_values = 1/(pulsation_W0 * res_values[0])
                    num = [1]
                    den = [res_values[0] * condo_values, 1]
                    return TransferFunction(num, den), {"R": res_values[0], "C": condo_values} 
                else:
                    raise ValueError("Veuillez insérer une valeur de résistance valable.")
            elif condo_values is not None:
                if isinstance(condo_values[0], (int,float)):
                    nbr_elements = len(condo_values)
                    if nbr_elements > 2:
                        raise ValueError("Pour le 1er ordre veuillez mettre qu'un seul élément.")
                    res_values = 1/(pulsation_W0 * condo_values[0])
                    num = [1]
                    den = [condo_values[0] * res_values, 1]
                    return TransferFunction(num, den), {"R": res_values, "C": condo_values[0]}
                else:
                    raise ValueError("Veuillez insérer une valeur valable.")