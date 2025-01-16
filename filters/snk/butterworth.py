import numpy as np
import matplotlib as plt

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
        Cette fonction calcule les composants manquants pour réaliser le filtre voulu ainsi que ça fonction de transfert.

        Pour les ordres impairs, les R1,C1 seront toujours les valeurs de composants du premier étage (1er ordre).
        Pour les ordres pairs, les composants seront par paire de 2.
        '''
        # Verifie si l'ordre du filtre est dans le dictionnaire
        if order not in self.BUTTERWORTH_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")   

        if cutoff_frequency is None:
            raise ValueError("Veuillez fournir une fréquence de coupure.")
        pulsation_W0 = 2 * np.pi * cutoff_frequency

        if order == 1:
            if res_values is not None:
                # Verifie que l'entrée est du type int ou float
                if isinstance(res_values[0], (int, float)):
                    # Verifie que la liste ne contient qu'1 seul élément
                    nbr_elements = len(res_values)
                    if nbr_elements > 1:
                        raise ValueError("Pour le 1er ordre veuillez mettre qu'une seule résistance.")
                    condo_values = 1/(pulsation_W0 * res_values[0])
                    return {"R": res_values, "C": condo_values} 
                else:
                    raise ValueError("Veuillez insérer une valeur de résistance valable.")
            elif condo_values is not None:
                # Verifie que l'entrée est du type int ou float
                if isinstance(condo_values[0], (int,float)):
                    nbr_elements = len(condo_values)
                    # Verifie que la liste ne contient qu'1 seul élément
                    if nbr_elements > 1 or nbr_elements == 0:
                        raise ValueError("Pour le 1er ordre veuillez mettre qu'un seul élément.")
                    res_values = 1/(pulsation_W0 * condo_values[0])
                    return {"R": res_values, "C": condo_values[0]}
                else:
                    raise ValueError("Veuillez insérer une valeur valable.")
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")
        '''    
        if order >= 2:
            if res_values is not None:
                if isinstance(res_values, (int,float)):
                    nbr_elements = len(res_values)
                    if nbr_elements > order:
                        raise ValueError(f"Votre liste contient trop de résistances.\n Inséré {nbr_elements} résistances dans votre liste.") 
                    elif nbr_elements < order:
                        raise ValueError(f"Votre liste ne contient pas assez de résistances.\n Inséré {nbr_elements} résistances dans votre liste.")
                    else:
                        pulsation_W0 = 2 * np.pi * cutoff_frequency
                        quality_q0.clear()
                        # Verifie si l'ordre est paire
                        if order % 2 == 0:
                            nbr_etages = order/2
                            for x in nbr_etages:
                                quality_q0.append(self.BUTTERWORTH_TABLE[order][x - 1]) # Au cas ou plusieurs etages
                                c2 = 1 / ((res_values[0] + res_values[1]) * pulsation_W0 * quality_q0[0])
                                c1 = ((res_values[0] + res_values[1]) * quality_q0[0]) / (res_values[0] * res_values[1] * pulsation_W0)
                                return {"R1": res_values[0], "R2": res_values[1], "C1": c1, "C2":c2}
                        # Verifie si l'ordre est impaire
                        if order % 2 != 0:
                            nbr_etages = order // 2 + 1
                else:
                    raise ValueError("Veuillez insérer des valeurs valables.")'''