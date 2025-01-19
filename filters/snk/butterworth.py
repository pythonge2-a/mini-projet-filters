import numpy as np
import matplotlib.pyplot as plt

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
        '''
        Cette fonction calcule les composants manquants pour réaliser le filtre Passe-Bas voulu.

        Si ordre impair alors premier élément de la liste doit être R1 du premier étage
        Format de la liste des résistances --> [R1, R2, ..., Rx]
        Format de la liste des condensateurs --> [C1, C2, ..., Cx]
        Format de sortie --> [Composant 1, Composant 2, ..., Composant x]
        '''
        # Verifie si l'ordre du filtre est dans le dictionnaire
        if order not in self.BUTTERWORTH_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")   

        # Verifie que une frequence de coupure à été donnée
        if cutoff_frequency is None:
            raise ValueError("Veuillez fournir une fréquence de coupure.")
        # Calcul de W0
        pulsation_W0 = 2 * np.pi * cutoff_frequency

        if order == 1:
            # Verifie si une liste à été entrée
            if res_values is not None:
                nbr_elements = len(res_values)
                # Verifie la taille de la liste
                if nbr_elements > 1 or nbr_elements == 0:
                    raise ValueError("Pour le 1er ordre veuillez mettre une résistance.")
                # Verifie le type des données de la liste
                if isinstance(res_values[0], (int, float)):
                    # Calcul du composant
                    condo_value = 1/(pulsation_W0 * res_values[0])
                    # Retourne les valeurs calculées
                    return {res_values[0], condo_value}
                else:
                    raise ValueError("Veuillez insérer des valeurs valables.")
            elif condo_values is not None:
                # Verifie que l'entrée est du type int ou float
                nbr_elements = len(condo_values)
                # Verifie que la liste ne contient qu'1 seul élément
                if nbr_elements > 1 or nbr_elements == 0:
                    raise ValueError("Pour le 1er ordre veuillez mettre qu'un seul élément.")
                if isinstance(condo_values[0], (int, float)):
                    res_value = 1/(pulsation_W0 * condo_values[0])
                    return {res_value, condo_values[0]}
                else:
                    raise ValueError("Veuillez insérer des int ou float.")
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")
            
        if order >= 2:
            # Verifie si une liste à été entrée
            if res_values is not None:
                nbr_elements = len(res_values)
                # Verifie si l'ordre est pair
                # Calcul le nombre d'étages du filtre
                if order % 2 == 0:
                    stage = order/2
                    # Verifie que la taille de la liste soit conforme
                    if nbr_elements < (stage * 2) or nbr_elements > (stage * 2):
                        raise ValueError(f"Veuillez mettre {stage * 2} éléments dans la liste.")
                    for x in res_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    condo_value = []
                    for y in range(int(stage)):
                        # Indices pour les résistances
                        res_start = y * 2
                        res_end = res_start + 2
                        current_res =  res_values[res_start:res_end]
                        # Calculs des composants
                        c = 1 / ((current_res[0] + current_res[1]) * pulsation_W0 * self.BUTTERWORTH_TABLE[order][y])
                        condo_value.append(c)

                        c = ((current_res[0] + current_res[1]) * self.BUTTERWORTH_TABLE[order][y]) / (current_res[0] * current_res[1] * pulsation_W0)
                        condo_value.append(c)

                # Verifie si l'ordre est impair
                elif order % 2 != 0:
                    stage = order // 2 + 1
                    # Verifie que la taille de la liste soit conforme 
                    if nbr_elements > (stage * 2 - 1) or nbr_elements < (stage * 2 - 1):
                        raise ValueError(f"Veuillez mettre {stage * 2 - 1} éléments dans la liste.")
                    for x in res_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    condo_value = []
                    for y in range(int(stage)):
                        # Indices pour les résistances
                        if y == 0: # Premier stage pour un ordre impair
                            # Le premier stage utilise une résistance et un condensateur
                            c1 = 1/(pulsation_W0 * res_values[0])
                            condo_value.append(c1)
                        else:
                            # Stages suivants avec deux résistances et deux condensateurs
                            res_start = (y - 1) * 2 + (1 if order % 2 == 1 else 0)
                            current_res = res_values[res_start:res_start + 2]
                            c1 = 1 / ((current_res[0] + current_res[1]) * pulsation_W0 * self.BUTTERWORTH_TABLE[order][y])
                            c2 = ((current_res[0] + current_res[1]) * self.BUTTERWORTH_TABLE[order][y]) / (current_res[0] * current_res[1] * pulsation_W0)
                            condo_value.extend([c1, c2])
                return condo_value
            
            # Verifie si une liste à été entrée
            if condo_values is not None:
                nbr_elements = len(condo_values)
                # Verifie si l'ordre est pair
                # Calcul le nombre d'étages du filtre
                if order % 2 == 0:  
                    stage = order/2
                    # Verifie que la taille de la liste soit conforme
                    if nbr_elements < (stage * 2) or nbr_elements > (stage * 2):
                        raise ValueError(f"Veuillez mettre {stage * 2} éléments dans la liste.")
                    for x in condo_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    res_value = []
                    for y in range(int(stage)):
                        c1_idx = y * 2
                        c2_idx = c1_idx + 1
                        q0 = self.BUTTERWORTH_TABLE[order][y]
                        c1, c2 = condo_values[c1_idx], condo_values[c2_idx]

                        # Vérification de la condition
                        if condo_values[c1_idx] < 4 * condo_values[c2_idx] * self.BUTTERWORTH_TABLE[order][y]**2:
                            raise ValueError(f"Condition non respectée au stage {y + 1}: C1 ({c1}) >= 4 * C2 ({c2}) * Q0^2 ({q0**2}).")

                        # Calcul des résistances
                        j = 1 / (pulsation_W0 * condo_values[c2_idx] * q0)
                        k = 1 / (pulsation_W0**2 * condo_values[c1_idx] * condo_values[c2_idx])
                        discriminant = j**2 - 4 * k

                        if discriminant < 0:
                            raise ValueError("Discriminant négatif : impossible de calculer les résistances.")

                        r2 = (j + np.sqrt(discriminant)) / 2
                        r1 = j - r2
                        res_value.extend([r1, r2])
                # Verifie si l'ordre est impair
                # Calcul le nombre d'étages du filtre
                elif order % 2 != 0:
                    stage = order // 2 + 1
                    # Verifie que la taille de la liste soit conforme 
                    if nbr_elements > (stage * 2 - 1) or nbr_elements < (stage * 2 - 1):
                        raise ValueError(f"Veuillez mettre {stage * 2 - 1} éléments dans la liste.")
                    for x in condo_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    res_value = []
                    # Calcul des étages d'ordre 2 suivants
                    for y in range(int(stage)):
                        c1_idx = y * 2 - 1
                        c2_idx = c1_idx + 1
                        q0 = self.BUTTERWORTH_TABLE[order][y]
                        c1, c2 = condo_values[c1_idx], condo_values[c2_idx]
                        # Calcul des résistances
                        if y == 0:
                            # Calcul du premier étage
                            r1 = 1/(pulsation_W0 * condo_values[0])
                            res_value.append(r1)
                        if y > 0:
                            # Vérification de la condition
                            if condo_values[c1_idx] < 4 * condo_values[c2_idx] * self.BUTTERWORTH_TABLE[order][y]**2:
                                raise ValueError(f"Condition non respectée au stage {y + 1}: C1 ({c1}) >= 4 * C2 ({c2}) * Q0^2 ({q0**2}).")
                            j = 1 / (pulsation_W0 * condo_values[c2_idx] * q0)
                            k = 1 / (pulsation_W0**2 * condo_values[c1_idx] * condo_values[c2_idx])
                            discriminant = j**2 - 4 * k

                            if discriminant < 0:
                                raise ValueError("Discriminant négatif : impossible de calculer les résistances.")

                            r2 = (j + np.sqrt(discriminant)) / 2
                            r1 = j - r2
                            res_value.extend([r1, r2])
                return res_value
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")

class Butterworth_HighPass:
    def __init__(self):    
        # Tableau des pulsations et facteurs de qualité des filtres passe-bas normalisés (Pour Butterworth rien ne change)
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
        '''
        Cette fonction calcule les composants manquants pour réaliser le filtre Passe-Haut voulu.

        Si ordre impair alors premier élément de la liste doit être R1 du premier étage
        Format de la liste des résistances --> [R1, R2, ..., Rx]
        Format de la liste des condensateurs --> [C1, C2, ..., Cx]
        Format de sortie --> [Composant 1, Composant 2, ..., Composant x]
        '''
        # Verifie si l'ordre du filtre est dans le dictionnaire
        if order not in self.BUTTERWORTH_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")   

        # Verifie que une frequence de coupure à été donnée
        if cutoff_frequency is None:
            raise ValueError("Veuillez fournir une fréquence de coupure.")
        # Calcul de W0
        pulsation_W0 = 2 * np.pi * cutoff_frequency

        if order == 1:
            # Verifie si une liste à été entrée
            if res_values is not None:
                nbr_elements = len(res_values)
                # Verifie la taille de la liste
                if nbr_elements > 1 or nbr_elements == 0:
                    raise ValueError("Pour le 1er ordre veuillez mettre une résistance.")
                # Verifie le type des données de la liste
                if isinstance(res_values[0], (int, float)):
                    # Calcul du composant
                    condo_value = 1/(pulsation_W0 * res_values[0])
                    # Retourne les valeurs calculées
                    return {res_values[0], condo_value}
                else:
                    raise ValueError("Veuillez insérer des valeurs valables.")
            elif condo_values is not None:
                # Verifie que l'entrée est du type int ou float
                nbr_elements = len(condo_values)
                # Verifie que la liste ne contient qu'1 seul élément
                if nbr_elements > 1 or nbr_elements == 0:
                    raise ValueError("Pour le 1er ordre veuillez mettre qu'un seul élément.")
                if isinstance(condo_values[0], (int, float)):
                    res_value = 1/(pulsation_W0 * condo_values[0])
                    return {res_value, condo_values[0]}
                else:
                    raise ValueError("Veuillez insérer des int ou float.")
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")
            
        if order >= 2:
            # Verifie si une liste à été entrée
            if res_values is not None:
                nbr_elements = len(res_values)
                # Verifie si l'ordre est pair
                # Calcul le nombre d'étages du filtre
                if order % 2 == 0:
                    stage = order/2
                    # Verifie que la taille de la liste soit conforme
                    if nbr_elements < (stage * 2) or nbr_elements > (stage * 2):
                        raise ValueError(f"Veuillez mettre {stage * 2} éléments dans la liste.")
                    for x in res_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    condo_value = []
                    for y in range(int(stage)):
                        r1_idx = y * 2
                        r2_idx = r1_idx + 1
                        q0 = self.BUTTERWORTH_TABLE[order][y]
                        r1, r2 = res_values[r1_idx], res_values[r2_idx]

                        # Vérification de la condition
                        if res_values[r2_idx] < 4 * res_values[r1_idx] * self.BUTTERWORTH_TABLE[order][y]**2:
                            raise ValueError(f"Condition non respectée au stage {y + 1}: R2 ({c1}) >= 4 * R1 ({c2}) * Q0^2 ({q0**2}).")

                        # Calcul des résistances
                        j = 1 / (pulsation_W0 * res_values[r1_idx] * q0)
                        k = 1 / (pulsation_W0**2 * res_values[r1_idx] * res_values[r2_idx])
                        discriminant = j**2 - 4 * k

                        if discriminant < 0:
                            raise ValueError("Discriminant négatif : impossible de calculer les résistances.")

                        c2 = (j + np.sqrt(discriminant)) / 2
                        c1 = j - r2
                        condo_value.extend([c1, c2])

                # Verifie si l'ordre est impair
                elif order % 2 != 0:
                    stage = order // 2 + 1
                    # Verifie que la taille de la liste soit conforme 
                    if nbr_elements > (stage * 2 - 1) or nbr_elements < (stage * 2 - 1):
                        raise ValueError(f"Veuillez mettre {stage * 2 - 1} éléments dans la liste.")
                    for x in res_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    condo_value = []
                    # Calcul des étages d'ordre 2 suivants
                    for y in range(int(stage)):
                        r1_idx = y * 2 - 1
                        r2_idx = r1_idx + 1
                        q0 = self.BUTTERWORTH_TABLE[order][y]
                        r1, r2 = res_values[r1_idx], res_values[r2_idx]
                        # Calcul des résistances
                        if y == 0:
                            # Calcul du premier étage
                            c1 = 1/(pulsation_W0 * res_values[0])
                            condo_value.append(c1)
                        if y > 0:
                            # Vérification de la condition
                            if res_values[r2_idx] < 4 * res_values[r1_idx] * self.BUTTERWORTH_TABLE[order][y]**2:
                                raise ValueError(f"Condition non respectée au stage {y + 1}: R2 ({c1}) >= 4 * R1 ({c2}) * Q0^2 ({q0**2}).")
                            j = 1 / (pulsation_W0 * res_values[r1_idx] * q0)
                            k = 1 / (pulsation_W0**2 * res_values[r1_idx] * res_values[r2_idx])
                            discriminant = j**2 - 4 * k

                            if discriminant < 0:
                                raise ValueError("Discriminant négatif : impossible de calculer les résistances.")

                            c2 = (j + np.sqrt(discriminant)) / 2
                            c1 = j - c2
                            condo_value.extend([c1, c2])
                return condo_value
            
            # Verifie si une liste à été entrée
            if condo_values is not None:
                nbr_elements = len(condo_values)
                # Verifie si l'ordre est pair
                # Calcul le nombre d'étages du filtre
                if order % 2 == 0:  
                    stage = order/2
                    # Verifie que la taille de la liste soit conforme
                    if nbr_elements < (stage * 2) or nbr_elements > (stage * 2):
                        raise ValueError(f"Veuillez mettre {stage * 2} éléments dans la liste.")
                    for x in condo_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    res_value = []
                    for y in range(int(stage)):
                        c1_idx = y * 2
                        c2_idx = c1_idx + 1
                        q0 = self.BUTTERWORTH_TABLE[order][y]
                        c1, c2 = condo_values[c1_idx], condo_values[c2_idx]

                        # Calcul des résistances
                        r2 = ((c1 + c2) * self.BUTTERWORTH_TABLE[order][y]) / (c1 * c2 * pulsation_W0)
                        r1 = 1/((c1 + c2) * pulsation_W0 * self.BUTTERWORTH_TABLE[order][y])
                        res_value.extend([r1, r2])
                # Verifie si l'ordre est impair
                # Calcul le nombre d'étages du filtre
                elif order % 2 != 0:
                    stage = order // 2 + 1
                    # Verifie que la taille de la liste soit conforme 
                    if nbr_elements > (stage * 2 - 1) or nbr_elements < (stage * 2 - 1):
                        raise ValueError(f"Veuillez mettre {stage * 2 - 1} éléments dans la liste.")
                    for x in condo_values:
                        # Verifier que chaque élément est du type int ou float
                        if isinstance(x, (int, float)):
                            None
                        else:
                            raise ValueError("Veuillez insérer des int ou float.")
                    res_value = []
                    # Calcul des étages d'ordre 2 suivants
                    for y in range(int(stage)):
                        c1_idx = y * 2 - 1
                        c2_idx = c1_idx + 1
                        q0 = self.BUTTERWORTH_TABLE[order][y]
                        c1, c2 = condo_values[c1_idx], condo_values[c2_idx]
                        # Calcul des résistances
                        if y == 0:
                            # Calcul du premier étage
                            r1 = 1/(pulsation_W0 * condo_values[0])
                            res_value.append(r1)
                        if y > 0:
                            # Vérification de la condition
                            r2 = ((c1 + c2) * self.BUTTERWORTH_TABLE[order][y]) / (c1 * c2 * pulsation_W0)
                            r1 = 1/((c1 + c2) * pulsation_W0 * self.BUTTERWORTH_TABLE[order][y])
                            res_value.extend([r1, r2])
                    return res_value
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")