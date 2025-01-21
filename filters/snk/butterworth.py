import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import TransferFunction, bode

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
                    return {"R": res_values[0], "C": condo_value}
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
                    return {"R": res_value, "C": condo_values[0]}
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
                return {"R": res_values, "C": condo_value}
            
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
                return {"R": res_value, "C": condo_values}
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")

    def graphs(self, order, cutoff_frequency=None,res_values=None, condo_values=None):
        # Calcul des Fct de transfert des differents etages
        if order > 10 or order < 1:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")   
        num_combined = []
        den_combined = []
        if cutoff_frequency is not None:    # Si une frequence de coupure est donnée
            # *Pour avoir un graphe qui est toujours dans les bonnes plages
            freq_min_hz = cutoff_frequency / 10_000  # 10^(-4) fois la fréquence de coupure
            freq_max_hz = cutoff_frequency * 10_000  # 10^(4) fois la fréquence de coupure
            if res_values is not None:
                condo = self.components(order=order, cutoff_frequency=cutoff_frequency, res_values=res_values)
                if order == 1:
                    num = [1.0]
                    den = [res_values[0] * condo["C"], 1.0]
                    tf = TransferFunction(num, den)
                    # *Pour avoir un graphe qui est toujours dans les bonnes plages
                    w = np.logspace(np.log10(2 * np.pi * freq_min_hz), np.log10(2 * np.pi * freq_max_hz), 500)
                    w, mag, phase = bode(tf, w=w)
                    freq_hz = w/(2 * np.pi)
                    fig_lp, (ax_mag_lp, ax_phase_lp) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

                    # Tracer les données
                    ax_mag_lp.semilogx(freq_hz, mag, 'b')
                    ax_mag_lp.set_ylabel('Magnitude (dB)')
                    ax_phase_lp.semilogx(freq_hz, phase, 'r')
                    ax_phase_lp.set_xlabel('Frequency (Hz)')
                    ax_phase_lp.set_ylabel('Phase (deg)')

                    # Activer les grilles pour les deux axes
                    ax_mag_lp.grid(True, which='both', axis='x')  # Ajout des lignes de grille entre les décennies
                    ax_phase_lp.grid(True, which='both', axis='x')
                    # Définir les ticks sur l'axe x pour que les grilles apparaissent entre chaque décennie
                    ax_mag_lp.set_xticks(np.logspace(np.log10(min(freq_hz)), np.log10(max(freq_hz)), num=9))  # ajuster les décades selon la plage de fréquence
                    ax_phase_lp.set_xticks(np.logspace(np.log10(min(freq_hz)), np.log10(max(freq_hz)), num=9))

                    # Ajuster la disposition pour éviter le chevauchement
                    plt.tight_layout()
                    # Afficher le graphique
                    plt.show()
            elif condo_values is not None:
                res = self.components(order=order, cutoff_frequency=cutoff_frequency, condo_values=condo_values)
                if order == 1:
                    num = [1.0]
                    den = [condo_values[0] * res["R"], 1.0]
                    tf = TransferFunction(num, den)
                    # *Pour avoir un graphe qui est toujours dans les bonnes plages
                    w = np.logspace(np.log10(2 * np.pi * freq_min_hz), np.log10(2 * np.pi * freq_max_hz), 500)
                    _, mag, phase = bode(tf, w=w)
                    freq_hz = w/(2 * np.pi)
                    fig_lp, (ax_mag_lp, ax_phase_lp) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

                    # Tracer les données
                    ax_mag_lp.semilogx(freq_hz, mag, 'b')
                    ax_mag_lp.set_ylabel('Magnitude (dB)')
                    ax_phase_lp.semilogx(freq_hz, phase, 'r')
                    ax_phase_lp.set_xlabel('Frequency (Hz)')
                    ax_phase_lp.set_ylabel('Phase (deg)')

                    # Activer les grilles pour les deux axes
                    ax_mag_lp.grid(True, which='both', axis='x')  # Ajout des lignes de grille entre les décennies
                    ax_phase_lp.grid(True, which='both', axis='x')
                    # Définir les ticks sur l'axe x pour que les grilles apparaissent entre chaque décennie
                    ax_mag_lp.set_xticks(np.logspace(np.log10(min(freq_hz)), np.log10(max(freq_hz)), num=9))  # ajuster les décades selon la plage de fréquence
                    ax_phase_lp.set_xticks(np.logspace(np.log10(min(freq_hz)), np.log10(max(freq_hz)), num=9))

                    # Ajuster la disposition pour éviter le chevauchement
                    plt.tight_layout()
                    # Afficher le graphique
                    plt.show()
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")
        elif res_values is not None and condo_values is not None:
            None
        else:
            raise KeyError("Veuillez inserer les listes de composants manquantes, " + 
                           "Ou ajouter une frequence de coupure.")
        
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
                    return {"R": res_values[0], "C": condo_value}
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
                    return {"R": res_value, "C": condo_values[0]}
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
                            raise ValueError(f"Condition non respectée au stage {y + 1}: R2 ({r2}) >= 4 * R1 ({r1}) * Q0^2 ({q0**2}).")

                        # Calcul des résistances
                        j = 1 / (pulsation_W0 * res_values[r1_idx] * q0)
                        k = 1 / (pulsation_W0**2 * res_values[r1_idx] * res_values[r2_idx])
                        discriminant = j**2 - 4 * k

                        if discriminant < 0:
                            raise ValueError("Discriminant négatif : impossible de calculer les résistances.")

                        c2 = (j + np.sqrt(discriminant)) / 2
                        c1 = j - c2
                        condo_value.extend([c1, c2])
                    return {"R": res_values, "C": condo_value}

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
                                raise ValueError(f"Condition non respectée au stage {y + 1}: R2 ({r2}) >= 4 * R1 ({r1}) * Q0^2 ({q0**2}).")
                            j = 1 / (pulsation_W0 * res_values[r1_idx] * q0)
                            k = 1 / (pulsation_W0**2 * res_values[r1_idx] * res_values[r2_idx])
                            discriminant = j**2 - 4 * k

                            if discriminant < 0:
                                raise ValueError("Discriminant négatif : impossible de calculer les résistances.")

                            c2 = (j + np.sqrt(discriminant)) / 2
                            c1 = j - c2
                            condo_value.extend([c1, c2])
                    return {"R": res_values, "C": condo_value}
            # Verifie si une liste à été entrée
            elif condo_values is not None:
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
                    return {"R": res_value, "C": condo_values}
            else:
                raise KeyError("Veuillez au moin insérer une liste de composants.")
            
    def graphs(self, order, cutoff_freq, r_vals, c_vals):
        # Calcul des Fct de transfert des differents etages
        if order > 10 or order < 1:
            raise ValueError("Ordre de filtre non Supporté, MAX : 10.")
        pulsation_w0 = 2 * np.pi * cutoff_freq
        num_stage = []
        den_stage = []

        if order % 2 == 0:
            stage = order / 2
            for x in range(int(stage)):
                r1_idx = x * 2
                r2_idx = r1_idx + 1
                c1_idx = x * 2
                c2_idx = r1_idx + 1
                num = r_vals[r1_idx] * r_vals[r2_idx] * c_vals[c1_idx] * c_vals[c2_idx]
                num_stage.append(num)
                den = 1 + c_vals[c2_idx] * (r_vals[c1_idx] + r_vals[r2_idx]) + c_vals[c1_idx] * c_vals[c2_idx] * r_vals[r1_idx] * r_vals[r2_idx]
                den_stage.append(den)
            combined_num = np.polymul(num_stage)
            combined_den = np.polymul(den_stage)
            combined_tf = TransferFunction(combined_num, combined_den)
        elif order % 2 != 0:
            stage = (order // 2) + 1
            # ordre paire
            for x in range(int(stage)):
                r1_idx = x * 2 - 1
                r2_idx = r1_idx + 1
                c1_idx = x * 2 - 1
                c2_idx = r1_idx + 1
                if x == 0:
                    num_stage1 = r_vals[x] * c_vals[x]
                    num_stage.append(num_stage1, 1)
                    den_stage1 = 1 + r_vals[x] * c_vals[x]
                    den_stage.append(den_stage1, 1)
                else:    
                    num = r_vals[r1_idx] * r_vals[r2_idx] * c_vals[c1_idx] * c_vals[c2_idx]
                    num_stage.append(num)
                    den = 1 + c_vals[c2_idx] * (r_vals[c1_idx] + r_vals[r2_idx]) + c_vals[c1_idx] * c_vals[c2_idx] * r_vals[r1_idx] * r_vals[r2_idx]
                    den_stage.append(den)
            combined_num = np.polymul(num_stage, num_stage1)
            combined_den = np.polymul(den_stage, den_stage1)
            combined_tf = TransferFunction(combined_num, combined_den)

        # Génère une plage de fréquences logarithmiques entre 10^2 et 10^6 rad/s avec 500 points
        w2 = np.logspace(2, 6, 500)

        # Calcule la réponse en fréquence (magnitude et phase) pour la fonction de transfert 'combined_tf'
        # à l'aide des fréquences définies dans 'w2'
        w2, mag2, phase2 = bode(combined_tf, w=w2)
        # Convertit les fréquences angulaires (rad/s) en fréquences linéaires (Hz)
        freq_hz2 = w2 / (2 * np.pi)

        # Crée une figure avec deux sous-graphes empilés pour la magnitude et la phase
        fig_hp, (ax_mag_hp, ax_phase_hp) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

        # Trace la magnitude (en dB) en fonction de la fréquence (Hz) sur une échelle logarithmique
        ax_mag_hp.semilogx(freq_hz2, mag2, 'b')  # Courbe en bleu
        ax_mag_hp.set_ylabel('Magnitude (dB)')    # Libellé de l'axe Y pour la magnitude

        # Trace la phase (en degrés) en fonction de la fréquence (Hz) sur une échelle logarithmique
        ax_phase_hp.semilogx(freq_hz2, phase2, 'r')  # Courbe en rouge
        ax_phase_hp.set_xlabel('Frequency (Hz)')     # Libellé de l'axe X pour la fréquence
        ax_phase_hp.set_ylabel('Phase (deg)')        # Libellé de l'axe Y pour la phase

        # Ajoute une grille aux deux sous-graphes pour une meilleure lisibilité
        ax_phase_hp.grid(True)
        ax_mag_hp.grid(True)

        # Ajuste automatiquement les espaces entre les sous-graphes pour éviter les chevauchements
        plt.tight_layout()

        # Affiche le graphique
        plt.show()

'''     w = np.logspace(2, 5, 400)
        w, mag, phase = bode(combined_tf, w=w)
        freq_hz = w/(2*np.pi)

        fig_lp, (ax_mag_lp, ax_phase_lp) = plt.subplots(2,1,figsize=(8,6), sharex=True)
        ax_mag_lp.semilogx(freq_hz, mag, 'b')
        ax_mag_lp.set_ylabel('Magnitude (dB)')

        ax_phase_lp.semilogx(freq_hz, phase, 'r')
        ax_phase_lp.set_xlabel('Frequency (Hz)')
        ax_phase_lp.set_ylabel('Phase (deg)')
        ax_phase_lp.grid(True)
        ax_mag_lp.grid(True)

        plt.tight_layout()
        plt.show()'''