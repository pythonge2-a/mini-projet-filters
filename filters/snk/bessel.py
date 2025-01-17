import numpy as np
from scipy.signal import TransferFunction,bode
import matplotlib.pyplot as plt
class lowpass:
# Initialisation de la classe Lowpass avec les pôles de Bessel.
    def __init__(self):    
        self.BESSEL_TABLE = {
            1: [(1.0, 0.0)],
            2: [(1.2723, 0.577)],
            3: [(1.3225, 0.0), (1.4474, 0.691)],
            4: [(1.431, 0.5219), (1.6043, 0.8055)],
            5: [(1.5015, 0.0000), (1.5555, 0.5635), (1.7545, 0.9165)],
            6: [(1.6030, 0.5103), (1.6882, 0.6112), (1.9037, 1.0233)],
            7: [(1.6840, 0.0000), (1.7160, 0.5324), (1.8221, 0.6608), (2.0491, 1.0233)],
            8: [(1.7772, 0.5060), (1.8308, 0.5596), (1.9518, 0.7109), (2.1872, 1.2257)],
            9: [(1.8570, 0.0000), (1.8788, 0.5197), (1.9483, 0.5895), (2.0808, 0.7606), (2.3228, 1.3219)],
           10: [(1.9412, 0.5039), (1.9790, 0.5376), (2.0606, 0.6205), (2.2021, 0.8098), (2.4487, 1.4153)],
        }
# Retourne les valeurs de pulsation normalisée et de facteur de qualité pour un ordre donné.
    def bessel_q0_omega0(self, order):
        if order not in self.BESSEL_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")
        return self.BESSEL_TABLE[order]
# Calcule un filtre passe-bas de premier ordre.
    def first_order_lowpass(self, cutoff_freq, r=None, c=None, omega0_norm=None):
        
        if omega0_norm is None:
            omega0_norm = 1.0
        omega0 = 2 * np.pi * cutoff_freq * omega0_norm

        if r is not None:
            c = 1 / (omega0 * r)
        elif c is not None:
            r = 1 / (omega0 * c)
        else:
            raise ValueError("Fournir R ou C pour le calcul du premier ordre.")

        num = [1.0]
        den = [r * c,1.0]
        return TransferFunction(num, den), {"R": r, "C": c}
 # Calcule un filtre passe-bas de second ordre avec la cellule Sallen-Key.
    def sallen_key_lowpass(self, order, cutoff_freq, r1=None, r2=None, c1=None, c2=None, omega0_norm=None, q0=None):
        if order != 2:
            raise ValueError("Seul l'ordre 2 est supporté pour Sallen-Key.")

        if omega0_norm is None or q0 is None:
            raise ValueError("Omega0 et Q0 doivent être fournis.")

        omega0 = 2 * np.pi * cutoff_freq * omega0_norm

        if c1 is not None and c2 is not None:
            r1_plus_r2 = 1 / (omega0 * c2 * q0)
            r1_times_r2 = 1 / (omega0**2 * c1 * c2)
            a = 1
            b = -r1_plus_r2
            c = r1_times_r2
            discriminant = b**2 - 4 * a*c
            if discriminant < 0:
                raise ValueError("Les paramètres fournis pour C1 et C2 ne permettent pas un calcul valide de R1 et R2.")
            r21 = (-b + np.sqrt(discriminant)) / 2
           # r22 = (-b - np.sqrt(discriminant)) / 2
            r11 = r1_plus_r2 - r21
           # r12 = r1_plus_r2 - r22

            num = [1.0]
            den = [r11 * r21 * c1 * c2,(r11 + r21) * c2,1.0]
           # den2 = [r12 * r22 * c1 * c2,(r12 + r22) * c2,1.0]

            tf = TransferFunction(num, den)
          #  tf2 = TransferFunction(num, den2)

            return [
                {"tf": tf, "params": {"R1": r11, "R2": r21, "C1": c1, "C2": c2}},
                #{"tf": tf2, "params": {"R1": r12, "R2": r22, "C1": c1, "C2": c2}},
            ]
        elif r1 is not None and r2 is not None:
            c2 = 1 / (omega0 * q0 * (r1 + r2))
            c1 = 1 / (omega0**2 * r1 * r2 * c2)
            if c1 <= 0 or c2 <= 0:
                raise ValueError("Les paramètres fournis pour R1 et R2 ne permettent pas un calcul valide de C1 et C2.")

            num = [1.0]
            den = [r1 * r2 * c1 * c2, (r1 + r2) * c2,1.0]

            tf = TransferFunction(num, den)
            return [{"tf": tf, "params": {"R1": r1, "R2": r2, "C1": c1, "C2": c2}}]
        else:
            raise ValueError("Veuillez fournir soit (C1, C2), soit (R1, R2).")
        
# Calcule un filtre passe-bas de n'importe quel ordre en utilisant des cellules en cascade.
# Permet de choisir entre spécifier les résistances ou les condensateurs.
    def components(self, order, cutoff_freq, r_vals=None, c_vals=None):
        if order not in self.BESSEL_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")

        poles = self.bessel_q0_omega0(order)
        stages = []
        num_combined, den_combined = [1], [1]

        # Calculer le nombre d'éléments nécessaires
        num_stages = order // 2 + (order % 2)
        num_elements = 2 * num_stages

        if r_vals is None and c_vals is None:
            raise ValueError("Veuillez fournir soit les résistances (r_vals), soit les condensateurs (c_vals).")

        if r_vals is None:
            r_vals = [None] * num_elements
        if c_vals is None:
            c_vals = [None] * num_elements

        # Vérification de la longueur des listes
        if len(r_vals) < num_elements:
            r_vals.extend([None] * (num_elements - len(r_vals)))
        if len(c_vals) < num_elements:
            c_vals.extend([None] * (num_elements - len(c_vals)))

        # Parcourir les pôles pour construire les étapes
        for i, (omega0_norm, q0) in enumerate(poles):
            if q0 == 0.0:  # Premier ordre
                tf, params = self.first_order_lowpass(
                    cutoff_freq, r=r_vals[i], c=c_vals[i], omega0_norm=omega0_norm
                )
            else:  # Deuxième ordre
                idx = 2 * i
                tf_data = self.sallen_key_lowpass(
                    order=2,
                    cutoff_freq=cutoff_freq,
                    r1=r_vals[idx],
                    r2=r_vals[idx + 1],
                    c1=c_vals[idx],
                    c2=c_vals[idx + 1],
                    omega0_norm=omega0_norm,
                    q0=q0,
                )
                tf = tf_data[0]["tf"]
                params = tf_data[0]["params"]

            # Mise à jour de la fonction de transfert combinée
            num_combined = np.polymul(num_combined, tf.num)
            den_combined = np.polymul(den_combined, tf.den)

            # Ajouter les informations de l'étape
            stages.append({"tf": tf, "params": params})

        # Fonction de transfert combinée
        combined_tf = TransferFunction(num_combined, den_combined)
        return combined_tf, stages
        # Affiche le diagramme de Bode pour un filtre donné.
    def graphs(self, order, cutoff_freq, r_vals=None, c_vals=None):
        combined_tf, _ = self.components(order, cutoff_freq, r_vals, c_vals)
        w = np.logspace(2, 5, 400)
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
        plt.show()

class highpass:
    # Initialisation de la classe Lowpass avec les pôles de Bessel.
    def __init__(self):
        
        self.BESSEL_TABLE = {
            1: [(1.0, 0.0)],
            2: [(1.2723, 0.577)],
            3: [(1.3225, 0.0), (1.4474, 0.691)],
            4: [(1.431, 0.5219), (1.6043, 0.8055)],
            5: [(1.5015, 0.0000), (1.5555, 0.5635), (1.7545, 0.9165)],
            6: [(1.6030, 0.5103), (1.6882, 0.6112), (1.9037, 1.0233)],
            7: [(1.6840, 0.0000), (1.7160, 0.5324), (1.8221, 0.6608), (2.0491, 1.0233)],
            8: [(1.7772, 0.5060), (1.8308, 0.5596), (1.9518, 0.7109), (2.1872, 1.2257)],
            9: [(1.8570, 0.0000), (1.8788, 0.5197), (1.9483, 0.5895), (2.0808, 0.7606), (2.3228, 1.3219)],
           10: [(1.9412, 0.5039), (1.9790, 0.5376), (2.0606, 0.6205), (2.2021, 0.8098), (2.4487, 1.4153)],
        }
# Retourne les valeurs de pulsation normalisée et de facteur de qualité pour un ordre donné.
    def bessel_q0_omega0(self, order):
        
        if order not in self.BESSEL_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")
        return self.BESSEL_TABLE[order]  # Retourne tous les pôles
    
    def first_order_highpass(self, cutoff_freq, r=None, c=None, omega0_norm=None):
        # Calcule un filtre du premier ordre.
        omega0_norm_hp = 1 / omega0_norm
        omega_c = 2 * np.pi * cutoff_freq * omega0_norm_hp

        if r is not None:
            c = 1 / (omega_c * r)
        elif c is not None:
            r = 1 / (omega_c * c)
        else:
            raise ValueError("Fournir R ou C.")
        re = r*c
        num = [r*c, 0]
        den = [r*c, 1.0]
        print(den)

        return TransferFunction(num, den), {"R": r, "C": c}
#Calcule la fonction de transfert Sallen-Key pour un filtre passe-haut Bessel. 
    def sallen_key_highpass(self, order, cutoff_freq, r1=None, r2=None, c1=None, c2=None, omega0_norm=None, q0=None):

        omega0_norm_hp = 1 / omega0_norm
        print(omega0_norm_hp)
        omega0 = 2 * np.pi * cutoff_freq * omega0_norm_hp

        if r1 is not None and r2 is not None:
            c1_plus_c2 = 1 / (omega0 * r1 * q0)
            c1_times_c2 = 1 / (omega0**2 * r1 * r2)

            a = 1
            b = -c1_plus_c2
            c = c1_times_c2
            discriminant = b**2 - 4 * a * c

            if discriminant < 0:
                raise ValueError("Les paramètres fournis pour R1 et R2 ne permettent pas un calcul valide de C1 et C2.")

            c21 = (-b + np.sqrt(discriminant)) / (2 * a)
            c22 = (-b - np.sqrt(discriminant)) / (2 * a)
            c11 = c1_plus_c2 - c21
            c12 = c1_plus_c2 - c22

            num1 = [r1 * r2 * c11 * c21,0,0]
            num2 = [r1 * r2 * c12 * c22,0,0]
            den1 = [r1 * r2 * c11 * c21,(r1 * c11 + r1 * c21),1.0]
            den2 = [r1 * r2 * c12 * c22,(r1 * c12 + r1 * c22),1.0]

            tf1 = TransferFunction(num1, den1)
            tf2 = TransferFunction(num2, den2)

            return [
                {"tf": tf1, "params": {"R1": r1, "R2": r2, "C1": c11, "C2": c21}},
                {"tf": tf2, "params": {"R1": r1, "R2": r2, "C1": c12, "C2": c22}},
            ]

        elif c1 is not None and c2 is not None:
            r1 = 1 / (omega0 * q0 * (c1 + c2))
            r2 = 1 / (omega0**2 * r1 * c1 * c2)
            print(q0)
            print(omega0)

            num = [r1 * r2 * c1 * c2,0,0]
            den = [r1 * r2 * c1 * c2,(r1 * c1 + r1 * c2),1.0]

            tf = TransferFunction(num, den)
            return [
                {"tf": tf, "params": {"R1": r1, "R2": r2, "C1": c1, "C2": c2}},
            ]

        else:
            raise ValueError("Veuillez fournir soit (R1, R2), soit (C1, C2).")
        
# Calcule un filtre passe-haut de n'importe quel ordre en utilisant des cellules en cascade.
    def components(self, order, cutoff_freq, r_vals=None, c_vals=None):
       
        if order not in self.BESSEL_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")

        poles = self.bessel_q0_omega0(order)
        stages = []
        num_combined, den_combined = [1], [1]

        num_stages = order // 2 + (order % 2)
        num_elements = 2 * num_stages

        if r_vals is None:
            r_vals = [None] * num_elements
        if c_vals is None:
            c_vals = [None] * num_elements

        if len(r_vals) < num_elements:
            r_vals.extend([None] * (num_elements - len(r_vals)))
        if len(c_vals) < num_elements:
            c_vals.extend([None] * (num_elements - len(c_vals)))

        for i, (omega0_norm, q0) in enumerate(poles):
            if q0 == 0.0:
                tf, params = self.first_order_highpass(
                    cutoff_freq, r=r_vals[i], c=c_vals[i], omega0_norm=omega0_norm
                )
            else:
                idx = 2 * i
                tf_data = self.sallen_key_highpass(
                    order=2,
                    cutoff_freq=cutoff_freq,
                    r1=r_vals[idx],
                    r2=r_vals[idx + 1],
                    c1=c_vals[idx],
                    c2=c_vals[idx + 1],
                    omega0_norm=omega0_norm,
                    q0=q0,
                )
                tf = tf_data[0]["tf"]
                params = tf_data[0]["params"]
                print(tf)
               

            num_combined = np.polymul(num_combined, tf.num)
            print(num_combined)
            den_combined = np.polymul(den_combined, tf.den)
            print(den_combined)
            stages.append({"tf": tf, "params": params})

        combined_tf = TransferFunction(num_combined, den_combined)
        print(combined_tf)
        return combined_tf, stages
    def graphs(self, order, cutoff_freq, r_vals=None, c_vals=None):
        combined_tf, _ = self.components(order, cutoff_freq, r_vals, c_vals)
        w2 = np.logspace(2, 6, 500)
        w2, mag2, phase2 = bode(combined_tf, w=w2)
        freq_hz2 = w2/(2*np.pi)

        fig_hp,(ax_mag_hp, ax_phase_hp) = plt.subplots(2,1,figsize=(8,6), sharex=True)
        ax_mag_hp.semilogx(freq_hz2, mag2, 'b')
        ax_mag_hp.set_ylabel('Magnitude (dB)')

        ax_phase_hp.semilogx(freq_hz2, phase2, 'r')
        ax_phase_hp.set_xlabel('Frequency (Hz)')
        ax_phase_hp.set_ylabel('Phase (deg)')
        ax_phase_hp.grid(True)
        ax_mag_hp.grid(True)
        plt.tight_layout()
        plt.show()