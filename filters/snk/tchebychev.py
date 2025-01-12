import numpy as np
from scipy.signal import TransferFunction, bode
import matplotlib.pyplot as plt

class Highpasstchebychev:
    """
    Librairie pour calculer un filtre Tchebychev passe-haut d'ordre n
    en utilisant la même table de pôles (omega0_norm, q0) que pour le passe-bas,
    MAIS en appliquant la relation directe :
        R1 (C1 + C2) = 1 / (Q * omega_HP)
        R1 * R2 * C1 * C2 = 1 / (omega_HP^2)
    pour la cellule d'ordre 2 (pas de discriminant).
    """

    def __init__(self):
        # Table de pôles Tchebychev : (omega0_norm, q0) pour la version passe-bas
        # On inversera la fréquence pour le passe-haut (omega_HP = (2*pi*fc)/omega0_norm).
        self.TCHEBYCHEV_TABLE = {
            1: [(1.9652, 0.0)],
            2: [(1.0500, 0.9565)],
            3: [(1.4942, 0.0), (0.9971, 2.0177)],
            4: [(0.5286, 0.7845), (0.9932, 3.5590)],
            5: [(0.2895, 0.0), (0.6552, 1.3988), (0.9941, 5.5564)],
            6: [(0.3531, 0.7609), (0.7468, 2.1980), (0.9954, 8.0037)],
            7: [(0.2054, 0.0), (0.4801, 1.2969), (0.8084, 3.1559), (0.9963, 10.8987)],
            8: [(0.2651, 0.7530), (0.5828, 1.9565), (0.8506, 4.2661), (0.9971, 14.2405)],
            9: [(0.1593, 0.0), (0.3773, 1.2600), (0.6622, 2.7129), (0.8806, 5.5266), (0.9976, 18.0286)],
           10: [(0.2121, 0.7495), (0.4761, 1.8645), (0.7215, 3.5605), (0.9025, 6.9367), (0.9980, 22.2630)],
        }

    def tchebychev_q0_omega0(self, order):
        """Retourne la liste des pôles (omega0_norm, q0) pour un ordre donné."""
        if order not in self.TCHEBYCHEV_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")
        return self.TCHEBYCHEV_TABLE[order]

    # -----------------------------------------------------------
    # 1) Cellule 1er ordre passe-haut
    # -----------------------------------------------------------
    def first_order_highpass(self, cutoff_freq, R=None, C=None, omega0_norm=1.0):
        """
        Construit un filtre passe-haut de 1er ordre :
           H_HP(s) = (s R C) / (1 + s R C).
        On définit : omega_HP = (2*pi*cutoff_freq) / omega0_norm.
        Puis, R*C = 1/omega_HP.
        """
        omega_hp = (2 * np.pi * cutoff_freq) / omega0_norm

        # Calcul direct
        if R is not None and C is None:
            C = 1 / (omega_hp * R)
        elif C is not None and R is None:
            R = 1 / (omega_hp * C)
        elif R is None and C is None:
            raise ValueError("Fournir R ou C (1er ordre HP).")

        num = [R*C, 0]
        den = [R*C, 1]
        tf = TransferFunction(num, den)
        return tf, {"R": R, "C": C}

    # -----------------------------------------------------------
    # 2) Cellule 2e ordre passe-haut (approche directe)
    # -----------------------------------------------------------
    def sallen_key_highpass_direct(self, cutoff_freq, C1, C2, omega0_norm=1.0, q0=1.0):
        """
        Construit un filtre passe-haut de 2e ordre via Sallen–Key,
        en utilisant les relations directes (sans discriminant) :

           1)  R1 (C1 + C2) = 1 / (Q * omega_HP)
           2)  R1 R2 C1 C2  = 1 / (omega_HP^2)

        => R1 = 1 / [ Q * omega_HP * (C1 + C2) ]
           R2 = 1 / [ omega_HP^2 * C1 * C2 * R1 ]

        La TF obtenue est :
            H(s) = (s^2 R1 R2 C1 C2) / [ s^2 R1 R2 C1 C2 + s R1 (C1 + C2) + 1 ].
        """
        omega_hp = (2 * np.pi * cutoff_freq) / omega0_norm

        # 1) Calcul de R1
        denom_R1 = q0 * omega_hp * (C1 + C2)
        if denom_R1 <= 0:
            raise ValueError("Impossible de calculer R1 (dénominateur <= 0).")
        R1 = 1 / denom_R1

        # 2) Calcul de R2
        denom_R2 = (omega_hp**2) * C1 * C2 * R1
        if denom_R2 <= 0:
            raise ValueError("Impossible de calculer R2 (dénominateur <= 0).")
        R2 = 1 / denom_R2

        # 3) Construction de la FT
        num = [R1*R2*C1*C2, 0, 0]
        den = [R1*R2*C1*C2, R1*(C1 + C2), 1]
        tf = TransferFunction(num, den)

        return tf, {"R1": R1, "R2": R2, "C1": C1, "C2": C2}

    # -----------------------------------------------------------
    # 3) Conception d'un filtre d'ordre n (cascade)
    # -----------------------------------------------------------
    def design_filter(self, order, cutoff_freq, c_vals=None, r_vals=None):
        """
        Conception d'un filtre Tchebychev passe-haut d'ordre 'order'.
        - On cascade (n//2) cellules d'ordre 2 et, si n est impair, 1 cellule d'ordre 1.
        - On récupère pour chaque pôle (omega0_norm, q0).
          * Si q0=0 => cellule d'ordre 1 (HP).
          * Si q0!=0 => cellule d'ordre 2.

        On impose ci-dessous que, si c_vals est fourni, on l'utilise
        pour chaque cellule d'ordre 2, et on calcule R1, R2 par la formule directe.
        (Ou l'inverse pour la 1re ordre.)

        c_vals et r_vals doivent avoir la longueur = order
        si on veut imposer explicitement les compos (pas obligatoire).
        """
        poles = self.tchebychev_q0_omega0(order)

        # Vérif longueurs
        if c_vals is not None and len(c_vals) != order:
            raise ValueError(f"c_vals doit avoir exactement {order} éléments.")
        if r_vals is not None and len(r_vals) != order:
            raise ValueError(f"r_vals doit avoir exactement {order} éléments.")

        if c_vals is None:
            c_vals = [None]*order
        if r_vals is None:
            r_vals = [None]*order

        num_combined = [1.]
        den_combined = [1.]
        stages = []
        idx = 0

        for (omega0_norm, q0) in poles:
            if q0 == 0.0:
                # ---- 1er ordre HP ----
                tf1, params1 = self.first_order_highpass(
                    cutoff_freq,
                    R=r_vals[idx],   # si imposé
                    C=c_vals[idx],   # si imposé
                    omega0_norm=omega0_norm
                )
                idx += 1
                stages.append({"tf": tf1, "params": params1})
                num_combined = np.polymul(num_combined, tf1.num)
                den_combined = np.polymul(den_combined, tf1.den)

            else:
                # ---- 2e ordre HP ----
                # On va APPLIQUER LA FORMULE DIRECTE:
                #    R1 (C1 + C2) = 1/(Q * omega_HP)
                #    R1 R2 C1 C2 = 1/(omega_HP^2)
                #
                # Suppose qu'on fixe c_vals[idx], c_vals[idx+1]
                # et qu'on calcule R1, R2 par la sallen_key_highpass_direct.
                #
                # c_vals[idx], c_vals[idx+1] -> C1, C2
                c1, c2 = c_vals[idx], c_vals[idx+1]
                idx += 2

                if (c1 is None) or (c2 is None):
                    raise ValueError("Pour la cellule d'ordre 2, il faut C1 et C2 (sinon impossible).")

                # Appel direct
                tf2, params2 = self.sallen_key_highpass_direct(
                    cutoff_freq,
                    c1, c2,
                    omega0_norm=omega0_norm,
                    q0=q0
                )
                stages.append({"tf": tf2, "params": params2})
                num_combined = np.polymul(num_combined, tf2.num)
                den_combined = np.polymul(den_combined, tf2.den)

        combined_tf = TransferFunction(num_combined, den_combined)
        return combined_tf, stages


# -----------------------------------------------------------------------
# Exemple d’utilisation
if __name__ == "__main__":
    # Filtre Tchebychev passe-haut d'ordre 4, Fc = 2.5 kHz
    order = 4
    fc = 2.5e3

    # On fixe 4 condensateurs => (1er ordre : 1 condo, 2ᵉ ordre : 2 condos, etc.)
    # ex. 10nF, 10nF, 10nF, 5nF
    # => la 1ʳᵉ cellule a q=0 => 1er ordre => c_vals[0] = 10nF
    # => la 2ʳᵉ cellule a q!=0 => c_vals[1], c_vals[2], c_vals[3], etc.
    #
    # Mais pour un ordre 4 Tchebychev, on a 2 pôles => chacun q != 0
    # => En réalité, selon votre table, vérifiez s'il y a un pôle q=0.
    # => S'il n'y en a pas, on aura DEUX cellules d'ordre 2 => on a besoin de 4 condos total.
    #
    # Pour simplifier, on indique 4 valeurs :
    c_vals = [10e-9, 10e-9, 10e-9, 5e-9]

    # On ne spécifie pas R => le code va calculer R1, R2 par la relation directe.
    r_vals = None

    hp_filter = Highpasstchebychev()
    tf_hp, stages_info = hp_filter.design_filter(
        order=order,
        cutoff_freq=fc,
        c_vals=c_vals,
        r_vals=r_vals
    )

    print("=== STAGES INFO ===")
    for i, st in enumerate(stages_info, start=1):
        print(f"Cellule {i}: {st['params']}")

    print("\n=== FONCTION DE TRANSFERT GLOBALE ===")
    print("Num =", tf_hp.num)
    print("Den =", tf_hp.den)

    # Bode plot
    w = np.logspace(2, 6, 500)
    w, mag, phase = bode(tf_hp, w=w)
    freq_hz = w / (2*np.pi)

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax_mag.semilogx(freq_hz, mag, 'b')
    ax_mag.set_ylabel('Magnitude (dB)')
    ax_mag.set_title(f'Tchebychev HP ordre={order} - Fc={fc/1e3} kHz (Relation directe)')

    ax_phase.semilogx(freq_hz, phase, 'r')
    ax_phase.set_xlabel('Fréquence (Hz)')
    ax_phase.set_ylabel('Phase (deg)')
    ax_phase.grid(True)
    ax_mag.grid(True)

    plt.tight_layout()
    plt.show()
