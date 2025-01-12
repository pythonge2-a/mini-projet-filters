import numpy as np
from scipy.signal import TransferFunction
import matplotlib.pyplot as plt

class BesselFilter:
    """
    Classe pour construire un filtre Bessel (passe-bas ou passe-haut) d'ordre n,
    avec table de pôles (omega0_norm, q0).

    Logique similaire à l'exemple unifié :
      - first_order_filter  (passe-bas / passe-haut)
      - second_order_filter (passe-bas / passe-haut) en formules directes
      - design_filter       (cascade)
      - bode_plot           (affichage)
    """

    def __init__(self):
        # Table Bessel : (omega0_norm, q0)
        self.BESSEL_TABLE = {
            1: [(1.0, 0.0)],
            2: [(1.2723, 0.577)],
            3: [(1.3225, 0.0), (1.4474, 0.691)],
            4: [(1.431, 0.5219), (1.6043, 0.8055)],
            5: [(1.5015, 0.0), (1.5555, 0.5635), (1.7545, 0.9165)],
            6: [(1.6030, 0.5103), (1.6882, 0.6112), (1.9037, 1.0233)],
            7: [(1.6840, 0.0), (1.7160, 0.5324), (1.8221, 0.6608), (2.0491, 1.0233)],
            8: [(1.7772, 0.5060), (1.8308, 0.5596), (1.9518, 0.7109), (2.1872, 1.2257)],
            9: [(1.8570, 0.0), (1.8788, 0.5197), (1.9483, 0.5895), (2.0808, 0.7606), (2.3228, 1.3219)],
           10: [(1.9412, 0.5039), (1.9790, 0.5376), (2.0606, 0.6205), (2.2021, 0.8098), (2.4487, 1.4153)],
        }

    def bessel_poles(self, order):
        """Retourne la liste (omega0_norm, q0) pour un ordre donné."""
        if order not in self.BESSEL_TABLE:
            raise ValueError(f"L'ordre {order} n'est pas supporté.")
        return self.BESSEL_TABLE[order]

    # ----------------------------------------------------------------
    # 1) Cellule 1er ordre
    # ----------------------------------------------------------------
    def first_order_filter(self, filter_type, cutoff_freq, R=None, C=None, omega0_norm=1.0):
        """
        Construit la FT d'une cellule 1er ordre (Bessel) en passe-bas ou passe-haut.
          - Passe-bas : H(s) = 1 / [1 + s R C]
          - Passe-haut : H(s) = s R C / [1 + s R C]

        On définit :
          w_LP = 2*pi*fc * omega0_norm    (passe-bas)
          w_HP = (2*pi*fc) / omega0_norm  (passe-haut)

        R*C = 1/w_{LP} ou 1/w_{HP} selon le type.
        """
        if filter_type not in ["lowpass", "highpass"]:
            raise ValueError("filter_type doit être 'lowpass' ou 'highpass'.")

        if filter_type == "lowpass":
            w_lp = 2*np.pi*cutoff_freq * omega0_norm
            # R*C = 1 / w_lp
            if R is not None and C is None:
                C = 1/(w_lp*R)
            elif C is not None and R is None:
                R = 1/(w_lp*C)
            else:
                if R is None and C is None:
                    raise ValueError("Fournir R ou C (1er ordre LP).")
            # H(s) = 1 / (1 + sRC)
            num = [1.0]
            den = [R*C, 1.0]

        else:
            # highpass
            w_hp = (2*np.pi*cutoff_freq)/omega0_norm
            # R*C = 1 / w_hp
            if R is not None and C is None:
                C = 1/(w_hp*R)
            elif C is not None and R is None:
                R = 1/(w_hp*C)
            else:
                if R is None and C is None:
                    raise ValueError("Fournir R ou C (1er ordre HP).")
            # H(s) = (sRC)/(1 + sRC)
            num = [R*C, 0]
            den = [R*C, 1.0]

        tf = TransferFunction(num, den)
        return tf, {"R": R, "C": C}

    # ----------------------------------------------------------------
    # 2) Cellule 2e ordre (formules directes)
    # ----------------------------------------------------------------
    def second_order_filter_bessel(self, filter_type, cutoff_freq, C1, C2, omega0_norm=1.0, q0=1.0):
        """
        Cellule Bessel d'ordre 2 (passe-bas ou passe-haut) via formules directes :

        Passe-bas :
          w_LP = 2*pi*fc * omega0_norm
          R1 (C1 + C2) = 1/(Q * w_LP)
          R1 R2 C1 C2  = 1/(w_LP^2)
          => FT(s) = 1 / [1 + s R1 (C1 + C2) + s^2 R1 R2 C1 C2]

        Passe-haut :
          w_HP = (2*pi*fc)/omega0_norm
          R1 (C1 + C2) = 1/(Q * w_HP)
          R1 R2 C1 C2  = 1/(w_HP^2)
          => FT(s) = s^2 R1 R2 C1 C2 / [ s^2 R1 R2 C1 C2 + s R1 (C1 + C2) + 1 ]
        """
        if filter_type not in ["lowpass", "highpass"]:
            raise ValueError("filter_type doit être 'lowpass' ou 'highpass'.")

        if filter_type == "lowpass":
            w_lp = 2*np.pi*cutoff_freq * omega0_norm
            if w_lp <= 0:
                raise ValueError("w_lp <= 0, check fc or omega0_norm.")
            # R1 = 1/(Q * w_lp * (C1 + C2))
            # R2 = 1/(w_lp^2 * C1 * C2 * R1)

            denom_R1 = q0 * w_lp * (C1 + C2)
            if denom_R1 <= 0:
                raise ValueError("Impossible calcul R1 (denom <=0).")
            R1 = 1/denom_R1

            denom_R2 = (w_lp**2)*C1*C2*R1
            if denom_R2 <= 0:
                raise ValueError("Impossible calcul R2 (denom <=0).")
            R2 = 1/denom_R2

            num = [1.0]
            den = [R1*R2*C1*C2, R1*(C1 + C2), 1.0]

        else:
            # highpass
            w_hp = (2*np.pi*cutoff_freq)/omega0_norm
            if w_hp <= 0:
                raise ValueError("w_hp <= 0, check fc or omega0_norm.")

            denom_R1 = q0 * w_hp * (C1 + C2)
            if denom_R1 <= 0:
                raise ValueError("Impossible calcul R1 (denom <=0).")
            R1 = 1/denom_R1

            denom_R2 = (w_hp**2)*C1*C2*R1
            if denom_R2 <= 0:
                raise ValueError("Impossible calcul R2 (denom <=0).")
            R2 = 1/denom_R2

            num = [R1*R2*C1*C2, 0, 0]
            den = [R1*R2*C1*C2, R1*(C1 + C2), 1.0]

        tf = TransferFunction(num, den)
        return tf, {"R1": R1, "R2": R2, "C1": C1, "C2": C2}

    # ----------------------------------------------------------------
    # 3) Conception d'un filtre Bessel d'ordre n
    # ----------------------------------------------------------------
    def design_filter(self, order, cutoff_freq, filter_type="lowpass", c_vals=None, r_vals=None):
        """
        Conception d'un filtre Bessel (passe-bas ou passe-haut) d'ordre 'order'.
        - On cascade (n//2) cellules d'ordre 2, plus éventuellement 1 cellule d'ordre 1 si n est impair.
        - On récupère (omega0_norm, q0) : si q0=0 => 1er ordre, sinon 2ᵉ ordre.
        - c_vals et r_vals (listes) : si vous voulez imposer explicitement 1 condo (ordre 1) ou 2 condos (ordre 2), etc.

        Retourne (tf_global, [liste de stages]).
        """
        if filter_type not in ["lowpass", "highpass"]:
            raise ValueError("filter_type doit être 'lowpass' ou 'highpass'.")

        poles = self.bessel_poles(order)

        # Vérification du nombre d'éléments
        # (Au maximum, pour un ordre n, on aura n composants si c_vals ou r_vals imposés.)
        if c_vals is not None and len(c_vals) < order:
            raise ValueError(f"c_vals doit avoir au moins {order} éléments.")
        if r_vals is not None and len(r_vals) < order:
            raise ValueError(f"r_vals doit avoir au moins {order} éléments.")

        # On complète si None
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
                # => cellule 1er ordre
                tf1, params1 = self.first_order_filter(
                    filter_type,
                    cutoff_freq,
                    R=r_vals[idx],
                    C=c_vals[idx],
                    omega0_norm=omega0_norm
                )
                idx += 1
                stages.append({"tf": tf1, "params": params1})
                # multiplication polynômes
                num_combined = np.polymul(num_combined, tf1.num)
                den_combined = np.polymul(den_combined, tf1.den)

            else:
                # => cellule 2ᵉ ordre
                # On prend c_vals[idx], c_vals[idx+1] => C1, C2
                c1, c2 = c_vals[idx], c_vals[idx+1]
                idx += 2

                tf2, params2 = self.second_order_filter_bessel(
                    filter_type,
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

    # ----------------------------------------------------------------
    # 4) Méthode pour tracer le Bode
    # ----------------------------------------------------------------
    def bode_plot(self, order, cutoff_freq, filter_type="lowpass", c_vals=None, r_vals=None):
        """
        Construit le filtre et affiche le Bode (amplitude + phase).
        """
        tf_global, stages = self.design_filter(order, cutoff_freq, filter_type, c_vals, r_vals)
        w, mag, phase = tf_global.bode()  # w en rad/s, mag en dB, phase en deg

        plt.figure(figsize=(10,6))

        plt.subplot(2,1,1)
        plt.semilogx(w, mag, 'b')
        plt.title(f"Bessel {filter_type} ordre={order}, Fc={cutoff_freq} Hz")
        plt.ylabel("Magnitude (dB)")
        plt.grid(True, which="both", ls="--")

        plt.subplot(2,1,2)
        plt.semilogx(w, phase, 'r')
        plt.xlabel("Pulsation (rad/s)")
        plt.ylabel("Phase (deg)")
        plt.grid(True, which="both", ls="--")

        plt.tight_layout()
        plt.show()

        return tf_global, stages

# -----------------------------------------------------------------------
# Exemple d'utilisation
if __name__ == "__main__":
    bf = BesselFilter()

    # 1) Passe-bas d'ordre 3, Fc=1 kHz
    #    On impose 3 condensateurs => ex. [10nF, 10nF, 4.7nF]
    order_lp = 3
    fc_lp = 1e3
    cvals_lp = [10e-9, 10e-9, 4.7e-9]
    bf.bode_plot(order_lp, fc_lp, filter_type="lowpass", c_vals=cvals_lp)

    # 2) Passe-haut d'ordre 4, Fc=2.5 kHz
    #    On impose 4 condensateurs => ex. [10nF, 10nF, 10nF, 5nF]
    order_hp = 4
    fc_hp = 2.5e3
    cvals_hp = [10e-9, 10e-9, 10e-9, 5e-9]
    bf.bode_plot(order_hp, fc_hp, filter_type="highpass", c_vals=cvals_hp)
