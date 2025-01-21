import numpy as np
from scipy.signal import TransferFunction, bode
import matplotlib.pyplot as plt

class TchebychevFilter:
    """
    Classe unique pour construire un filtre Tchebychev passe-bas OU passe-haut
    d'ordre n, avec ou sans cellule d'ordre 1 (si un pôle a q=0).

    - On utilise des formules directes pour la cellule ordre 2.
    - On stocke une table de pôles (omega0_norm, q0) adaptée à un certain ripple.
    """

    def __init__(self):
        # Table d'exemple (à adapter selon votre ripple).
        # On suppose qu'on a la version "passe-bas" de la table.
        self.TCHEBYCHEV_TABLE = {
            1: [(1.9652, 0.0)],
            2: [(1.0500, 0.9565)],
            3: [(1.4942, 0.0), (0.9971, 2.0177)],
            4: [(0.5286, 0.7845), (0.9932, 3.5590)],
            5: [(0.2895, 0.0), (0.6552, 1.3988), (0.9941, 5.5564)],
            # etc.
        }

    def tchebychev_poles(self, order):
        """Retourne la liste (omega0_norm, q0) pour l'ordre donné."""
        if order not in self.TCHEBYCHEV_TABLE:
            raise ValueError(f"Table indisponible pour l'ordre {order}.")
        return self.TCHEBYCHEV_TABLE[order]

    # ----------------------------------------------------------------
    # 1) Cellule 1er ordre
    # ----------------------------------------------------------------
    def first_order_filter(self, filter_type, cutoff_freq, R=None, C=None, omega0_norm=1.0):
        """
        Construit la FT d'une cellule 1er ordre :
          - Passe-bas : H(s) = 1 / [1 + s R C]
          - Passe-haut : H(s) = (s R C) / [1 + s R C]

        On définit la pulsation selon le type :
          - w_LP = 2*pi*fc * omega0_norm
          - w_HP = (2*pi*fc) / omega0_norm

        R*C = 1 / w  (selon LP ou HP).
        """
        if filter_type not in ["lowpass", "highpass"]:
            raise ValueError("filter_type doit être 'lowpass' ou 'highpass'.")

        if filter_type == "lowpass":
            omega = 2*np.pi*cutoff_freq * omega0_norm
            # R*C = 1/omega
        else:
            # highpass
            omega = (2*np.pi*cutoff_freq) / omega0_norm
            # R*C = 1/omega

        if R is not None and C is None:
            C = 1/(omega * R)
        elif C is not None and R is None:
            R = 1/(omega * C)
        elif R is None and C is None:
            raise ValueError("Fournir R ou C pour la cellule 1er ordre.")

        if filter_type == "lowpass":
            # H(s) = 1 / [1 + sRC]
            num = [1.0]
            den = [R*C, 1.0]
        else:
            # highpass => H(s) = sRC / [1 + sRC]
            num = [R*C, 0]
            den = [R*C, 1.0]

        tf = TransferFunction(num, den)
        return tf, {"R": R, "C": C}

    # ----------------------------------------------------------------
    # 2) Cellule 2e ordre "directe"
    # ----------------------------------------------------------------
    def second_order_filter_direct(self, filter_type, cutoff_freq, C1, C2, omega0_norm=1.0, q0=1.0):
        """
        Construit une cellule 2e ordre, passe-bas ou passe-haut,
        via les formules directes:

        - Passe-bas :
            R1 (C1 + C2) = 1 / [Q * w_LP],
            R1 R2 C1 C2  = 1 / w_LP^2,
            H(s) = 1 / [1 + s(R1)(C1+C2) + s^2(R1R2C1C2)].

        - Passe-haut :
            R1 (C1 + C2) = 1 / [Q * w_HP],
            R1 R2 C1 C2  = 1 / w_HP^2,
            H(s) = (s^2 R1R2C1C2)/[ s^2 R1R2C1C2 + s R1(C1+C2) + 1].
        """
        if filter_type not in ["lowpass", "highpass"]:
            raise ValueError("filter_type doit être 'lowpass' ou 'highpass'.")

        if filter_type == "lowpass":
            omega = 2*np.pi*cutoff_freq * omega0_norm
        else:
            omega = (2*np.pi*cutoff_freq) / omega0_norm

        # 1) Calcul R1
        denom_R1 = q0 * omega * (C1 + C2)
        if denom_R1 <= 0:
            raise ValueError("Impossible de calculer R1 (dénominateur <= 0).")
        R1 = 1/denom_R1

        # 2) Calcul R2
        denom_R2 = (omega**2)*C1*C2*R1
        if denom_R2 <= 0:
            raise ValueError("Impossible de calculer R2 (dénominateur <= 0).")
        R2 = 1/denom_R2

        if filter_type == "lowpass":
            # num = 1
            # den = [R1R2C1C2, R1(C1+C2), 1]
            num = [1.0]
            den = [R1*R2*C1*C2, R1*(C1 + C2), 1.0]
        else:
            # HP => num = s^2 R1R2C1C2
            # den = s^2 R1R2C1C2 + s R1(C1+C2) + 1
            num = [R1*R2*C1*C2, 0, 0]
            den = [R1*R2*C1*C2, R1*(C1 + C2), 1.0]

        tf = TransferFunction(num, den)
        return tf, {"R1": R1, "R2": R2, "C1": C1, "C2": C2}

    # ----------------------------------------------------------------
    # 3) Conception d'un filtre d'ordre n
    # ----------------------------------------------------------------
    def design_filter(self, order, cutoff_freq, filter_type="lowpass", c_vals=None, r_vals=None):
        """
        - order : ordre du filtre
        - cutoff_freq : fréquence de coupure (Hz)
        - filter_type : 'lowpass' ou 'highpass'
        - c_vals / r_vals : listes de longueur = order (facultatives)
                            pour imposer des composants.

        On cascade chaque pôle (omega0_norm, q0) :
         - si q0=0 => cellule 1er ordre
         - sinon => cellule 2e ordre
        """
        if filter_type not in ["lowpass", "highpass"]:
            raise ValueError("filter_type doit être 'lowpass' ou 'highpass'.")

        poles = self.tchebychev_poles(order)

        # Vérif longueur
        if c_vals is not None and len(c_vals) != order:
            raise ValueError(f"c_vals doit avoir {order} éléments.")
        if r_vals is not None and len(r_vals) != order:
            raise ValueError(f"r_vals doit avoir {order} éléments.")

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
                # => 1er ordre
                tf1, params1 = self.first_order_filter(
                    filter_type,
                    cutoff_freq,
                    R=r_vals[idx],
                    C=c_vals[idx],
                    omega0_norm=omega0_norm
                )
                idx += 1
                stages.append({"tf": tf1, "params": params1})
                num_combined = np.polymul(num_combined, tf1.num)
                den_combined = np.polymul(den_combined, tf1.den)

            else:
                # => 2e ordre
                # On va consommer 2 valeurs (C1, C2) et pas de R imposé => calcul direct
                C1, C2 = c_vals[idx], c_vals[idx+1]
                idx += 2
                if (C1 is None) or (C2 is None):
                    raise ValueError("Cellule 2e ordre => il faut C1 et C2 (ou on revoit la logique).")

                tf2, params2 = self.second_order_filter_direct(
                    filter_type,
                    cutoff_freq,
                    C1, C2,
                    omega0_norm=omega0_norm,
                    q0=q0
                )
                stages.append({"tf": tf2, "params": params2})
                num_combined = np.polymul(num_combined, tf2.num)
                den_combined = np.polymul(den_combined, tf2.den)

        # TF globale
        tf_global = TransferFunction(num_combined, den_combined)
        return tf_global, stages


# ----------------------------------------------------------------
# Exemple d'utilisation
# ----------------------------------------------------------------
if __name__ == "__main__":
    # On crée notre objet
    filter_designer = TchebychevFilter()

    # 1) Ex: un filtre passe-bas d'ordre 3, fc=1 kHz
    order_lp = 3
    fc_lp = 1e3
    # On fixe 3 condensateurs, p.ex. c_vals=[10nF, 10nF, 4.7nF],
    # sachant que la 1ʳᵉ cellule (q=0) n'utilisera qu'un seul condo,
    # et la 2ᵉ (q!=0) en utilisera 2.
    cvals_lp = [10e-9, 10e-9, 4.7e-9]

    tf_lp, stages_lp = filter_designer.design_filter(
        order=order_lp,
        cutoff_freq=fc_lp,
        filter_type="lowpass",
        c_vals=cvals_lp,
        r_vals=None  # on laisse calculer R
    )
    print("\n=== PASSE-BAS Ordre 3 ===")
    for i, st in enumerate(stages_lp, start=1):
        print(f"Cellule {i} => {st['params']}")
    print("TF globale (num, den):", tf_lp.num, tf_lp.den)

    # Bode plot LP
    w = np.logspace(2, 5, 400)
    w, mag, phase = bode(tf_lp, w=w)
    freq_hz = w/(2*np.pi)

    fig_lp, (ax_mag_lp, ax_phase_lp) = plt.subplots(2,1,figsize=(8,6), sharex=True)
    ax_mag_lp.semilogx(freq_hz, mag, 'b')
    ax_mag_lp.set_ylabel('Magnitude (dB)')
    ax_mag_lp.set_title(f'Tchebychev LP ordre={order_lp}, Fc={fc_lp/1e3} kHz')

    ax_phase_lp.semilogx(freq_hz, phase, 'r')
    ax_phase_lp.set_xlabel('Frequency (Hz)')
    ax_phase_lp.set_ylabel('Phase (deg)')
    ax_phase_lp.grid(True)
    ax_mag_lp.grid(True)
    plt.tight_layout()
    plt.show()


    # 2) Ex: un filtre passe-haut d'ordre 4, fc=2.5 kHz
    order_hp = 4
    fc_hp = 2.5e3
    # 4 condensateurs => 2 cellules d'ordre 2, par ex. [10nF, 10nF, 10nF, 5nF]
    cvals_hp = [10e-9, 10e-9, 10e-9, 5e-9]

    tf_hp, stages_hp = filter_designer.design_filter(
        order=order_hp,
        cutoff_freq=fc_hp,
        filter_type="highpass",
        c_vals=cvals_hp,
        r_vals=None
    )
    print("\n=== PASSE-HAUT Ordre 4 ===")
    for i, st in enumerate(stages_hp, start=1):
        print(f"Cellule {i} => {st['params']}")
    print("TF globale (num, den):", tf_hp.num, tf_hp.den)

    # Bode plot HP
    w2 = np.logspace(2, 6, 500)
    w2, mag2, phase2 = bode(tf_hp, w=w2)
    freq_hz2 = w2/(2*np.pi)

    fig_hp, (ax_mag_hp, ax_phase_hp) = plt.subplots(2,1,figsize=(8,6), sharex=True)
    ax_mag_hp.semilogx(freq_hz2, mag2, 'b')
    ax_mag_hp.set_ylabel('Magnitude (dB)')
    ax_mag_hp.set_title(f'Tchebychev HP ordre={order_hp}, Fc={fc_hp/1e3} kHz')

    ax_phase_hp.semilogx(freq_hz2, phase2, 'r')
    ax_phase_hp.set_xlabel('Frequency (Hz)')
    ax_phase_hp.set_ylabel('Phase (deg)')
    ax_phase_hp.grid(True)
    ax_mag_hp.grid(True)
    plt.tight_layout()
    plt.show()
