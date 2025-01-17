import pytest
import numpy as np
from scipy.signal import bode
from filters.snk.tchebychev import TchebychevFilter

def test_lowpass_order3():
    """Test pour le filtre passe-bas d'ordre 3."""
    filter_designer = TchebychevFilter()
    order_lp = 3
    fc_lp = 1e3
    cvals_lp = [10e-9, 10e-9, 4.7e-9]

    tf_lp, stages_lp = filter_designer.design_filter(
        order=order_lp,
        cutoff_freq=fc_lp,
        filter_type="lowpass",
        c_vals=cvals_lp,
        r_vals=None
    )

    print(f"\n=== Test Passe-Bas Ordre {order_lp} ===")
    print(f"Fréquence de coupure : {fc_lp} Hz")
    print(f"Nombre de cellules obtenues : {len(stages_lp)}")
    for i, stage in enumerate(stages_lp, start=1):
        print(f"Cellule {i} paramètres : {stage['params']}")
        print(f"Cellule {i} fonction de transfert : {stage['tf']}")

    assert len(stages_lp) == 2, "Nombre incorrect de cellules dans le filtre passe-bas."
    for stage in stages_lp:
        if "R" in stage["params"] and "C" in stage["params"]:
            assert "R" in stage["params"] and "C" in stage["params"], "Paramètres R et C manquants dans une cellule d'ordre 1."
        elif "R1" in stage["params"] and "R2" in stage["params"]:
            assert "R1" in stage["params"] and "R2" in stage["params"], "Paramètres R1 et R2 manquants dans une cellule d'ordre 2."

def test_highpass_order4():
    """Test pour le filtre passe-haut d'ordre 4."""
    filter_designer = TchebychevFilter()
    order_hp = 4
    fc_hp = 2.5e3
    cvals_hp = [10e-9, 10e-9, 10e-9, 5e-9]

    tf_hp, stages_hp = filter_designer.design_filter(
        order=order_hp,
        cutoff_freq=fc_hp,
        filter_type="highpass",
        c_vals=cvals_hp,
        r_vals=None
    )

    print(f"\n=== Test Passe-Haut Ordre {order_hp} ===")
    print(f"Fréquence de coupure : {fc_hp} Hz")
    print(f"Nombre de cellules obtenues : {len(stages_hp)}")
    for i, stage in enumerate(stages_hp, start=1):
        print(f"Cellule {i} paramètres : {stage['params']}")
        print(f"Cellule {i} fonction de transfert : {stage['tf']}")

    assert len(stages_hp) == 2, "Nombre incorrect de cellules dans le filtre passe-haut."
    for stage in stages_hp:
        assert "R1" in stage["params"] and "R2" in stage["params"], "Paramètres R1 et R2 manquants dans une cellule."

def test_bode_response():
    """Test pour vérifier la réponse en fréquence."""
    filter_designer = TchebychevFilter()
    order_lp = 3
    fc_lp = 1e3
    cvals_lp = [10e-9, 10e-9, 4.7e-9]

    tf_lp, _ = filter_designer.design_filter(
        order=order_lp,
        cutoff_freq=fc_lp,
        filter_type="lowpass",
        c_vals=cvals_lp,
        r_vals=None
    )

    w = np.logspace(2, 5, 400)
    w, mag, phase = bode(tf_lp, w=w)

    print(f"\n=== Test Réponse en Fréquence Passe-Bas Ordre {order_lp} ===")
    print(f"Fréquence de coupure : {fc_lp} Hz")
    print(f"Longueur du vecteur fréquence : {len(w)}")
    print(f"Premières valeurs de magnitude (dB) : {mag[:5]}")
    print(f"Premières valeurs de phase (degrés) : {phase[:5]}")

    assert len(mag) == len(w), "La longueur de la réponse magnitude ne correspond pas."
    assert len(phase) == len(w), "La longueur de la réponse phase ne correspond pas."
