import math


class FilterBase:
    """Base pour les filtres Sallen-Key."""

    # Table des paramètres pour Tchebychev 1 dB
    FILTER_PARAMS = {
    "Tchebychev (r=1dB)": {
        1: [(1.9652, None)],  # Une seule cellule d'ordre 1, pas de Qk nécessaire
        2: [(1.0500, 0.9565)],
        3: [(0.4942, None), (0.9971, 2.0177)],  # Une cellule d'ordre 1 et une cellule d'ordre 2
        4: [(0.5286, 0.7845), (0.9932, 3.5590)],
        5: [(0.2895, None), (0.6552, 1.3988), (0.9941, 5.5564)],  # Une cellule d'ordre 1 et deux cellules d'ordre 2
        6: [(0.3531, 0.7609), (0.7468, 2.1980), (0.9954, 8.0037)],
        7: [(0.2054, None), (0.4801, 1.2969), (0.8084, 3.1559), (0.9963, 10.8987)],
        8: [(0.2651, 0.7530), (0.5828, 1.9565), (0.8506, 4.2661), (0.9971, 14.2405)],
        9: [(0.1593, None), (0.3773, 1.2600), (0.6622, 2.7129), (0.8806, 5.5266), (0.9976, 18.0286)],
        10: [(0.2121, 0.7495), (0.4761, 1.8645), (0.7215, 3.5605), (0.9025, 6.9367), (0.9980, 22.2630)],
    },
}

    @staticmethod
    def get_params(filter_type, order, mode):
        """Récupère les paramètres du filtre et ajuste selon le mode."""
        if filter_type not in FilterBase.FILTER_PARAMS:
            raise ValueError(f"Type de filtre non supporté : {filter_type}")
        if order not in FilterBase.FILTER_PARAMS[filter_type]:
            raise ValueError(f"Ordre {order} non supporté pour le filtre {filter_type}")

        params = FilterBase.FILTER_PARAMS[filter_type][order]

        # Ajustement des wk/wr pour passe-haut (1 / wk/wr)
        adjusted_params = []
        for wk_wr, qk in params:
            original_wk_wr = wk_wr
            if mode == "hp" and wk_wr is not None:
                wk_wr = 1 / wk_wr
            adjusted_params.append((wk_wr, qk, original_wk_wr))  # Inclut l'original pour affichage
        return adjusted_params

    @staticmethod
    def calculate_resistances(f, capacitors, params, order):
        """
        Calcule les résistances pour un filtre de Sallen-Key pour un ordre donné.
        Gère les ordres pairs et impairs.
        """
        resistances = {}
        num_cells = len(params)  # Chaque cellule d'ordre 2 a un jeu de (wk/wr, Q)
        capacitor_index = 0  # Index des condensateurs utilisés

        print("\n--- Détails des calculs pour chaque cellule ---")
        for i, (wk_wr, q, original_wk_wr) in enumerate(params):
            if q is not None:  # Cellule d'ordre 2
                try:
                    c1, c2 = capacitors[capacitor_index], capacitors[capacitor_index + 1]
                except IndexError:
                    raise ValueError(
                        f"La liste des condensateurs est trop courte pour l'ordre {order}. "
                        f"Au moins {2 * num_cells} condensateurs sont nécessaires."
                    )
                capacitor_index += 2
                w = wk_wr * f * 2 * math.pi

                r1 = 1 / (q * w * (c1 + c2))
                r2 = 1 / (w**2 * c1 * c2 * r1)

                resistances[f"R{2 * i + 1}"] = r1
                resistances[f"R{2 * i + 2}"] = r2

                print(
                    f"Cellule {i + 1} (ordre 2) : original wk/wr={original_wk_wr:.4f}, utilisé wk/wr={wk_wr:.4f}, "
                    f"Q={q:.4f}, C1={c1:.1e}, C2={c2:.1e}, R1={r1:.1f}, R2={r2:.1f}"
                )
            else:  # Cellule d'ordre 1
                try:
                    c1 = capacitors[capacitor_index]
                except IndexError:
                    raise ValueError(
                        f"La liste des condensateurs est trop courte pour l'ordre {order}. "
                        f"Au moins {num_cells + 1} condensateurs sont nécessaires."
                    )
                capacitor_index += 1
                w = wk_wr * f * 2 * math.pi  # Utilisation correcte de wk_wr pour les cellules d'ordre 1

                r1 = 1 / (w * c1)
                resistances[f"R{2 * i + 1}"] = r1

                print(
                    f"Cellule {i + 1} (ordre 1) : original wk/wr={original_wk_wr:.4f}, utilisé wk/wr={wk_wr:.4f}, "
                    f"C1={c1:.1e}, R1={r1:.1f}"
                )

        print("--- Fin des calculs pour les cellules ---")
        return resistances



class FilterModule:
    """Module pour appeler des fonctions comme `filters.snk.tchebychev.hp.order3`."""

    def __init__(self, filter_type, mode):
        self.filter_type = filter_type
        self.mode = mode

    def order(self, f, capacitors, order):
        """Calcule les résistances pour un filtre d'un ordre donné."""
        params = FilterBase.get_params(self.filter_type, order, self.mode)
        return FilterBase.calculate_resistances(f, capacitors, params, order)


# Initialisation des modules pour tous les types de filtres
class Filters:
    tchebychev = {
        "lp": FilterModule("Tchebychev (r=1dB)", "lp"),
        "hp": FilterModule("Tchebychev (r=1dB)", "hp"),
    }


filters = Filters()

if __name__ == "__main__":
    print("\n--- Test Tchebychev ordre 4 passe-haut ---")
    capacitors = [10e-9, 10e-9, 10e-9, 5e-9]
    params_hp_tchebychev_order4 = FilterBase.get_params("Tchebychev (r=1dB)", 4, "hp")
    print("\n--- Paramètres récupérés depuis le tableau ---")
    for i, (wk_wr, qk, original_wk_wr) in enumerate(params_hp_tchebychev_order4, start=1):
        wk_wr_str = f"{wk_wr:.4f}" if wk_wr is not None else "N/A"
        original_wk_wr_str = f"{original_wk_wr:.4f}" if original_wk_wr is not None else "N/A"
        qk_str = f"{qk:.4f}" if qk is not None else "None"
        print(
            f"  Cellule {i}: original wk/wr = {original_wk_wr_str}, utilisé wk/wr = {wk_wr_str}, Qk = {qk_str}"
        )

    result_hp_tchebychev = filters.tchebychev["hp"].order(f=2500, capacitors=capacitors, order=4)

    print("\n--- Résistances calculées (arrondies au demi-ohm) ---")
    for res_name, res_value in result_hp_tchebychev.items():
        print(f"  {res_name} = {res_value:.1f} Ω")






    print("\n--- Test Tchebychev ordre 3 passe-bas ---")

    # Capacitors utilisés
    capacitors = [10e-9, 10e-9, 5e-9]

    # Récupérer les paramètres depuis le tableau
    params_lp_tchebychev_order3 = FilterBase.get_params("Tchebychev (r=1dB)", 3, "lp")

    print("\n--- Paramètres récupérés depuis le tableau ---")
    for i, (wk_wr, qk, original_wk_wr) in enumerate(params_lp_tchebychev_order3, start=1):
        wk_wr_str = f"{wk_wr:.4f}" if wk_wr is not None else "N/A"
        original_wk_wr_str = f"{original_wk_wr:.4f}" if original_wk_wr is not None else "N/A"
        qk_str = f"{qk:.4f}" if qk is not None else "None"
        print(
            f"  Cellule {i}: original wk/wr = {original_wk_wr_str}, utilisé wk/wr = {wk_wr_str}, Qk = {qk_str}"
        )

    # Calcul des résistances
    result_lp_tchebychev_order3 = filters.tchebychev["lp"].order(f=2500, capacitors=capacitors, order=3)

    print("\n--- Résistances calculées (arrondies au demi-ohm) ---")
    for res_name, res_value in result_lp_tchebychev_order3.items():
        print(f"  {res_name} = {res_value:.1f} Ω")