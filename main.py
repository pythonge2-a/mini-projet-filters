from filters.snk.bessel import lowpass


def main():
    # Initialisation de Filters
    highpass_filters = lowpass()
    # Appeler le calcul des composants pour Tchebychev passe-bas d'ordre 2
    c_vals = [1e-6, 1e-9, 1e-6, 1e-9]
    tf, stages = highpass_filters.components(order=4, cutoff_freq=1000, c_vals=c_vals)
    print("Fonction de transfert combinée :", tf)
    for i, stage in enumerate(stages):
        print(f"Stage {i+1}:")
        print("Fonction de transfert :", stage["tf"])
        print("Paramètres :", stage["params"])


if __name__ == "__main__":
    main()
