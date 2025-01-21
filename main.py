from filters import Filters

def main():
    # Initialisation de Filters
    filters = Filters()
        # Appeler le calcul des composants pour Tchebychev passe-bas d'ordre 2
    result = filters.snk.tchebychev.2.components(1000, 2000, 1e-6, 2e-6)
    print("Résultat des composants :", result)

if __name__ == "__main__":
    main()
