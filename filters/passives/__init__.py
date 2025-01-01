from snk.bessel import Lowpass

if __name__ == "__main__":
    lowpass_filter = Lowpass()
    tf, params = lowpass_filter.sallen_key_lowpass(2, 1000, r1=1000, r2=1000)
    print("Fonction de transfert :", tf)
    print("Paramètres calculés :", params)

