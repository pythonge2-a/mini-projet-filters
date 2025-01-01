# __init__.py
import butterworth as bt
import bessel as bl

lowpass = bl.Lowpass()
highpass = bl.Highpass()

# Spécifier les résistances et calculer les condensateurs
r_vals = [3300,1000,1000,2200,1000,1000]
tf, stages = lowpass.multiple_order_lowpass(order=6, cutoff_freq=20000, r_vals=r_vals)

print("Fonction de transfert combinée :", tf)
for i, stage in enumerate(stages):
    print(f"Stage {i+1}:")
    print("Fonction de transfert :", stage['tf'])
    print("Paramètres :", stage['params'])