# __init__.py
import butterworth as bt
import bessel as bl


lowpass = bl.lowpass()
highpass = bl.highpass()

# Spécifier les résistances et calculer les condensateurs
r_vals = [1000,10000]
highpass.graphs(order=2, cutoff_freq=1000, r_vals=r_vals)




