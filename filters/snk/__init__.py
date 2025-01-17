import bessel as bt
import butterworth as btw

'''
lowpass = bt.lowpass()
highpass=bt.highpass()
r_vals = [1000,10000,1000,10000]
#lowpass.graphs(order =1, cutoff_freq=1000,r_vals=r_vals)

tf,stages = highpass.multiple_order_highpass(order=4, cutoff_freq=1000, r_vals=r_vals)
print("Fonction de transfert combinée :", tf)
for i, stage in enumerate(stages):
    print(f"Stage {i+1}:")
    print("Fonction de transfert :", stage['tf'])
    print("Paramètres :", stage['params'])
'''
r_vals2 = [2.2e-6, 2.2e-6, 2.2e-9, 2.2e-6, 2.2e-9]
lowpass=btw.Butterworth_LowPass()
values = lowpass.components(order=5, cutoff_frequency=10000, condo_values=r_vals2)
print(values)