import bessel as bt

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