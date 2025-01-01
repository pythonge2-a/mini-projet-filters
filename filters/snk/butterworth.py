import math
#
class Butterworth_LowPass:
    def __init__(self,order,cutoff_frequency,r1,r2,c1,c2):
        self.order = order
        self.cutoff_frequency = cutoff_frequency
        self.r1 = r1
        self.r2 = r2
        self.c1 = c1
        self.c2 = c2

    def components(order, r1, r2, c1, c2):
        print(r1)
        return None
    
    def graphs(order, cutoff_frequncy):
        return None
    
class HighPass:
    def __init__(self,order,cutoff_frequency):
        self.order = order
        self.cutoff_frequency = cutoff_frequency

    def transfer_function(self, frequency):
        w0 = 2 * math.pi * self.cutoff_frequency
        return 1 / math.sqrt(1 + (frequency / w0) ** (2*self.order))

class BandPass:
    def __init__(self,order,cutoff_frequency):
        self.order = order
        self.cutoff_frequency = cutoff_frequency

    def transfer_function(self, frequency):
        w0 = 2 * math.pi * self.cutoff_frequency
        return 1 / math.sqrt(1 + (frequency / w0) ** (2*self.order))