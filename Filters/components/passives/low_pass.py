import math


class LowPassFilter:
    @staticmethod
    def cutoff_frequency(R, L, C):
        """
        Calculate the cutoff frequency of a low-pass RLC filter.

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Cutoff frequency in Hz.
        """
        if C == 0 or L == 0:
            raise ValueError("Capacitance and inductance must be non-zero.")
        return 1 / (2 * math.pi * math.sqrt(L * C))
