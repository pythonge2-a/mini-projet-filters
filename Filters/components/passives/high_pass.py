import math


class HighPassFilter:
    @staticmethod
    def cutoff_frequency(R, L, C):
        """
        Calculate the cutoff frequency of a high-pass RLC filter.

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Cutoff frequency in Hz.
        """
        if L == 0:
            raise ValueError("Inductance must be non-zero.")
        return R / (2 * math.pi * L)
