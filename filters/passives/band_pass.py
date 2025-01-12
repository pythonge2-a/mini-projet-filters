import math

class BandPassFilter:
    @staticmethod
    def order_1_resonant_frequency(R, L):
        """
        Calculate the resonant frequency of a band-pass RLC filter order 1.

        Parameters:
            R (float): Resistance in ohms.
            C (float): Capacitance in farads.

        Returns:
            float: Resonant frequency in Hz.
        """
        if L == 0 or R == 0:
            raise ValueError("Inductance and resistance must be non-zero.")
        return R / (2 * math.pi * L)

    #The following code is for a Band Pass Filter of order 2
    @staticmethod
    def resonant_frequency(L, C):
        """
        Calculate the resonant frequency of a band-pass RLC filter.

        Parameters:
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Resonant frequency in Hz.
        """
        if L == 0 or C == 0:
            raise ValueError("Inductance and capacitance must be non-zero.")
        return 1 / (2 * math.pi * math.sqrt(L * C))

    @staticmethod
    def quality_factor(R, L, C):
        """
        Calculate the quality factor of a band-pass RLC filter.

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Quality factor (unitless).
        """
        if L == 0 or C == 0:
            raise ValueError("Inductance and capacitance must be non-zero.")
        return math.sqrt(L / C) / R

    @staticmethod
    def bandwidth(resonant_frequency, quality_factor):
        """
        Calculate the bandwidth of a band-pass RLC filter.

        Parameters:
            resonant_frequency (float): Resonant frequency in Hz.
            quality_factor (float): Quality factor (unitless).

        Returns:
            float: Bandwidth in Hz.
        """
        if quality_factor == 0:
            raise ValueError("Quality factor must be non-zero.")
        return resonant_frequency / quality_factor
