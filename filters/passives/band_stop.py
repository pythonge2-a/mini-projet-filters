import math


class BandStopFilter:

    @staticmethod
    def order_1_resonant_frequency(R, L):
        """
        Calculate the resonant frequency of a band-stop RL filter (order 1).

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.

        Returns:
            float: Resonant frequency in Hz.
        """
        if R == 0 or L == 0:
            raise ValueError("Resistance and inductance must be non-zero.")
        return R / (2 * math.pi * L)

    @staticmethod
    def order_1_bandwidth(R, L):
        """
        Calculate the bandwidth of a band-stop RL filter (order 1).

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.

        Returns:
            float: Bandwidth in Hz.
        """
        if R == 0 or L == 0:
            raise ValueError("Resistance and inductance must be non-zero.")
        return R / (2 * math.pi * L)

    #The following code is for a Band Pass Filter of order 2
    @staticmethod
    def resonant_frequency(L, C):
        """
        Calculate the resonant frequency of a band-stop RLC filter.

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
    def bandwidth(R, L, C):
        """
        Calculate the bandwidth of a band-stop RLC filter.

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Bandwidth in Hz.
        """
        if L == 0 or C == 0:
            raise ValueError("Inductance and capacitance must be non-zero.")
        return R / (2 * math.pi * L)
