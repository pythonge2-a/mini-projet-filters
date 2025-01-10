import math


class BandStopFilter:
    @staticmethod
    def resonant_frequency_order_1(R, L):
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

    # order 2
    @staticmethod
    def resonant_frequency_order_2(L, C):
        """
        Calculate the resonant frequency of a band-stop RLC filter (order 2).

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
        Calculate the bandwidth of a band-stop RLC filter (order 2).

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Bandwidth in Hz.
        """
        if R == 0 or L == 0 or C == 0:
            raise ValueError(
                "Resistance, inductance, and capacitance must be non-zero."
            )
        return R / (2 * math.pi * L)

    @staticmethod
    def quality_factor_order_2(R, L, C):
        """
        Calculate the quality factor of a band-stop RLC filter (order 2).

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Quality factor (unitless).
        """
        if L == 0 or C == 0 or R == 0:
            raise ValueError(
                "Resistance, inductance, and capacitance must be non-zero."
            )
        resonant_frequency = BandStopFilter.resonant_frequency_order_2(L, C)
        bandwidth = BandStopFilter.bandwidth_order_2(R, L, C)
        return resonant_frequency / bandwidth

    @staticmethod
    def cutoff_frequencies_order_2(L, C, R):
        """
        Calculate the cutoff frequencies (f_c1 and f_c2) for a band-stop RLC filter (order 2).

        Parameters:
            L (float): Inductance in henries.
            C (float): Capacitance in farads.
            R (float): Resistance in ohms.

        Returns:
            tuple: Lower cutoff frequency (f_c1) and upper cutoff frequency (f_c2).
        """
        if L == 0 or C == 0 or R == 0:
            raise ValueError(
                "Inductance, capacitance, and resistance must be non-zero."
            )

        f0 = BandStopFilter.resonant_frequency_order_2(L, C)
        BW = BandStopFilter.bandwidth_order_2(R, L, C)
        f_c1 = f0 - (BW / 2)
        f_c2 = f0 + (BW / 2)
        return f_c1, f_c2
