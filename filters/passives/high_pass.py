import math


class HighPassFilter:
    @staticmethod
    def cutoff_frequency_rc(R, C):
        """
        Calculate the cutoff frequency of a high-pass RC filter.

        Parameters:
            R (float): Resistance in ohms.
            C (float): Capacitance in farads.

        Returns:
            float: Cutoff frequency in Hz.
        """
        if R <= 0 or C <= 0:
            raise ValueError("Resistance and capacitance must be greater than zero.")
        return 1 / (2 * math.pi * R * C)

    @staticmethod
    def cutoff_frequency_rl(R, L):
        """
        Calculate the cutoff frequency of a high-pass RL filter.

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.

        Returns:
            float: Cutoff frequency in Hz.
        """
        if R <= 0 or L <= 0:
            raise ValueError("Resistance and inductance must be greater than zero.")
        return R / (2 * math.pi * L)

    @staticmethod
    def cutoff_frequency_rlc(L, C):
        """
        Calculate the cutoff frequency of a high-pass RLC filter.

        Parameters:
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Cutoff frequency in Hz.
        """
        if L <= 0 or C <= 0:
            raise ValueError("Inductance and capacitance must be greater than zero.")
        return 1 / (2 * math.pi * math.sqrt(L * C))

    @staticmethod
    def quality_factor(R, L, C):
        """
        Calculate the quality factor of a high-pass RLC filter.

        Parameters:
            R (float): Resistance in ohms.
            L (float): Inductance in henries.
            C (float): Capacitance in farads.

        Returns:
            float: Quality factor (unitless).
        """
        if R <= 0 or L <= 0 or C <= 0:
            raise ValueError(
                "Resistance, inductance, and capacitance must be greater than zero."
            )
        return math.sqrt(L / C) / R

    @staticmethod
    def bandwidth(cutoff_frequency, quality_factor):
        """
        Calculate the bandwidth of a high-pass RLC filter.

        Parameters:
            cutoff_frequency (float): Cutoff frequency in Hz.
            quality_factor (float): Quality factor (unitless).

        Returns:
            float: Bandwidth in Hz.
        """
        if quality_factor <= 0:
            raise ValueError("Quality factor must be greater than zero.")
        return cutoff_frequency / quality_factor
