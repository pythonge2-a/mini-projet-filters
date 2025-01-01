import math

class LowPassFilter:
    @staticmethod
    def cutoff_frequency_order_1(R, C):
        """
        Calculate the cutoff frequency of a low-pass RLC filter order 1.

        Parameters:
        R (float): Resistance
        C (float): Capacitance

        Return:
        float: Cutoff frequency
        """

        if C <= 0 or R <= 0:
            raise ValueError("Capacitance and resistance must be positive.")
        return 1 / (2 * math.pi * R * C)

    @staticmethod
    def cutoff_frequency_order_2(R1, R2, C1, C2):
        """
        Calculate the cutoff frequency of a low-pass RLC filter order 2.

        Parameters:
        R1 (float): Resistance 1
        R2 (float): Resistance 1
        C1 (float): Capacitance 1
        C2 (float): Capacitance 2

        Return:
        float: Cutoff frequency
        """

        if R1 <= 0 or R2 <= 0 or C1 <= 0 or C2 <= 0:
            raise ValueError("Resistances and capacitances must be positive.")
        return 1 / (2 * math.pi * math.sqrt(R1 * R2 * C1 * C2))
