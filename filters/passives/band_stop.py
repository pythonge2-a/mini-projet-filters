import math


class BandStopFilter:
    @staticmethod
    def bandstop_rlc_series(center_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a band-stop RLC filter in series configuration.

        Parameters:
        center_frequency (float): Desired center frequency in Hz
        quality_factor (float): Desired quality factor (Q)
        resistance (float, optional): Resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value, "C": value}
        """
        omega_0 = 2 * math.pi * center_frequency

        if not resistance:
            raise ValueError("Resistance must be provided for RLC calculations.")

        L = quality_factor * resistance / omega_0
        C = 1 / (omega_0**2 * L)

        return {"R": resistance, "L": L, "C": C}

    @staticmethod
    def bandstop_rlc_parallel(center_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a band-stop RLC filter in parallel configuration.

        Parameters:
        center_frequency (float): Desired center frequency in Hz
        quality_factor (float): Desired quality factor (Q)
        resistance (float, optional): Resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value, "C": value}
        """
        omega_0 = 2 * math.pi * center_frequency

        if not resistance:
            raise ValueError("Resistance must be provided for RLC calculations.")

        C = quality_factor / (omega_0 * resistance)
        L = resistance / (omega_0 * quality_factor)

        return {"R": resistance, "L": L, "C": C}
