import math


class HighPassFilter:
    @staticmethod
    def highpass_rc(cutoff_frequency, resistance=None, capacitance=None):
        """
        Calculate the components (R or C) for a high-pass RC filter.

        Parameters:
        cutoff_frequency (float): Desired cutoff frequency in Hz
        resistance (float, optional): Resistance in ohms (if known)
        capacitance (float, optional): Capacitance in farads (if known)

        Return:
        dict: Calculated component values {"R": value, "C": value}
        """
        if not resistance and not capacitance:
            raise ValueError("Either resistance or capacitance must be provided.")

        if resistance:
            capacitance = 1 / (2 * math.pi * resistance * cutoff_frequency)
        elif capacitance:
            resistance = 1 / (2 * math.pi * capacitance * cutoff_frequency)

        return {"R": resistance, "C": capacitance}

    @staticmethod
    def highpass_rl(cutoff_frequency, resistance=None, inductance=None):
        """
        Calculate the components (R or L) for a high-pass RL filter.

        Parameters:
        cutoff_frequency (float): Desired cutoff frequency in Hz
        resistance (float, optional): Resistance in ohms (if known)
        inductance (float, optional): Inductance in henries (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value}
        """
        if not resistance and not inductance:
            raise ValueError("Either resistance or inductance must be provided.")

        if resistance:
            inductance = resistance / (2 * math.pi * cutoff_frequency)
        elif inductance:
            resistance = 2 * math.pi * cutoff_frequency * inductance

        return {"R": resistance, "L": inductance}

    @staticmethod
    def highpass_rlc_series(cutoff_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a high-pass RLC filter in series configuration.

        Parameters:
        cutoff_frequency (float): Desired cutoff frequency in Hz
        quality_factor (float): Desired quality factor (Q)
        resistance (float, optional): Resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value, "C": value}
        """
        omega_c = 2 * math.pi * cutoff_frequency

        if not resistance:
            raise ValueError("Resistance must be provided for RLC calculations.")

        L = resistance / (omega_c * quality_factor)
        C = quality_factor / (omega_c * resistance)

        return {"R": resistance, "L": L, "C": C}

    @staticmethod
    def highpass_rlc_parallel(cutoff_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a high-pass RLC filter in parallel configuration.

        Parameters:
        cutoff_frequency (float): Desired cutoff frequency in Hz
        quality_factor (float): Desired quality factor (Q)
        resistance (float, optional): Resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value, "C": value}
        """
        omega_c = 2 * math.pi * cutoff_frequency

        if not resistance:
            raise ValueError("Resistance must be provided for RLC calculations.")

        L = quality_factor * resistance / omega_c
        C = 1 / (omega_c * resistance * quality_factor)

        return {"R": resistance, "L": L, "C": C}
