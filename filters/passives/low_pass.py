import math
import numpy as np
import matplotlib.pyplot as plt


class LowPassFilter:
    @staticmethod
    def lowpass_rc(cutoff_frequency, resistance=None, capacitance=None):
        """
        Calculate the components (R or C) for a low-pass RC filter.

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
    def lowpass_rl(cutoff_frequency, resistance=None, inductance=None):
        """
        Calculate the components (R or L) for a low-pass RL filter.

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
    def lowpass_rlc(cutoff_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a low-pass RLC filter.

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

    @staticmethod
    def lowpass_double_rc(cutoff_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R1, R2, C1, C2) for a low-pass RC filter of order 2.

        Parameters:
        cutoff_frequency (float): Desired cutoff frequency in Hz
        quality_factor (float): Desired quality factor (Q)
        resistance (float, optional): Shared resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R1": value, "R2": value, "C1": value, "C2": value}
        """
        omega_c = 2 * math.pi * cutoff_frequency

        if not resistance:
            raise ValueError("Resistance must be provided for RC calculations.")

        R1 = resistance
        R2 = resistance
        C1 = 1 / (omega_c * R1 * quality_factor)
        C2 = quality_factor / (omega_c * R2)

        return {"R1": R1, "R2": R2, "C1": C1, "C2": C2}

    @staticmethod
    def bode_plot(
        cutoff_frequency,
        resistance,
        inductance=None,
        capacitance=None,
        filter_type="RC",
    ):
        """
        Plot the Bode diagram for low-pass RC, RL, and RLC filters.

        Parameters:
        cutoff_frequency (float): Desired cutoff frequency in Hz
        resistance (float): Resistance in ohms
        inductance (float, optional): Inductance in henries (for RL or RLC filters)
        capacitance (float, optional): Capacitance in farads (for RC or RLC filters)
        filter_type (str): Type of the filter ("RC", "RL", "RLC")
        """
        frequencies = np.logspace(1, 6, 500)  # Frequency range: 10 Hz to 1 MHz
        omega = 2 * np.pi * frequencies

        if filter_type == "RC":
            if capacitance is None:
                raise ValueError("Capacitance must be provided for RC filter.")
            response = 1 / np.sqrt(1 + (omega * resistance * capacitance) ** 2)
            title = "Filtre Passe-Bas RC"

        elif filter_type == "RL":
            if inductance is None:
                raise ValueError("Inductance must be provided for RL filter.")
            response = 1 / np.sqrt(1 + ((omega * inductance) / resistance) ** 2)
            title = "Filtre Passe-Bas RL"

        elif filter_type == "RLC":
            if inductance is None or capacitance is None:
                raise ValueError(
                    "Inductance and capacitance must be provided for RLC filter."
                )
            numerator = 1
            denominator = np.sqrt(
                1
                + (resistance * capacitance * omega) ** 2
                - (omega**2) * inductance * capacitance
            )
            response = numerator / denominator
            title = "Filtre Passe-Bas RLC Série"

        else:
            raise ValueError("Invalid filter type. Choose 'RC', 'RL', or 'RLC'.")

        # Plot the amplitude response
        plt.figure(figsize=(10, 6))
        plt.semilogx(frequencies, 20 * np.log10(response))
        plt.axvline(
            cutoff_frequency, color="red", linestyle="--", label="Fréquence de coupure"
        )
        plt.title(f"Diagramme de Bode - {title}")
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("Gain (dB)")
        plt.grid(which="both", linestyle="--", linewidth=0.5)
        plt.legend()
        plt.show()
