import math
import numpy as np
import matplotlib.pyplot as plt


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
    def highpass_rlc(cutoff_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a high-pass RLC filter.

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
    def bode_plot(
        cutoff_frequency,
        resistance,
        inductance=None,
        capacitance=None,
        filter_type="RC",
    ):
        """
        Plot the Bode diagram (gain and phase) for high-pass RC, RL, and RLC filters.

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
            gain = (omega * resistance * capacitance) / np.sqrt(
                1 + (omega * resistance * capacitance) ** 2
            )
            phase = np.arctan(
                1 / (omega * resistance * capacitance)
            )  # Phase in radians
            title = "Filtre Passe-Haut RC"

        elif filter_type == "RL":
            if inductance is None:
                raise ValueError("Inductance must be provided for RL filter.")
            gain = (omega * inductance / resistance) / np.sqrt(
                1 + (omega * inductance / resistance) ** 2
            )
            phase = np.arctan(resistance / (omega * inductance))  # Phase in radians
            title = "Filtre Passe-Haut RL"

        elif filter_type == "RLC":
            if inductance is None or capacitance is None:
                raise ValueError(
                    "Inductance and capacitance must be provided for RLC filter."
                )
            gain = (
                omega**2
                * inductance
                * capacitance
                / np.sqrt(
                    (1 - (omega**2) * inductance * capacitance) ** 2
                    + (omega * resistance * capacitance) ** 2
                )
            )
            phase = np.arctan(
                (omega * resistance * capacitance)
                / (1 - omega**2 * inductance * capacitance)
            )  # Phase in radians
            title = "Filtre Passe-Haut RLC"

        else:
            raise ValueError("Invalid filter type. Choose 'RC', 'RL', or 'RLC'.")

        # Plot the amplitude response
        plt.figure(figsize=(10, 6))
        plt.subplot(2, 1, 1)
        plt.semilogx(frequencies, 20 * np.log10(gain), label="Gain")
        plt.axvline(
            cutoff_frequency, color="red", linestyle="--", label="Fréquence de coupure"
        )
        plt.title(f"Diagramme de Bode - {title}")
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("Gain (dB)")
        plt.grid(which="both", linestyle="--", linewidth=0.5)
        plt.legend()

        # Plot the phase response
        plt.subplot(2, 1, 2)
        plt.semilogx(frequencies, np.degrees(phase), label="Phase", color="orange")
        plt.axvline(
            cutoff_frequency, color="red", linestyle="--", label="Fréquence de coupure"
        )
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("Phase (degrés)")
        plt.grid(which="both", linestyle="--", linewidth=0.5)
        plt.legend()
        plt.tight_layout()
        plt.show()
