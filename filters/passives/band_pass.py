import math
import numpy as np
import matplotlib.pyplot as plt


class BandPassFilter:
    @staticmethod
    def bandpass_rc(resonant_frequency, bandwidth, resistance=None, capacitance=None):
        """
        Calculate the components (R or C) for a band-pass RC filter of order 1.

        Parameters:
        resonant_frequency (float): Desired resonant frequency in Hz
        bandwidth (float): Desired bandwidth in Hz
        resistance (float, optional): Resistance in ohms (if known)
        capacitance (float, optional): Capacitance in farads (if known)

        Return:
        dict: Calculated component values {"R": value, "C": value}
        """
        omega_bw = 2 * math.pi * bandwidth

        if not resistance and not capacitance:
            raise ValueError("Either resistance or capacitance must be provided.")

        if resistance:
            capacitance = 1 / (omega_bw * resistance)
        elif capacitance:
            resistance = 1 / (omega_bw * capacitance)

        return {"R": resistance, "C": capacitance}

    @staticmethod
    def bandpass_rl(resonant_frequency, bandwidth, resistance=None, inductance=None):
        """
        Calculate the components (R or L) for a band-pass RL filter of order 1.

        Parameters:
        resonant_frequency (float): Desired resonant frequency in Hz
        bandwidth (float): Desired bandwidth in Hz
        resistance (float, optional): Resistance in ohms (if known)
        inductance (float, optional): Inductance in henries (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value}
        """
        omega_bw = 2 * math.pi * bandwidth

        if not resistance and not inductance:
            raise ValueError("Either resistance or inductance must be provided.")

        if resistance:
            inductance = resistance / omega_bw
        elif inductance:
            resistance = omega_bw * inductance

        return {"R": resistance, "L": inductance}

    @staticmethod
    def bandpass_rlc(resonant_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a band-pass RLC filter.

        Parameters:
        resonant_frequency (float): Desired resonant frequency in Hz
        quality_factor (float): Desired quality factor (Q)
        resistance (float, optional): Resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R": value, "L": value, "C": value}
        """
        omega_0 = 2 * math.pi * resonant_frequency

        if not resistance:
            raise ValueError("Resistance must be provided for RLC calculations.")

        C = quality_factor / (omega_0 * resistance)
        L = resistance / (omega_0 * quality_factor)

        return {"R": resistance, "L": L, "C": C}

    @staticmethod
    def bandpass_double_rc(resonant_frequency, bandwidth, resistance=None):
        """
        Calculate the components (R1, R2, C1, C2) for a band-pass RC filter of order 2.

        Parameters:
        resonant_frequency (float): Desired resonant frequency in Hz
        bandwidth (float): Desired bandwidth in Hz
        resistance (float, optional): Shared resistance in ohms (if known)

        Return:
        dict: Calculated component values {"R1": value, "R2": value, "C1": value, "C2": value}
        """
        omega_bw = 2 * math.pi * bandwidth

        if not resistance:
            raise ValueError("Resistance must be provided for RC calculations.")

        R1 = resistance
        R2 = resistance
        C1 = 1 / (omega_bw * R1)
        C2 = 1 / (omega_bw * R2)

        return {"R1": R1, "R2": R2, "C1": C1, "C2": C2}

    @staticmethod
    def bode_plot(resonant_frequency, quality_factor, resistance):
        """
        Plot the Bode diagram (gain and phase) for band-pass filter.

        Parameters:
        resonant_frequency (float): Resonant frequency in Hz
        quality_factor (float): Quality factor (Q)
        resistance (float): Resistance in ohms
        """
        frequencies = np.logspace(1, 6, 500)  # Fréquences de 10 Hz à 1 MHz
        omega = 2 * np.pi * frequencies
        omega_0 = 2 * np.pi * resonant_frequency

        C = quality_factor / (omega_0 * resistance)
        L = resistance / (omega_0 * quality_factor)

        gain = (
            omega
            * resistance
            * C
            / np.sqrt((1 - (omega**2) * L * C) ** 2 + (omega * resistance * C) ** 2)
        )
        phase = -np.arctan((omega * L - 1 / (omega * C)) / resistance) * (180 / np.pi)

        # Plot the amplitude response
        plt.figure(figsize=(10, 8))
        plt.subplot(2, 1, 1)
        plt.semilogx(frequencies, 20 * np.log10(gain))
        plt.axvline(
            resonant_frequency,
            color="red",
            linestyle="--",
            label="Fréquence de résonance",
        )
        plt.title("Diagramme de Bode - Filtre Passe-Bande")
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("Gain (dB)")
        plt.grid(which="both", linestyle="--", linewidth=0.5)
        plt.legend()

        # Plot the phase response
        plt.subplot(2, 1, 2)
        plt.semilogx(frequencies, phase, color="orange", label="Phase")
        plt.axvline(
            resonant_frequency,
            color="red",
            linestyle="--",
            label="Fréquence de résonance",
        )
        plt.xlabel("Fréquence (Hz)")
        plt.ylabel("Phase (degrés)")
        plt.grid(which="both", linestyle="--", linewidth=0.5)
        plt.legend()
        plt.tight_layout()
        plt.show()
