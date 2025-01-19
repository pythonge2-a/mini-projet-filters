import math
import numpy as np
import matplotlib.pyplot as plt


class BandStopFilter:
    @staticmethod
    def bandstop_rlc(center_frequency, quality_factor, resistance=None):
        """
        Calculate the components (R, L, C) for a band-stop RLC filter.

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

    @staticmethod
    def bode_plot_bandstop(resonant_frequency, resistance, inductance, capacitance):
        """
        Plot the Bode diagram for a band-stop filter (RLC Parallel or Series).

        Parameters:
        resonant_frequency (float): Resonant frequency in Hz
        resistance (float): Resistance in ohms
        inductance (float): Inductance in henries
        capacitance (float): Capacitance in farads
        """
        frequencies = np.logspace(1, 6, 500)  # Fréquences de 10 Hz à 1 MHz
        omega = 2 * np.pi * frequencies
        gain = np.sqrt(1 + (omega**2 * inductance * capacitance) ** 2) / np.sqrt(
            (1 + omega**2 * inductance * capacitance) ** 2
            + (omega * resistance * capacitance) ** 2
        )
        phase = np.arctan2(
            omega * resistance * capacitance,
            1 - omega**2 * inductance * capacitance,
        ) * (180 / np.pi)

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
        plt.title("Diagramme de Bode - Filtre Coupe-Bande")
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
