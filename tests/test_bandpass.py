import unittest
from filters.passives.band_pass import BandPassFilter
import math


class TestBandPassFilter(unittest.TestCase):
    def setUp(self):
        self.filter = BandPassFilter()

    def test_bandpass_rl(self):
        resonant_frequency = 1000  # Hz
        bandwidth = 200  # Hz

        # Test avec résistance donnée
        resistance = 1000  # Ohms
        result = self.filter.bandpass_rl(
            resonant_frequency, bandwidth, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertAlmostEqual(
            result["L"], resistance / (2 * math.pi * bandwidth), places=3
        )
        print("Band-pass RL Test (Resistance Given):", result)

        # Test avec inductance donnée
        inductance = 1e-3  # Henries
        result = self.filter.bandpass_rl(
            resonant_frequency, bandwidth, inductance=inductance
        )
        self.assertIn("R", result)
        self.assertAlmostEqual(
            result["R"], 2 * math.pi * bandwidth * inductance, places=3
        )
        print("Band-pass RL Test (Inductance Given):", result)

    def test_bandpass_rlc(self):
        resonant_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.bandpass_rlc(
            resonant_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertIn("C", result)
        omega_0 = 2 * math.pi * resonant_frequency
        self.assertAlmostEqual(
            result["C"], quality_factor / (omega_0 * resistance), places=3
        )
        self.assertAlmostEqual(
            result["L"], resistance / (omega_0 * quality_factor), places=3
        )
        print("Band-pass RLC Test:", result)

    def test_bandpass_double_rc(self):
        resonant_frequency = 1000  # Hz
        bandwidth = 200  # Hz
        resistance = 1000  # Ohms

        result = self.filter.bandpass_double_rc(
            resonant_frequency, bandwidth, resistance=resistance
        )
        self.assertIn("R1", result)
        self.assertIn("R2", result)
        self.assertIn("C1", result)
        self.assertIn("C2", result)
        omega_bw = 2 * math.pi * bandwidth
        self.assertAlmostEqual(result["C1"], 1 / (omega_bw * resistance), places=10)
        self.assertAlmostEqual(result["C2"], 1 / (omega_bw * resistance), places=10)
        print("Band-pass Double RC Test:", result)

    def test_bandpass_rc(self):
        resonant_frequency = 1000  # Hz
        bandwidth = 200  # Hz

        # Test avec résistance donnée
        resistance = 1000  # Ohms
        result = self.filter.bandpass_rc(
            resonant_frequency, bandwidth, resistance=resistance
        )
        self.assertIn("C", result)
        omega_bw = 2 * math.pi * bandwidth
        self.assertAlmostEqual(result["C"], 1 / (omega_bw * resistance), places=10)
        print("Band-pass RC Test (Resistance Given):", result)

        # Test avec capacité donnée
        capacitance = 1e-6  # Farads
        result = self.filter.bandpass_rc(
            resonant_frequency, bandwidth, capacitance=capacitance
        )
        self.assertIn("R", result)
        self.assertAlmostEqual(result["R"], 1 / (omega_bw * capacitance), places=6)
        print("Band-pass RC Test (Capacitance Given):", result)

    def test_exceptions(self):
        # Test exception for missing inputs in RL
        with self.assertRaises(ValueError):
            self.filter.bandpass_rl(1000, 200)


if __name__ == "__main__":
    unittest.main()
