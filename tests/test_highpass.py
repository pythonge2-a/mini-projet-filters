import unittest
from filters.passives.high_pass import HighPassFilter


class TestHighPassFilter(unittest.TestCase):
    def setUp(self):
        self.filter = HighPassFilter()

    def test_highpass_rc(self):
        cutoff_frequency = 1000  # Hz

        # Test with known resistance
        resistance = 1000  # Ohms
        result = self.filter.highpass_rc(cutoff_frequency, resistance=resistance)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["C"], 1.59e-7, places=3)
        print("High-pass RC Test (Resistance Given):", result)

        # Test with known capacitance
        capacitance = 1e-6  # Farads
        result = self.filter.highpass_rc(cutoff_frequency, capacitance=capacitance)
        self.assertIn("R", result)
        self.assertAlmostEqual(result["R"], 159.155, places=3)
        print("High-pass RC Test (Capacitance Given):", result)

    def test_highpass_rl(self):
        cutoff_frequency = 1000  # Hz

        # Test with known resistance
        resistance = 1000  # Ohms
        result = self.filter.highpass_rl(cutoff_frequency, resistance=resistance)
        self.assertIn("L", result)
        self.assertAlmostEqual(result["L"], 0.159, places=3)
        print("High-pass RL Test (Resistance Given):", result)

        # Test with known inductance
        inductance = 1e-3  # Henries
        result = self.filter.highpass_rl(cutoff_frequency, inductance=inductance)
        self.assertIn("R", result)
        self.assertAlmostEqual(result["R"], 6.283, places=3)
        print("High-pass RL Test (Inductance Given):", result)

    def test_highpass_rlc(self):
        cutoff_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.highpass_rlc(
            cutoff_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["L"], 0.225, places=3)
        self.assertAlmostEqual(result["C"], 2.25e-7, places=3)
        print("High-pass RLC Test:", result)

    def test_exceptions(self):
        # Test exception for missing inputs in RC
        with self.assertRaises(ValueError):
            self.filter.highpass_rc(1000)

        # Test exception for missing inputs in RL
        with self.assertRaises(ValueError):
            self.filter.highpass_rl(1000)

        # Test exception for missing resistance in RLC Series
        with self.assertRaises(ValueError):
            self.filter.highpass_rlc(1000, 0.707)


if __name__ == "__main__":
    unittest.main()
