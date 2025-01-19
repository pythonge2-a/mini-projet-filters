import unittest
from filters.passives.low_pass import LowPassFilter


class TestLowPassFilter(unittest.TestCase):
    def setUp(self):
        self.filter = LowPassFilter()

    def test_lowpass_rc(self):
        cutoff_frequency = 1000  # Hz

        # Test avec résistance donnée
        resistance = 1000  # Ohms
        result = self.filter.lowpass_rc(cutoff_frequency, resistance=resistance)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["C"], 1.59e-7, places=3)
        print("Low-pass RC Test (Resistance Given):", result)

        # Test avec capacité donnée
        capacitance = 1e-6  # Farads
        result = self.filter.lowpass_rc(cutoff_frequency, capacitance=capacitance)
        self.assertIn("R", result)
        self.assertAlmostEqual(result["R"], 159.155, places=3)
        print("Low-pass RC Test (Capacitance Given):", result)

    def test_lowpass_rl(self):
        cutoff_frequency = 1000  # Hz

        # Test avec résistance donnée
        resistance = 1000  # Ohms
        result = self.filter.lowpass_rl(cutoff_frequency, resistance=resistance)
        self.assertIn("L", result)
        self.assertAlmostEqual(result["L"], 0.159, places=3)
        print("Low-pass RL Test (Resistance Given):", result)

        # Test avec inductance donnée
        inductance = 1e-3  # Henries
        result = self.filter.lowpass_rl(cutoff_frequency, inductance=inductance)
        self.assertIn("R", result)
        self.assertAlmostEqual(result["R"], 6.283, places=3)
        print("Low-pass RL Test (Inductance Given):", result)

    def test_lowpass_rlc_series(self):
        cutoff_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.lowpass_rlc_series(
            cutoff_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["L"], 0.1125, places=3)
        self.assertAlmostEqual(result["C"], 2.25e-7, places=3)
        print("Low-pass RLC Series Test:", result)

    def test_lowpass_rlc_parallel(self):
        cutoff_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.lowpass_rlc_parallel(
            cutoff_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["L"], 0.225, places=3)
        self.assertAlmostEqual(result["C"], 1.59e-7, places=3)
        print("Low-pass RLC Parallel Test:", result)

    def test_lowpass_double_rc(self):
        cutoff_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.lowpass_double_rc(
            cutoff_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("R1", result)
        self.assertIn("R2", result)
        self.assertIn("C1", result)
        self.assertIn("C2", result)
        self.assertAlmostEqual(result["C1"], 2.25e-7, places=3)
        self.assertAlmostEqual(result["C2"], 1.59e-7, places=3)
        print("Low-pass Double RC Test:", result)

    def test_exceptions(self):
        # Test exception for missing inputs in RC
        with self.assertRaises(ValueError):
            self.filter.lowpass_rc(1000)

        # Test exception for missing inputs in RL
        with self.assertRaises(ValueError):
            self.filter.lowpass_rl(1000)

        # Test exception for missing resistance in RLC Series
        with self.assertRaises(ValueError):
            self.filter.lowpass_rlc_series(1000, 0.707)

        # Test exception for missing resistance in RLC Parallel
        with self.assertRaises(ValueError):
            self.filter.lowpass_rlc_parallel(1000, 0.707)

        # Test exception for missing resistance in Double RC
        with self.assertRaises(ValueError):
            self.filter.lowpass_double_rc(1000, 0.707)


if __name__ == "__main__":
    unittest.main()
