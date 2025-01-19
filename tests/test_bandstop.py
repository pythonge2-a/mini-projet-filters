import unittest
from filters.passives.band_stop import BandStopFilter


class TestBandStopFilter(unittest.TestCase):
    def setUp(self):
        self.filter = BandStopFilter()

    def test_bandstop_rlc_series(self):
        center_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.bandstop_rlc_series(
            center_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["L"], 0.1125, places=3)
        self.assertAlmostEqual(result["C"], 1.1252e-7, places=3)
        print("Band-stop RLC Series Test:", result)

    def test_bandstop_rlc_parallel(self):
        center_frequency = 1000  # Hz
        quality_factor = 0.707
        resistance = 1000  # Ohms

        result = self.filter.bandstop_rlc_parallel(
            center_frequency, quality_factor, resistance=resistance
        )
        self.assertIn("L", result)
        self.assertIn("C", result)
        self.assertAlmostEqual(result["L"], 0.2251, places=3)
        self.assertAlmostEqual(result["C"], 1.5916e-7, places=3)
        print("Band-stop RLC Parallel Test:", result)

    def test_exceptions(self):
        center_frequency = 1000  # Hz
        quality_factor = 0.707

        with self.assertRaises(ValueError):
            self.filter.bandstop_rlc_series(center_frequency, quality_factor)

        with self.assertRaises(ValueError):
            self.filter.bandstop_rlc_parallel(center_frequency, quality_factor)


if __name__ == "__main__":
    unittest.main()
