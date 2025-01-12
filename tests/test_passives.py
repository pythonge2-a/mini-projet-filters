import unittest
from filters.passives import (
    BandPassFilter,
    BandStopFilter,
    LowPassFilter,
    HighPassFilter,
)


class TestPassiveFilters(unittest.TestCase):
    # Tests pour BandPassFilter
    def test_band_pass(self):
        L = 10e-3
        C = 1e-6
        R = 1000

        print("\nTest BandPassFilter:")
        print(f"Entrées : L={L}, C={C}, R={R}")

        resonant_frequency_order1 = BandPassFilter.order_1_resonant_frequency(R, L)
        resonant_frequency_order2 = BandPassFilter.resonant_frequency(L, C)
        quality_factor = BandPassFilter.quality_factor(R, L, C)
        bandwidth = BandPassFilter.bandwidth(resonant_frequency_order2, quality_factor)

        print(f"Fréquence de résonance ordre 1 : {resonant_frequency_order1:.2f} Hz")
        print(f"Fréquence de résonance ordre 2 : {resonant_frequency_order2:.2f} Hz")
        print(f"Facteur de qualité : {quality_factor:.2f}")
        print(f"Bande passante : {bandwidth:.2f} Hz")

        self.assertAlmostEqual(resonant_frequency_order1, 15915.49, places=2)
        self.assertAlmostEqual(resonant_frequency_order2, 1591.55, places=2)
        self.assertAlmostEqual(quality_factor, 0.1, places=2)
        self.assertAlmostEqual(bandwidth, 15915.49, places=2)

    # Tests pour BandStopFilter
    def test_band_stop(self):
        L = 10e-3
        C = 1e-6
        R = 1000

        print("\nTest BandStopFilter:")
        print(f"Entrées : L={L}, C={C}, R={R}")

        resonant_frequency_order1 = BandStopFilter.resonant_frequency_order_1(R, L)
        resonant_frequency_order2 = BandStopFilter.resonant_frequency_order_2(L, C)
        bandwidth_order2 = BandStopFilter.bandwidth(R, L, C)

        print(f"Fréquence de résonance ordre 1 : {resonant_frequency_order1:.2f} Hz")
        print(f"Fréquence de résonance ordre 2 : {resonant_frequency_order2:.2f} Hz")
        print(f"Bande passante ordre 2 : {bandwidth_order2:.2f} Hz")

        self.assertAlmostEqual(resonant_frequency_order1, 15915.49, places=2)
        self.assertAlmostEqual(resonant_frequency_order2, 1591.55, places=2)
        self.assertAlmostEqual(bandwidth_order2, 15915.49, places=2)

    # Tests pour HighPassFilter
    def test_high_pass(self):
        R = 1000
        C = 1e-6
        L = 10e-3

        print("\nTest HighPassFilter:")
        print(f"Entrées RC : R={R}, C={C}")
        print(f"Entrées RL : R={R}, L={L}")

        cutoff_frequency_rc = HighPassFilter.cutoff_frequency_rc(R, C)
        cutoff_frequency_rl = HighPassFilter.cutoff_frequency_rl(R, L)
        cutoff_frequency_rlc = HighPassFilter.cutoff_frequency_rlc(L, C)
        quality_factor = HighPassFilter.quality_factor(R, L, C)
        bandwidth = HighPassFilter.bandwidth(cutoff_frequency_rlc, quality_factor)

        print(f"Fréquence de coupure (RC) : {cutoff_frequency_rc:.2f} Hz")
        print(f"Fréquence de coupure (RL) : {cutoff_frequency_rl:.2f} Hz")
        print(f"Fréquence de coupure (RLC) : {cutoff_frequency_rlc:.2f} Hz")
        print(f"Facteur de qualité : {quality_factor:.2f}")
        print(f"Bande passante : {bandwidth:.2f} Hz")

        self.assertAlmostEqual(cutoff_frequency_rc, 159.15, places=2)
        self.assertAlmostEqual(cutoff_frequency_rl, 15915.49, places=2)
        self.assertAlmostEqual(cutoff_frequency_rlc, 1591.55, places=2)
        self.assertAlmostEqual(quality_factor, 0.1, places=3)
        self.assertAlmostEqual(bandwidth, 15915.49, places=2)

    # Tests pour LowPassFilter
    def test_low_pass(self):
        R = 1000
        C = 1e-6
        R1 = 1000
        R2 = 2000
        C1 = 1e-6
        C2 = 2e-6
        L = 10e-3

        print("\nTest LowPassFilter:")
        print(f"Entrées RC : R={R}, C={C}")
        print(f"Entrées R1R2C1C2 : R1={R1}, R2={R2}, C1={C1}, C2={C2}")
        print(f"Entrées RLC : L={L}, C={C}")

        cutoff_frequency_rc = LowPassFilter.cutoff_frequency_rc(R, C)
        cutoff_frequency_rlc = LowPassFilter.cutoff_frequency_rlc(L, C)
        quality_factor = LowPassFilter.quality_factor(R, L, C)
        bandwidth_rlc = LowPassFilter.bandwidth(cutoff_frequency_rlc, quality_factor)

        print(f"Fréquence de coupure (RC) : {cutoff_frequency_rc:.2f} Hz")
        print(f"Fréquence de coupure (ordre 2) : {cutoff_frequency_rlc:.2f} Hz")
        print(f"Fréquence de coupure (RLC) : {cutoff_frequency_rlc:.2f} Hz")
        print(f"Facteur de qualité (RLC) : {quality_factor:.2f}")
        print(f"Bande passante (RLC) : {bandwidth_rlc:.2f} Hz")

        self.assertAlmostEqual(cutoff_frequency_rc, 159.15, places=2)
        self.assertAlmostEqual(cutoff_frequency_rlc, 1591.55, places=2)
        self.assertAlmostEqual(quality_factor, 0.1, places=3)
        self.assertAlmostEqual(bandwidth_rlc, 15915.49, places=2)


if __name__ == "__main__":
    unittest.main()
