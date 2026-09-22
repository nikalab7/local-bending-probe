"""Regression for the hybrid's formerly discontinuous twist handover."""
import unittest

import numpy as np

from metrics import HybridAxisBend, helix_on_arc
from perturbation import half_axis_angle, jacobian, rigid_half_bend


class HybridHandoverChecks(unittest.TestCase):
    def test_straight_axis_abstains_between_branches(self):
        metric = HybridAxisBend(9)

        def straight(turn):
            rise = float(np.clip(1.5 + (turn - 100) / 80.0 * 1.8, 1.1, 3.4))
            radius = np.sqrt(3.8 ** 2 - rise ** 2) / (
                2 * np.sin(np.radians(turn) / 2))
            return helix_on_arc(9, radius, rise, turn, 0)

        self.assertEqual(metric.branch(straight(171.35)), "bisector")
        self.assertEqual(metric.branch(straight(171.40)), "abstain")
        self.assertTrue(np.isnan(metric(straight(171.40))))
        self.assertEqual(metric.branch(straight(175)), "abstain")
        self.assertEqual(metric.branch(straight(179)), "smoothed")
        self.assertLess(abs(metric(straight(180))), 1e-6)
        self.assertEqual(metric.branch(straight(183)), "abstain")

    def test_jacobian_uses_measured_axis_change(self):
        P = helix_on_arc(9, 2.3, 1.5, 100.0, 12.0)
        bent = rigid_half_bend(P, 5.0, 45.0)
        self.assertTrue(np.isfinite(half_axis_angle(P)))
        self.assertTrue(np.isfinite(half_axis_angle(bent)))
        result = jacobian(HybridAxisBend(9), P)
        self.assertGreaterEqual(result["n_ok"], 3)
        self.assertTrue(np.isfinite(result["slope"]))


if __name__ == "__main__":
    unittest.main()
