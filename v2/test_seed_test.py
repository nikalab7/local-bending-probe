"""Focused regressions for seed-test noise and input registration."""
import contextlib
import csv
import io
import os
import tempfile
import unittest

import numpy as np

from seed_test import (bootstrap_tau_ci, calibrated_threshold, parse_ca, run)


class SeedTestChecks(unittest.TestCase):
    def test_systematic_floor_survives_more_crystals(self):
        # One offset per biological variant must survive median aggregation.
        one, _ = calibrated_threshold(1.0, 1, 1, T=80_000, sigma_syst=2.0)
        many, _ = calibrated_threshold(1.0, 20, 20, T=80_000, sigma_syst=2.0)
        iid_only, _ = calibrated_threshold(1.0, 20, 20, T=80_000)
        self.assertGreater(many, 0.75 * one)
        self.assertGreater(many, 3 * iid_only)
        shared, _ = calibrated_threshold(1.0, 20, 20, T=80_000,
                                         sigma_syst=2.0,
                                         wt_forms=["A"] * 20,
                                         mut_forms=["A"] * 20)
        self.assertLess(shared, many / 2)
        self.assertTrue(np.isnan(calibrated_threshold(1.0, 5, 2,
                                                    sigma_syst=float("nan"))[0]))

    def test_cif_is_not_silently_parsed_as_pdb(self):
        with self.assertRaisesRegex(ValueError, "mmCIF parsing is not implemented"):
            parse_ca("example.cif")

    def test_mismatched_mutant_cannot_enter_seed(self):
        with tempfile.TemporaryDirectory() as directory:
            for name, mutant_aa in (("wt1", None), ("wt2", None), ("mut", "GLY")):
                with open(os.path.join(directory, name + ".pdb"), "w") as out:
                    for i in range(10, 19):
                        aa = mutant_aa if i == 14 and mutant_aa else "ALA"
                        x = 0.9 * (-1) ** i
                        z = 3.3 * (i - 10)
                        out.write(f"ATOM   {i:4d}  CA  {aa} A{i:4d}    "
                                  f"{x:8.3f}{0.0:8.3f}{z:8.3f}"
                                  f"{1.0:6.2f}{20.0:6.2f}           C\n")
            manifest = os.path.join(directory, "variants.csv")
            with open(manifest, "w", newline="") as out:
                writer = csv.writer(out)
                writer.writerow(("pdb_id", "chain", "variant", "resnum",
                                 "wt_aa", "mut_aa", "crystal_form"))
                writer.writerow(("wt1", "A", "WT", "", "", "", "P1"))
                writer.writerow(("wt2", "A", "WT", "", "", "", "P1"))
                writer.writerow(("mut", "A", "A14V", "14", "A", "V", "P1"))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                result = run(directory, manifest)
            self.assertIn("1 rejected", output.getvalue())
            self.assertIn("no scorable variants", output.getvalue())
            self.assertEqual(result, {})

            # Correcting the deposited residue restores the pair, but one
            # crystal form cannot identify a systematic noise floor.
            pdb_path = os.path.join(directory, "mut.pdb")
            with open(pdb_path) as source:
                contents = source.read()
            with open(pdb_path, "w") as out:
                out.write(contents.replace("GLY A  14", "VAL A  14"))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                result = run(directory, manifest)
            self.assertIn("0 rejected", output.getvalue())
            self.assertIn("INCONCLUSIVE", output.getvalue())
            self.assertEqual(result["v2_hybrid"]["n"], 1)

    def test_multi_protein_manifest_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = os.path.join(directory, "variants.csv")
            with open(manifest, "w", newline="") as out:
                writer = csv.writer(out)
                writer.writerow(("pdb_id", "chain", "variant", "resnum",
                                 "wt_aa", "mut_aa", "crystal_form", "protein_id"))
                writer.writerow(("a", "A", "WT", "", "", "", "P1", "A"))
                writer.writerow(("b", "A", "WT", "", "", "", "P1", "B"))
            with self.assertRaisesRegex(SystemExit, "single-protein only"):
                run(directory, manifest)

    def test_tau_ci_requires_independent_clusters(self):
        low, high = bootstrap_tau_ci([3, 4, 5], [1, 1, 1], [1, 1, 1],
                                      ["site-a", "site-b", "site-c"],
                                      n_boot=100)
        self.assertLessEqual(low, high)
        self.assertTrue(np.isnan(bootstrap_tau_ci([3, 4], [1, 1], [1, 1],
                                                   ["same", "same"],
                                                   n_boot=100)[0]))


if __name__ == "__main__":
    unittest.main()
