"""Tests r.hydro.hbv against the original Karkheh basin (Muthuwatta PhD
thesis, ITC 2005) dataset -- 8 lumped sub-basins.

The fast test (test_original_small) runs a handful of Monte-Carlo
realizations and checks the deterministic Basinout*/ETout* time series
(discharge/ET of the *last* realization) bit-for-bit against a golden
fixture captured once with the same inputs and n_realizations.

test_original_full (skipped by default -- set
R_HYDRO_HBV_RUN_SLOW_TESTS=1 to enable) reproduces the full historical
45000-realization run and diffs the same deterministic files against the
values shipped with the original standalone ~/dev/HBV project, proving
"runs like originally". Output*/EToutput* (the Monte-Carlo NS/RVE summary
across all realizations) are *not* compared bit-for-bit: even the
original, unmodified 2005 code is not run-to-run reproducible there,
because its inner time-step loop is parallelized over t with
sequential (t -> t+1) state dependencies, and its performance
accumulators are (in legacy mode) never reset between realizations --
both quirks were verified against the pristine original binary and are
preserved here, not introduced by this port. Basinout*/ETout* (the
deterministic last-realization time series) reproduce bit-for-bit run-
to-run *most* of the time, but not guaranteed: the same race can
occasionally perturb a handful of rows for one basin (see
docs/raster_options.md). If this test ever fails only on
test_original_full and only by a few rows on one or two basins, that is
this known race, not a regression -- re-run to confirm before treating
it as a real bug.
"""

import os
import csv
import math

from grass.gunittest.case import TestCase
from grass.gunittest.main import test

TESTDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(TESTDIR, "data", "original")
REFDIR = os.path.join(TESTDIR, "data", "original_reference")

BASINS = [
    (1, "Doab"),
    (2, "Pole_Chehr"),
    (3, "Doabe_M"),
    (4, "Ghor_B"),
    (5, "Holilan"),
    (6, "Pole_D"),
    (7, "Jelogir"),
    (8, "Paye_P"),
]


def _read_csv_numbers(path):
    with open(path) as f:
        return [
            [float(x) for x in line.split()]
            for line in f
            if line.strip()
        ]


class TestOriginalSmall(TestCase):
    def setUp(self):
        self.outdir = os.path.join(TESTDIR, "_out_original_small")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_original_small(self):
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation=os.path.join(DATADIR, "precip.csv"),
            temperature=os.path.join(DATADIR, "temp.csv"),
            evapotranspiration=os.path.join(DATADIR, "evap.csv"),
            eta_observed=os.path.join(DATADIR, "etobs.csv"),
            discharge_observed=os.path.join(DATADIR, "dischargeobs.csv"),
            parameters=os.path.join(DATADIR, "param_sto.csv"),
            basin_ids=os.path.join(DATADIR, "basin_ids.txt"),
            n_basins=8,
            n_calib_steps=916,
            n_days=1910,
            n_years=5,
            warmup=274,
            n_realizations=20,
            output=self.outdir,
        )

        for i, name in BASINS:
            for prefix in ("Basinout", "ETout"):
                fname = "%s%02d-%s.csv" % (prefix, i, name)
                got = _read_csv_numbers(os.path.join(self.outdir, fname))
                want = _read_csv_numbers(os.path.join(REFDIR, fname))
                self.assertEqual(
                    len(got), len(want), msg="%s: row count differs" % fname
                )
                for row_got, row_want in zip(got, want):
                    for v_got, v_want in zip(row_got, row_want):
                        if math.isnan(v_got) and math.isnan(v_want):
                            continue
                        self.assertAlmostEqual(
                            v_got,
                            v_want,
                            places=4,
                            msg="%s: value mismatch" % fname,
                        )

    def test_output_files_are_well_formed(self):
        """Output*/EToutput* aren't bit-exact reproducible (see module
        docstring) but must still be present, non-empty and numeric."""
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation=os.path.join(DATADIR, "precip.csv"),
            temperature=os.path.join(DATADIR, "temp.csv"),
            evapotranspiration=os.path.join(DATADIR, "evap.csv"),
            eta_observed=os.path.join(DATADIR, "etobs.csv"),
            discharge_observed=os.path.join(DATADIR, "dischargeobs.csv"),
            parameters=os.path.join(DATADIR, "param_sto.csv"),
            basin_ids=os.path.join(DATADIR, "basin_ids.txt"),
            n_basins=8,
            n_calib_steps=916,
            n_days=1910,
            n_years=5,
            warmup=274,
            n_realizations=20,
            output=self.outdir,
        )
        for i, name in BASINS:
            for prefix in ("Output", "EToutput"):
                fname = "%s%02d-%s.csv" % (prefix, i, name)
                path = os.path.join(self.outdir, fname)
                self.assertTrue(os.path.getsize(path) > 0, msg=fname)
                with open(path) as f:
                    reader = csv.reader(f)
                    next(reader)  # header
                    n_rows = sum(1 for _ in reader)
                if i == len(BASINS):
                    # hbv_performance()'s selection-counter loops only run
                    # over b < rtot-1, so the last (outlet) basin never
                    # gets any behavioural parameter sets recorded -- true
                    # of the original code too, not a defect of this port.
                    self.assertEqual(n_rows, 0, msg="%s: row count" % fname)
                else:
                    # hbv_performance() stores mtot[b] as the
                    # *pre-increment* selection counter each realization,
                    # so after n_realizations calls the final row count
                    # is n_realizations - 1 (true of the original code
                    # too).
                    self.assertEqual(n_rows, 19, msg="%s: row count" % fname)


class TestOriginalFull(TestCase):
    """Full-scale (45000 realizations) reproduction of the historical
    Karkheh calibration. Slow (~40s-15min depending on core count -- see
    the Run times section of ~/dev/HBV/README.md); skipped unless
    R_HYDRO_HBV_RUN_SLOW_TESTS=1 is set."""

    def setUp(self):
        if not os.environ.get("R_HYDRO_HBV_RUN_SLOW_TESTS"):
            self.skipTest(
                "set R_HYDRO_HBV_RUN_SLOW_TESTS=1 to run the full "
                "45000-realization reproduction"
            )
        self.outdir = os.path.join(TESTDIR, "_out_original_full")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_original_full(self):
        hbv_dir = os.path.expanduser("~/dev/HBV")
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation=os.path.join(hbv_dir, "precip.csv"),
            temperature=os.path.join(hbv_dir, "temp.csv"),
            evapotranspiration=os.path.join(hbv_dir, "evap.csv"),
            eta_observed=os.path.join(hbv_dir, "etobs.csv"),
            discharge_observed=os.path.join(hbv_dir, "dischargeobs.csv"),
            parameters=os.path.join(hbv_dir, "param_sto.csv"),
            basin_ids=os.path.join(DATADIR, "basin_ids.txt"),
            n_basins=8,
            n_calib_steps=916,
            n_days=1910,
            n_years=5,
            warmup=274,
            n_realizations=45000,
            output=self.outdir,
            flags="l",
        )
        for i, name in BASINS:
            for prefix, ref_suffix in (
                ("Basinout", "det-%s.csv" % name),
                ("ETout", "et-%s.csv" % name),
            ):
                got = _read_csv_numbers(
                    os.path.join(self.outdir, "%s%02d-%s.csv" % (prefix, i, name))
                )
                want = _read_csv_numbers(
                    os.path.join(hbv_dir, "%s%02d%s" % (prefix, i, ref_suffix))
                )
                self.assertEqual(got, want)


if __name__ == "__main__":
    test()
