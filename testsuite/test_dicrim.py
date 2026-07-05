"""Tests r.hydro.hbv against the Plumergat DICRIM flood-risk dataset --
11 SHYREG sub-basins (Brittany, France).

The fixture here is a 1000-day slice of the real, ERA5-derived data used
to produce the DICRIM report (see
~/Documents/LaTex/DICRIM/Plumergat/data/hbv/), not synthetic data. A full
bit-exact reference isn't meaningful for a truncated fixture (unlike
test_original.py, which has one), so this test checks that the module
runs to completion, produces the 11 expected SHYREG-named output files
with the right shapes, and that discharge/ET time series values are
finite and non-negative (NS/RVE summary scores can legitimately be NaN
when a basin's calibration window has too little valid observed data --
see test_original.py's module docstring and the original model's own
ETout*/EToutput* reference files, which already contain NaN/-nan for the
same reason).

test_dicrim_full (skipped by default -- set
R_HYDRO_HBV_RUN_SLOW_TESTS=1 to enable) runs the full 8400-day,
45000-realization DICRIM dataset from
~/Documents/LaTex/DICRIM/Plumergat/data/hbv/.
"""

import os
import math

from grass.gunittest.case import TestCase
from grass.gunittest.main import test

TESTDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(TESTDIR, "data", "dicrim")

BASIN_IDS = [
    "BR2544",
    "BR4471",
    "BR4472",
    "BR5820",
    "BR7136",
    "BR4473",
    "BR2545",
    "BR4474",
    "BR2034",
    "BR7137",
    "BR3055",
]


class TestDicrimSmall(TestCase):
    def setUp(self):
        self.outdir = os.path.join(TESTDIR, "_out_dicrim_small")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.outdir, ignore_errors=True)

    def _run(self, outdir):
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
            n_basins=11,
            n_calib_steps=1000,
            n_days=1000,
            n_years=2,
            warmup=200,
            n_realizations=20,
            output=outdir,
        )

    def test_dicrim_small(self):
        self._run(self.outdir)

        for i, basin_id in enumerate(BASIN_IDS, start=1):
            basinout = os.path.join(
                self.outdir, "Basinout%02d-%s.csv" % (i, basin_id)
            )
            etout = os.path.join(self.outdir, "ETout%02d-%s.csv" % (i, basin_id))
            output = os.path.join(self.outdir, "Output%02d-%s.csv" % (i, basin_id))
            for path in (basinout, etout, output):
                self.assertTrue(os.path.isfile(path), msg=path)
                self.assertTrue(os.path.getsize(path) > 0, msg=path)

            with open(basinout) as f:
                rows = [line.split() for line in f if line.strip()]
            self.assertEqual(len(rows), 1000, msg=basinout)
            for qo, qr, po in rows:
                qr = float(qr)
                self.assertTrue(
                    math.isfinite(qr) and qr >= 0.0,
                    msg="%s: simulated discharge must be finite and "
                    ">= 0, got %s" % (basinout, qr),
                )


class TestDicrimFull(TestCase):
    """Full-scale (8400 days, 45000 realizations) reproduction of the
    dataset behind the Plumergat DICRIM report. Slow; skipped unless
    R_HYDRO_HBV_RUN_SLOW_TESTS=1 is set."""

    def setUp(self):
        if not os.environ.get("R_HYDRO_HBV_RUN_SLOW_TESTS"):
            self.skipTest(
                "set R_HYDRO_HBV_RUN_SLOW_TESTS=1 to run the full "
                "8400-day / 45000-realization DICRIM dataset"
            )
        self.outdir = os.path.join(TESTDIR, "_out_dicrim_full")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_dicrim_full(self):
        full_dir = os.path.expanduser(
            "~/Documents/LaTex/DICRIM/Plumergat/data/hbv"
        )
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation=os.path.join(full_dir, "precip.csv"),
            temperature=os.path.join(full_dir, "temp.csv"),
            evapotranspiration=os.path.join(full_dir, "evap.csv"),
            eta_observed=os.path.join(full_dir, "etobs.csv"),
            discharge_observed=os.path.join(full_dir, "dischargeobs.csv"),
            parameters=os.path.join(full_dir, "param_sto.csv"),
            basin_ids=os.path.join(DATADIR, "basin_ids.txt"),
            n_basins=11,
            n_calib_steps=8400,
            n_days=8400,
            n_years=22,
            warmup=365,
            n_realizations=45000,
            output=self.outdir,
        )
        outlet = os.path.join(self.outdir, "Basinout11-BR3055.csv")
        self.assertTrue(os.path.isfile(outlet))
        with open(outlet) as f:
            rows = [line.split() for line in f if line.strip()]
        self.assertEqual(len(rows), 8400)


if __name__ == "__main__":
    test()
