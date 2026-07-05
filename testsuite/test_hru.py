"""Tests r.hydro.hbv's semi-distributed elevation-band HRU support
(hru_bands= on r.hydro.hbv.basins, read through basins_input.c's
station_of[]/rtot grouping): each basin's HRUs get independent state and
Monte-Carlo parameter draws (btot computation units), but discharge/ETa
are aggregated back to the reporting-station level (rtot, matching
observed data) before scoring and reporting -- see
docs/raster_options.md and the Context section of the approved plan for
why this split exists and what hbv_performance.c's own btot-1-vs-rtot-1
loops require.

hru_bands=1 (the default, exercised by test_original.py/test_dicrim.py/
test_basins_integration.py/test_table_io.py) must remain an exact no-op;
this file only covers the hru_bands>1 path specifically.
"""

import os
import shutil

from grass.gunittest.case import TestCase
from grass.gunittest.gmodules import call_module
from grass.gunittest.main import test

TESTDIR = os.path.dirname(os.path.abspath(__file__))

N_DAYS = 20
STATIONS = ["BasinA", "BasinB"]

PRECIP = {(s, d): 2.0 + 0.05 * d + 0.3 * i for i, s in enumerate(STATIONS) for d in range(N_DAYS)}
TEMP = {(s, d): 10.0 + 0.1 * d - 0.5 * i for i, s in enumerate(STATIONS) for d in range(N_DAYS)}
EVAP = {(s, d): 3.0 + 0.02 * d for s in STATIONS for d in range(N_DAYS)}
ETA_OBS = {(s, d): 2.5 + 0.02 * d for s in STATIONS for d in range(N_DAYS)}
Q_OBS = {(s, d): 5.0 + 0.1 * d + i for i, s in enumerate(STATIONS) for d in range(N_DAYS)}

BOUNDS_16ROW = [
    "200.0", "400.0", "1.5", "3.0", "0.3", "0.7", "2.0", "4.0",
    "0.05", "0.2", "0.01", "0.03", "0.2", "1.0", "0.0", "0.3",
]


def _write_csv(path, values):
    with open(path, "w") as f:
        for d in range(N_DAYS):
            f.write(",".join("%.6f" % values[(s, d)] for s in STATIONS) + "\n")


class TestHRUIntegration(TestCase):
    dem = "hbvhrui_dem"
    outlets = "hbvhrui_outlets"
    basins = "hbvhrui_basins"
    basins_vector = "hbvhrui_basins_v"

    @classmethod
    def setUpClass(cls):
        cls.datadir = os.path.join(TESTDIR, "_hru_integration_data")
        os.makedirs(cls.datadir, exist_ok=True)

        with open(os.path.join(cls.datadir, "bounds16.csv"), "w") as f:
            for row in BOUNDS_16ROW:
                f.write(row + "\n")

        _write_csv(os.path.join(cls.datadir, "precip.csv"), PRECIP)
        _write_csv(os.path.join(cls.datadir, "temp.csv"), TEMP)
        _write_csv(os.path.join(cls.datadir, "evap.csv"), EVAP)
        _write_csv(os.path.join(cls.datadir, "etobs.csv"), ETA_OBS)
        _write_csv(os.path.join(cls.datadir, "dischargeobs.csv"), Q_OBS)

        cls.use_temp_region()
        cls.runModule("g.region", n=1.0, s=0.0, e=1.0, w=0.0, res=0.002)
        cls.runModule(
            "r.mapcalc",
            expression="%s = 100 + 0.05*col() + 0.03*row()"
            " + 10*sin(col()/8.0) + 8*cos(row()/6.0)" % cls.dem,
            overwrite=True,
        )
        call_module(
            "v.in.ascii",
            input="-",
            stdin="0.95,0.02,BasinA\n0.02,0.95,BasinB\n",
            output=cls.outlets,
            separator="comma",
            columns="x double precision, y double precision, id varchar(20)",
            overwrite=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()
        cls.runModule(
            "g.remove",
            flags="f",
            type=["raster", "vector"],
            name=[cls.dem, cls.outlets, cls.basins, cls.basins_vector],
        )
        shutil.rmtree(cls.datadir, ignore_errors=True)

    def setUp(self):
        self.outdir = os.path.join(TESTDIR, "_out_hru_integration")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_hru_bands_end_to_end(self):
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation=os.path.join(self.datadir, "precip.csv"),
            temperature=os.path.join(self.datadir, "temp.csv"),
            evapotranspiration=os.path.join(self.datadir, "evap.csv"),
            eta_observed=os.path.join(self.datadir, "etobs.csv"),
            discharge_observed=os.path.join(self.datadir, "dischargeobs.csv"),
            elevation=self.dem,
            outlets=self.outlets,
            id_column="id",
            threshold=20,
            snap_radius=15,
            hru_bands=2,
            parameters_template=os.path.join(self.datadir, "bounds16.csv"),
            basins=self.basins,
            basins_vector=self.basins_vector,
            n_calib_steps=15,
            n_days=N_DAYS,
            n_years=1,
            warmup=5,
            n_realizations=10,
            output=self.outdir,
            overwrite=True,
        )

        # 2 basins x 2 bands = 4 HRUs -> 4 Output/EToutput file groups,
        # named after each HRU (<basin_id>-b<band_index>), one per
        # basin's own first-band representative -- but named by the
        # basin_id r.hydro.hbv.basins assigned each HRU row, which for
        # hru_bands mode is the *parent* basin_id repeated per band;
        # r.hydro.hbv itself only reports the first rtot (=2) rows'
        # worth of Output/EToutput (see hbv_performance.c's own
        # btot-1-vs-rtot-1 split, documented in main.c/docs), so exactly
        # 2 Output/EToutput groups exist, matching the 2 Basinout/ETout
        # (station-level) groups.
        found_basinout = sorted(
            f for f in os.listdir(self.outdir) if f.startswith("Basinout")
        )
        found_etout = sorted(
            f for f in os.listdir(self.outdir) if f.startswith("ETout")
            and not f.startswith("EToutput")
        )
        self.assertEqual(len(found_basinout), 2)
        self.assertEqual(len(found_etout), 2)

        for fname in found_basinout:
            path = os.path.join(self.outdir, fname)
            with open(path) as f:
                rows = [line.split() for line in f if line.strip()]
            self.assertEqual(len(rows), 15, msg="%s: row count" % fname)
            for qo, qr, po in rows:
                qr = float(qr)
                self.assertTrue(
                    __import__("math").isfinite(qr) and qr >= 0.0,
                    msg="%s: qr must be finite and >=0, got %s" % (fname, qr),
                )


if __name__ == "__main__":
    test()
