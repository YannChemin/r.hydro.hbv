"""Tests r.hydro.hbv's GRASS-native basin-delineation interface: instead
of a "parameters" CSV + "basin_ids" text file, r.hydro.hbv can spawn
r.hydro.hbv.basins itself (elevation+outlets given) to delineate basins
and write physiography/parameter bounds straight into a basins vector's
attribute table, or read an already-delineated basins vector directly
(basins_vector given alone) -- no CSV round-trip for basin setup, only
the meteorological time series stay as CSV (raster forcing is out of
scope for this addon, see docs/raster_options.md).

Uses a small synthetic DEM/outlets fixture (not a real catchment -- see
testsuite/test_basins.py in the r.hydro.hbv.basins module for the same
caveat) paired with the *real* original-mode time series fixture
(testsuite/data/original/*.csv) so the model actually has meaningful
forcing/observed data to run against.
"""

import csv
import io
import os

from grass.gunittest.case import TestCase
from grass.gunittest.gmodules import call_module
from grass.gunittest.main import test

TESTDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(TESTDIR, "data", "original")

BOUNDS_16ROW = [
    "1.0", "1000.0",  # fc lo/hi
    "0.5", "5.0",  # beta lo/hi
    "0.2", "0.9",  # lp lo/hi
    "1.0", "5.0",  # alpha lo/hi
    "0.01", "0.3",  # kf lo/hi
    "0.005", "0.05",  # ks lo/hi
    "0.1", "2.0",  # perc lo/hi
    "0.0", "0.5",  # cflux lo/hi
]


class TestBasinsIntegration(TestCase):
    dem = "hbvint_dem"
    landcover = "hbvint_landcover"
    outlets = "hbvint_outlets"
    basins = "hbvint_basins"
    basins_vector = "hbvint_basins_v"

    @classmethod
    def setUpClass(cls):
        cls.use_temp_region()
        cls.runModule("g.region", n=1.0, s=0.0, e=1.0, w=0.0, res=0.002)
        cls.runModule(
            "r.mapcalc",
            expression="%s = 100 + 0.05*col() + 0.03*row()"
            " + 10*sin(col()/8.0) + 8*cos(row()/6.0)" % cls.dem,
            overwrite=True,
        )
        cls.runModule(
            "r.mapcalc",
            expression="%s = if((row()+col()) %% 3 == 0, 1, 2)" % cls.landcover,
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

        cls.bounds_path = os.path.join(TESTDIR, "_bounds16.csv")
        with open(cls.bounds_path, "w") as f:
            for i in range(0, 16, 2):
                f.write("%s\n%s\n" % (BOUNDS_16ROW[i], BOUNDS_16ROW[i + 1]))

    @classmethod
    def tearDownClass(cls):
        cls.del_temp_region()
        cls.runModule(
            "g.remove",
            flags="f",
            type=["raster", "vector"],
            name=[
                cls.dem,
                cls.landcover,
                cls.outlets,
                cls.basins,
                cls.basins_vector,
            ],
        )
        if os.path.exists(cls.bounds_path):
            os.remove(cls.bounds_path)

    def setUp(self):
        self.outdir = os.path.join(TESTDIR, "_out_basins_integration")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.outdir, ignore_errors=True)

    def _time_series_options(self):
        return dict(
            dataset="custom",
            precipitation=os.path.join(DATADIR, "precip.csv"),
            temperature=os.path.join(DATADIR, "temp.csv"),
            evapotranspiration=os.path.join(DATADIR, "evap.csv"),
            eta_observed=os.path.join(DATADIR, "etobs.csv"),
            discharge_observed=os.path.join(DATADIR, "dischargeobs.csv"),
            n_calib_steps=916,
            n_days=1910,
            n_years=5,
            warmup=274,
            n_realizations=20,
        )

    def _assert_basin_outputs(self, outdir, expected_ids):
        for basin_id in expected_ids:
            for prefix in ("Output", "Basinout", "ETout", "EToutput"):
                path = None
                for i in range(1, len(expected_ids) + 1):
                    candidate = os.path.join(
                        outdir, "%s%02d-%s.csv" % (prefix, i, basin_id)
                    )
                    if os.path.exists(candidate):
                        path = candidate
                        break
                self.assertIsNotNone(
                    path,
                    msg="no %s file found for basin %s in %s"
                    % (prefix, basin_id, outdir),
                )
                self.assertGreater(os.path.getsize(path), 0)

    def test_delineate_from_elevation_and_outlets(self):
        """r.hydro.hbv spawns r.hydro.hbv.basins itself when given
        elevation+outlets, and reads the resulting vector directly --
        no parameters/basin_ids CSV files involved."""
        opts = self._time_series_options()
        opts.update(
            elevation=self.dem,
            outlets=self.outlets,
            id_column="id",
            threshold=20,
            snap_radius=15,
            landcover=self.landcover,
            forest_cats=[1],
            parameters_template=self.bounds_path,
            basins=self.basins,
            basins_vector=self.basins_vector,
            output=self.outdir,
        )
        self.assertModule("r.hydro.hbv", **opts, overwrite=True)

        self.assertRasterExists(self.basins)
        self.assertVectorExists(self.basins_vector)
        self._assert_basin_outputs(self.outdir, ["BasinA", "BasinB"])

        rows = call_module(
            "v.db.select",
            map=self.basins_vector,
            columns=["basin_id", "area_km2", "fc_lo", "fc_hi"],
            format="csv",
        )
        reader = csv.reader(io.StringIO(rows.strip()))
        next(reader)
        data = list(reader)
        self.assertEqual(len(data), 2)
        for _, area, fc_lo, fc_hi in data:
            self.assertGreater(float(area), 0.0)
            self.assertEqual(float(fc_lo), 1.0)
            self.assertEqual(float(fc_hi), 1000.0)

    def test_reuse_existing_basins_vector(self):
        """r.hydro.hbv can also just read an already-delineated basins
        vector (basins_vector alone, no elevation/outlets) -- e.g. one
        built by a separate r.hydro.hbv.basins run."""
        self.assertModule(
            "r.hydro.hbv.basins",
            elevation=self.dem,
            outlets=self.outlets,
            id_column="id",
            threshold=20,
            snap_radius=15,
            landcover=self.landcover,
            forest_cats=[1],
            parameters_template=self.bounds_path,
            basins=self.basins,
            basins_vector=self.basins_vector,
            overwrite=True,
        )

        opts = self._time_series_options()
        opts.update(basins_vector=self.basins_vector, output=self.outdir)
        self.assertModule("r.hydro.hbv", **opts, overwrite=True)

        self._assert_basin_outputs(self.outdir, ["BasinA", "BasinB"])

    def test_elevation_without_outlets_fails(self):
        opts = self._time_series_options()
        opts.update(elevation=self.dem, output=self.outdir)
        self.assertModuleFail("r.hydro.hbv", **opts, overwrite=True)

    def test_elevation_outlets_without_basins_vector_fails(self):
        opts = self._time_series_options()
        opts.update(
            elevation=self.dem,
            outlets=self.outlets,
            parameters_template=self.bounds_path,
            basins=self.basins,
            output=self.outdir,
        )
        self.assertModuleFail("r.hydro.hbv", **opts, overwrite=True)

    def test_elevation_outlets_without_parameters_template_fails(self):
        opts = self._time_series_options()
        opts.update(
            elevation=self.dem,
            outlets=self.outlets,
            basins=self.basins,
            basins_vector=self.basins_vector,
            output=self.outdir,
        )
        self.assertModuleFail("r.hydro.hbv", **opts, overwrite=True)


if __name__ == "__main__":
    test()
