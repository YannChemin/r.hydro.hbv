"""Long-running, network-dependent "version 2" of the original Karkheh
testsuite: instead of the bundled precip/temp/evap CSVs (dataset=original,
see test_original.py), this fetches real ERA5(-Land) reanalysis data for
the 8 named Karkheh sub-basins (Doab, Pole Chehr, Doabe Merek, Ghor
Baghestan, Holilan, Pole Dokhtar, Jelogir, Paye Pol -- names and rough
bounding region confirmed against muthuwatta.pdf / the thesis this
dataset is from) via t.in.era5 + r.hydro.hbv.forcing +
r.hydro.hbv's era5_start=/era5_end= integration, end to end.

Two things this deliberately does NOT reproduce:

1. Real sub-basin geometry. No DEM covering the Karkheh basin (Iran) is
   available in this environment, so basin "delineation" here uses a
   smooth synthetic elevation surface over the real Karkheh bounding
   region purely so r.hydro.hbv.basins has *something* to trace flow
   paths on -- it is not a hydrologically meaningful watershed boundary.
   Likewise the 8 stations' coordinates (KARKHEH_STATIONS below) are
   this session's best-effort placement from general geography (Pol-e
   Dokhtar's is a real, reasonably well-known town; the rest are
   approximate), not authoritative gauge coordinates -- there is no
   station coordinate table in this repository or the thesis PDF's
   extracted text to draw from.
2. eta_observed/discharge_observed. ERA5 is a meteorological reanalysis,
   not an observation network -- it has no equivalent of the gauged
   discharge or SEBS-derived ETa the original calibration used. This
   test reuses the *existing* dataset=original discharge/ETa CSVs
   (testsuite/data/original/*.csv) as calibration targets, exactly as
   test_original.py does, changing only precipitation/temperature/
   potential_evaporation to come from live ERA5(-Land) instead of the
   bundled precip.csv/temp.csv/evap.csv.

What this *does* validate for real: the full t.in.era5 -> zonal-mean ->
long-table -> r.hydro.hbv pipeline against actual network calls to the
Copernicus Climate Data Store, over a real multi-month period within the
thesis's own stated calibration window (April 2000 - March 2002), for
all 8 basins in one run, with month-chunked fetching/caching exercised
across several calendar months (not just the single-month smoke test
covered manually during development).

Requires network access and a working ~/.cdsapirc (see t.in.era5's own
docs for how to set that up) -- skipped unless
R_HYDRO_HBV_RUN_ERA5_TESTS=1 is set, since (unlike the other slow tests
in this suite) it also depends on an external service, not just wall-clock
time.
"""

import math
import os
import shutil
import unittest

from grass.gunittest.case import TestCase
from grass.gunittest.gmodules import call_module
from grass.gunittest.main import test

TESTDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(TESTDIR, "data", "original")

# Best-effort approximate coordinates (WGS84) for the 8 named Karkheh
# sub-basins -- see module docstring point 1. Pol-e Dokhtar is a real,
# identifiable town in Lorestan province; the others are approximate
# placements within the Upper Karkheh region described in the thesis
# (roughly 31-34.5N, 46.5-49E), spread out in plausible upstream-to-
# downstream order, not surveyed gauge locations.
KARKHEH_STATIONS = [
    ("Doab", 46.85, 34.05),
    ("Pole_Chehr", 47.20, 34.15),
    ("Doabe_M", 46.95, 33.95),
    ("Ghor_B", 47.40, 34.25),
    ("Holilan", 47.65, 33.65),
    ("Pole_D", 47.72, 33.17),  # Pol-e Dokhtar, Lorestan -- real town
    ("Jelogir", 48.05, 33.35),
    ("Paye_P", 48.15, 33.80),
]

# Within the thesis's own stated calibration window (April 2000 -
# March 2002); ~90 days is long enough to exercise month-chunked
# fetching/caching across 4 calendar months for 3 variables x 8
# (synthetic) basins, while remaining practical to actually run.
ERA5_START = "2001-01-01"
ERA5_END = "2001-03-31"
N_DAYS = 90


class TestKarkhehEra5V2(TestCase):
    dem = "karkheh_v2_dem"
    outlets = "karkheh_v2_outlets"
    basins = "karkheh_v2_basins"
    basins_vector = "karkheh_v2_basins_v"

    @classmethod
    def setUpClass(cls):
        if not os.environ.get("R_HYDRO_HBV_RUN_ERA5_TESTS"):
            raise unittest.SkipTest(
                "set R_HYDRO_HBV_RUN_ERA5_TESTS=1 to run the live-ERA5 "
                "Karkheh v2 testsuite (needs network + ~/.cdsapirc)"
            )

        cls.era5_cache_dir = os.environ.get(
            "R_HYDRO_HBV_ERA5_CACHE", os.path.join(TESTDIR, "_era5_cache")
        )
        os.makedirs(cls.era5_cache_dir, exist_ok=True)

        cls.bounds_path = os.path.join(TESTDIR, "_karkheh_v2_bounds16.csv")
        # reuse the real original-mode parameter bounds (16 rows, 8
        # columns -- one per Karkheh sub-basin, same order as
        # KARKHEH_STATIONS/basin_ids.txt) rather than inventing new ones
        with open(os.path.join(DATADIR, "param_sto.csv")) as f:
            rows = [line.rstrip("\n") for line in f if line.strip()]
        with open(cls.bounds_path, "w") as f:
            for row in rows[:16]:
                f.write(row + "\n")

        cls.use_temp_region()
        # covers the whole approximate Karkheh region with margin
        cls.runModule("g.region", n=35, s=31, e=49.5, w=46, res=0.02)
        cls.runModule(
            "r.mapcalc",
            expression="%s = 800 + 400*sin(col()/12.0) + 300*cos(row()/9.0)"
            % cls.dem,
            overwrite=True,
        )

        stdin = "\n".join(
            "%f,%f,%s" % (lon, lat, name) for name, lon, lat in KARKHEH_STATIONS
        )
        call_module(
            "v.in.ascii",
            input="-",
            stdin=stdin,
            output=cls.outlets,
            separator="comma",
            columns="x double precision, y double precision, id varchar(20)",
            overwrite=True,
        )

        cls.runModule(
            "r.hydro.hbv.basins",
            elevation=cls.dem,
            outlets=cls.outlets,
            id_column="id",
            threshold=50,
            snap_radius=20,
            parameters_template=cls.bounds_path,
            basins=cls.basins,
            basins_vector=cls.basins_vector,
            overwrite=True,
        )

    @classmethod
    def tearDownClass(cls):
        if not os.environ.get("R_HYDRO_HBV_RUN_ERA5_TESTS"):
            return
        cls.del_temp_region()
        cls.runModule(
            "g.remove",
            flags="f",
            type=["raster", "vector"],
            name=[cls.dem, cls.outlets, cls.basins, cls.basins_vector],
        )
        if os.path.exists(cls.bounds_path):
            os.remove(cls.bounds_path)

    def setUp(self):
        self.outdir = os.path.join(TESTDIR, "_out_karkheh_era5_v2")
        os.makedirs(self.outdir, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_karkheh_era5_v2(self):
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            era5_start=ERA5_START,
            era5_end=ERA5_END,
            era5_cache_dir=self.era5_cache_dir,
            eta_observed=os.path.join(DATADIR, "etobs.csv"),
            discharge_observed=os.path.join(DATADIR, "dischargeobs.csv"),
            basins=self.basins,
            basins_vector=self.basins_vector,
            n_calib_steps=N_DAYS,
            n_years=1,
            warmup=10,
            n_realizations=20,
            output=self.outdir,
            overwrite=True,
        )

        for i, (name, _, _) in enumerate(KARKHEH_STATIONS, start=1):
            for prefix in ("Output", "Basinout", "ETout", "EToutput"):
                path = os.path.join(
                    self.outdir, "%s%02d-%s.csv" % (prefix, i, name)
                )
                self.assertTrue(os.path.isfile(path), msg=path)

        with open(
            os.path.join(self.outdir, "Basinout01-Doab.csv")
        ) as f:
            rows = [line.split() for line in f if line.strip()]
        self.assertEqual(len(rows), N_DAYS)
        for qo, qr, po in rows:
            self.assertTrue(math.isfinite(float(qr)))
            self.assertTrue(math.isfinite(float(po)))


if __name__ == "__main__":
    test()
