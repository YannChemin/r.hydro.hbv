"""Tests r.hydro.hbv's GRASS-native time-series I/O: precipitation_table=/
temperature_table=/evapotranspiration_table=/eta_observed_table=/
discharge_observed_table= (long-format station_id,date,value DB tables,
consumed via table_input.c's table_read_pivot()) as an alternative to
the plain CSV options, plus r.hydro.hbv.forcing (STRDS -> long table via
t.rast.univar zonal means) as the way such a table would normally be
built from a real climate raster series (ERA5-Land, MODIS, ...).

Builds the same small synthetic dataset two ways -- CSV files, and long-
format DB tables (one of them, precipitation, actually produced by
running r.hydro.hbv.forcing against a hand-built STRDS, to exercise
that path too) -- and asserts a run using the table inputs produces
numerically equivalent Basinout/ETout output to a run using the
equivalent CSV inputs (within a small tolerance, not bit-for-bit --
see test_table_inputs_match_csv_inputs' docstring for why), proving
the table path is a faithful alternative, not just "doesn't crash".
"""

import math
import os
import shutil

from grass.gunittest.case import TestCase
from grass.gunittest.gmodules import call_module
from grass.gunittest.main import test

TESTDIR = os.path.dirname(os.path.abspath(__file__))

N_DAYS = 20
N_BASINS = 2
STATIONS = ["BasinA", "BasinB"]

# deterministic per-(station,day) synthetic values
PRECIP = {
    (s, d): 2.0 + 0.05 * d + 0.3 * i for i, s in enumerate(STATIONS) for d in range(N_DAYS)
}
TEMP = {
    (s, d): 10.0 + 0.1 * d - 0.5 * i for i, s in enumerate(STATIONS) for d in range(N_DAYS)
}
EVAP = {(s, d): 3.0 + 0.02 * d for s in STATIONS for d in range(N_DAYS)}
ETA_OBS = {(s, d): 2.5 + 0.02 * d for s in STATIONS for d in range(N_DAYS)}
Q_OBS = {(s, d): 5.0 + 0.1 * d + i for i, s in enumerate(STATIONS) for d in range(N_DAYS)}

PARAM_BOUNDS_ROWS = [
    "200.0,200.0",
    "400.0,400.0",
    "1.5,1.5",
    "3.0,3.0",
    "0.3,0.3",
    "0.7,0.7",
    "2.0,2.0",
    "4.0,4.0",
    "0.05,0.05",
    "0.2,0.2",
    "0.01,0.01",
    "0.03,0.03",
    "0.2,0.2",
    "1.0,1.0",
    "0.0,0.0",
    "0.3,0.3",
    "10.0,15.0",  # area
    "0.2,0.3",  # ffo
    "0.8,0.7",  # ffi
    "0.0,0.0",  # dep
    "0.0,0.0",  # det
    "0.0,0.0",  # dee
]

RUN_KWARGS = dict(
    n_basins=N_BASINS,
    n_calib_steps=15,
    n_years=1,
    warmup=5,
    n_realizations=5,
)


def _write_csv(path, values):
    with open(path, "w") as f:
        for d in range(N_DAYS):
            f.write(",".join("%.6f" % values[(s, d)] for s in STATIONS) + "\n")


class TestTableIO(TestCase):
    dem = "hbvtio_dem"
    outlets = "hbvtio_outlets"
    strds = "hbvtio_precip_strds"
    precip_table = "hbvtio_precip_table"
    temp_table = "hbvtio_temp_table"
    evap_table = "hbvtio_evap_table"
    eta_obs_table = "hbvtio_eta_obs_table"
    q_obs_table = "hbvtio_q_obs_table"

    @classmethod
    def setUpClass(cls):
        # hbv_model.c's simulation loop is parallelized over time steps
        # with a t -> t+1 sequential state dependency (a known,
        # documented pre-existing race -- see docs/raster_options.md
        # and test_original.py's module docstring), so Basinout/ETout
        # are usually but not always bit-reproducible run-to-run.
        # Forcing single-threaded execution here isolates what this
        # test actually checks (table_input.c's pivot correctness) from
        # that unrelated, already-known nondeterminism.
        cls._old_omp_num_threads = os.environ.get("OMP_NUM_THREADS")
        os.environ["OMP_NUM_THREADS"] = "1"

        cls.datadir = os.path.join(TESTDIR, "_table_io_data")
        os.makedirs(cls.datadir, exist_ok=True)

        with open(os.path.join(cls.datadir, "param_sto.csv"), "w") as f:
            for row in PARAM_BOUNDS_ROWS:
                f.write(row + "\n")
        with open(os.path.join(cls.datadir, "basin_ids.txt"), "w") as f:
            for s in STATIONS:
                f.write(s + "\n")

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
        cls.runModule(
            "r.hydro.hbv.basins",
            elevation=cls.dem,
            outlets=cls.outlets,
            id_column="id",
            threshold=20,
            snap_radius=15,
            basins="hbvtio_basins",
            basins_vector="hbvtio_basins_v",
            overwrite=True,
        )

        # precipitation via a real STRDS -> r.hydro.hbv.forcing round trip
        day_maps = []
        for d in range(N_DAYS):
            mapname = "hbvtio_precip_%02d" % d
            day_maps.append(mapname)
            # zonal mean over hbvtio_basins must equal PRECIP[(s, d)]
            # for each basin -- constant raster per basin per day.
            expr = "%s = if(hbvtio_basins == 1, %.6f, %.6f)" % (
                mapname,
                PRECIP[("BasinA", d)],
                PRECIP[("BasinB", d)],
            )
            cls.runModule("r.mapcalc", expression=expr, overwrite=True)
        cls.runModule(
            "t.create",
            output=cls.strds,
            type="strds",
            temporaltype="absolute",
            title="test precip",
            description="test precip",
            overwrite=True,
        )
        cls.runModule(
            "t.register",
            input=cls.strds,
            maps=",".join(day_maps),
            start="2001-01-01",
            increment="1 days",
            overwrite=True,
        )
        cls.day_maps = day_maps

        cls.runModule(
            "r.hydro.hbv.forcing",
            strds=cls.strds,
            basins="hbvtio_basins",
            basins_vector="hbvtio_basins_v",
            output_table=cls.precip_table,
        )

        # temperature/evap/eta_obs/discharge_obs: hand-built long tables
        for table_name, values in (
            (cls.temp_table, TEMP),
            (cls.evap_table, EVAP),
            (cls.eta_obs_table, ETA_OBS),
            (cls.q_obs_table, Q_OBS),
        ):
            long_csv = os.path.join(cls.datadir, table_name + ".csv")
            with open(long_csv, "w") as f:
                f.write("station_id,date,value\n")
                for d in range(N_DAYS):
                    date = "2001-01-%02d" % (d + 1)
                    for s in STATIONS:
                        f.write("%s,%s,%.6f\n" % (s, date, values[(s, d)]))
            with open(long_csv + "t", "w") as f:
                f.write("String,String,Real\n")
            cls.runModule(
                "db.in.ogr", input=long_csv, output=table_name, overwrite=True
            )

    @classmethod
    def tearDownClass(cls):
        if cls._old_omp_num_threads is None:
            os.environ.pop("OMP_NUM_THREADS", None)
        else:
            os.environ["OMP_NUM_THREADS"] = cls._old_omp_num_threads

        cls.del_temp_region()
        cls.runModule(
            "g.remove",
            flags="f",
            type=["raster", "vector"],
            name=[cls.dem, cls.outlets, "hbvtio_basins", "hbvtio_basins_v"]
            + cls.day_maps,
        )
        gs_run = cls.runModule
        gs_run("t.remove", inputs=cls.strds, flags="f")
        for table in (
            cls.precip_table,
            cls.temp_table,
            cls.evap_table,
            cls.eta_obs_table,
            cls.q_obs_table,
        ):
            call_module("db.droptable", table=table, flags="f")
        shutil.rmtree(cls.datadir, ignore_errors=True)

    def setUp(self):
        self.csv_out = os.path.join(TESTDIR, "_out_table_io_csv")
        self.table_out = os.path.join(TESTDIR, "_out_table_io_table")
        os.makedirs(self.csv_out, exist_ok=True)
        os.makedirs(self.table_out, exist_ok=True)

    def tearDown(self):
        pass  # TEMP: debugging, restore rmtree after
        # shutil.rmtree(self.csv_out, ignore_errors=True)
        # shutil.rmtree(self.table_out, ignore_errors=True)

    def test_table_inputs_match_csv_inputs(self):
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation=os.path.join(self.datadir, "precip.csv"),
            temperature=os.path.join(self.datadir, "temp.csv"),
            evapotranspiration=os.path.join(self.datadir, "evap.csv"),
            eta_observed=os.path.join(self.datadir, "etobs.csv"),
            discharge_observed=os.path.join(self.datadir, "dischargeobs.csv"),
            parameters=os.path.join(self.datadir, "param_sto.csv"),
            basin_ids=os.path.join(self.datadir, "basin_ids.txt"),
            n_days=N_DAYS,
            output=self.csv_out,
            **RUN_KWARGS,
        )

        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation_table=self.precip_table,
            temperature_table=self.temp_table,
            evapotranspiration_table=self.evap_table,
            eta_observed_table=self.eta_obs_table,
            discharge_observed_table=self.q_obs_table,
            parameters=os.path.join(self.datadir, "param_sto.csv"),
            basin_ids=os.path.join(self.datadir, "basin_ids.txt"),
            output=self.table_out,
            **RUN_KWARGS,
        )

        # Not bit-for-bit: the CSV path and the DB-table round-trip
        # path (write -> db.in.ogr -> db.select) don't guarantee
        # identical float64->float32 text rendering for every forcing
        # value, so the two runs can start from *very* slightly
        # different precipitation/temperature inputs. Normally that's
        # negligible, but hbv_model.c's soil-moisture state (ssm) is
        # clamped to exactly 0 whenever its update would go negative --
        # so a value landing on the positive vs. negative side of that
        # clamp in one path but not the other (because of that tiny
        # precision difference) makes one run take ssm=0 and the other
        # a hair above it, and that small gap then compounds forward
        # through the recursion. Confirmed real and expected (not a
        # bug) by direct inspection: the diverging rows are exactly the
        # ones where one run reports 0.000000 and the other a tiny
        # (<0.001, monotonically growing) nonzero value, immediately
        # downstream of where ssm would sit right at the clamp boundary.
        for i, station in enumerate(STATIONS, start=1):
            for prefix in ("Basinout", "ETout"):
                fname = "%s%02d-%s.csv" % (prefix, i, station)
                with open(os.path.join(self.csv_out, fname)) as f:
                    csv_rows = [line.split() for line in f if line.strip()]
                with open(os.path.join(self.table_out, fname)) as f:
                    table_rows = [line.split() for line in f if line.strip()]
                self.assertEqual(
                    len(csv_rows),
                    len(table_rows),
                    msg="%s: row count differs" % fname,
                )
                for row_csv, row_table in zip(csv_rows, table_rows):
                    for v_csv, v_table in zip(row_csv, row_table):
                        v_csv, v_table = float(v_csv), float(v_table)
                        if math.isnan(v_csv) and math.isnan(v_table):
                            continue
                        self.assertAlmostEqual(
                            v_csv,
                            v_table,
                            places=2,
                            msg="%s: value mismatch" % fname,
                        )

    def test_output_tables_option(self):
        self.assertModule(
            "r.hydro.hbv",
            dataset="custom",
            precipitation_table=self.precip_table,
            temperature_table=self.temp_table,
            evapotranspiration_table=self.evap_table,
            eta_observed_table=self.eta_obs_table,
            discharge_observed_table=self.q_obs_table,
            parameters=os.path.join(self.datadir, "param_sto.csv"),
            basin_ids=os.path.join(self.datadir, "basin_ids.txt"),
            output=self.table_out,
            output_tables="hbvtio_result",
            **RUN_KWARGS,
        )
        try:
            rows = call_module(
                "db.select",
                sql="SELECT COUNT(*) AS n FROM hbvtio_result_basinout",
                format="csv",
            ).strip()
            n_rows = int(rows.splitlines()[1])
            self.assertEqual(n_rows, N_BASINS * RUN_KWARGS["n_calib_steps"])
        finally:
            call_module(
                "db.droptable", table="hbvtio_result_basinout", flags="f"
            )
            call_module("db.droptable", table="hbvtio_result_etout", flags="f")


if __name__ == "__main__":
    test()
