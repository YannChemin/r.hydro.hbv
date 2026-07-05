# r.hydro.hbv

A [GRASS GIS](https://grass.osgeo.org/) addon implementing the HBV
(Hydrologiska Byråns Vattenbalansavdelning) rainfall-runoff model with
Monte-Carlo/GLUE calibration, OpenMP-parallelized. Ported from an
existing standalone OpenMP HBV implementation into a full GRASS
workflow: real terrain and climate data in, calibrated discharge/ETa
simulations out, with no CSV file prep required if you don't want one.

```
g.region n=34 s=31 e=49 w=47 res=0:06

r.hydro.hbv dataset=custom \
  dem_source=copernicus_glo30 outlets=outlets id_column=id \
  threshold=1000 snap_radius=30 parameters_template=bounds16.csv \
  basins=basins basins_vector=basins_v \
  era5_start=2001-01-01 era5_end=2001-03-31 \
  eta_observed=etobs.csv discharge_observed=dischargeobs.csv \
  n_calib_steps=90 n_years=1 warmup=10 n_realizations=45000 \
  output=results
```

The one-liner above delineates real sub-basins from a DEM fetched on
the fly, drives the model with real ERA5(-Land) precipitation/
temperature/PET, and calibrates against observed discharge/ETa — for
any catchment on Earth with only outlet coordinates and parameter
bounds as local input.

## Modes

- **`dataset=original`** — the bundled Karkheh basin (Iran) dataset
  from Muthuwatta's PhD thesis (ITC, 2005): 8 sub-basins, CSV
  parameters/forcing/observations.
- **`dataset=dicrim`** — the bundled Plumergat flood-risk dataset
  (Brittany, France).
- **`dataset=custom`** — bring your own inputs, CSV or GRASS-native
  (below), lumped or semi-distributed.

## GRASS-native basin delineation (instead of CSV)

Give **elevation** + **outlets** (a vector point per station) instead
of a "parameters" CSV and "basin_ids" text file, and `r.hydro.hbv`
spawns *r.hydro.hbv.basins* itself to delineate sub-basins
(`r.watershed`/`r.stream.snap`/`r.stream.basins`), compute physiography,
and write it all directly into a basins vector's attribute table.
Already have a delineated basins vector from a previous run? Pass
**basins_vector** alone and skip delineation entirely.

No DEM handy either? **dem_source=copernicus_glo30** (or `_glo90`)
fetches one via [r.in.dem](https://github.com/YannChemin/r.in.dem) —
a real, no-API-key-needed global DEM — for the current region before
delineating.

## GRASS-native time-series I/O (instead of CSV)

Precipitation/temperature/evapotranspiration/observed-ETa/observed-
discharge each accept a `*_table=` alternative: a long-format
(`station_id, date, value`) GRASS DB table instead of a CSV file — e.g.
one built by *r.hydro.hbv.forcing* from a climate STRDS (ERA5-Land,
MODIS, ...) via zonal means.

Don't have that STRDS yet? **era5_start=**/**era5_end=** replace
precipitation/temperature/evapotranspiration entirely: `r.hydro.hbv`
spawns [t.in.era5](https://github.com/YannChemin/t.in.era5) itself
(ERA5-Land, falling back to ERA5) and *r.hydro.hbv.forcing* to build
those tables automatically. Observed discharge/ETa are real
measurements, not reanalysis output, so they always still need a CSV
or a pre-built table.

`output_tables=<prefix>` additionally writes the discharge/ETa results
as GRASS DB tables, alongside the CSV outputs.

## Semi-distributed (HRU) mode

`hru_bands=N` (default 1, today's lumped behavior) splits each basin
into N elevation-band hydrological response units, each with
independent state and its own Monte-Carlo parameter draw — the
rainfall-runoff engine itself (`hbv_model.c`/`hbv_performance.c`) is
untouched; only forcing distribution and result aggregation change.

## Testing

```
testsuite/test_original.py               # bundled Karkheh dataset
testsuite/test_dicrim.py                 # bundled Plumergat dataset
testsuite/test_basins_integration.py     # elevation+outlets delineation path
testsuite/test_table_io.py               # GRASS DB table forcing/output
testsuite/test_hru.py                    # semi-distributed HRU mode
testsuite/test_karkheh_era5_v2.py        # real ERA5 forcing, real network
```

The last one needs live network access and a working `~/.cdsapirc` (see
*t.in.era5*'s README), so it's gated behind an environment variable:

```
R_HYDRO_HBV_RUN_ERA5_TESTS=1 python3 -m grass.gunittest.main testsuite/test_karkheh_era5_v2.py
```

The other five run offline against synthetic terrain/data and pass by
default.

## Requirements

- GRASS GIS core (`r.watershed`, `r.stream.snap`, `r.stream.basins`,
  `r.to.vect`, the temporal framework)
- OpenMP-capable C compiler
- [r.hydro.hbv.basins](https://github.com/YannChemin/r.hydro.hbv.basins)
  and [r.hydro.hbv.forcing](https://github.com/YannChemin/r.hydro.hbv.forcing)
  for delineation/forcing-table features -- see
  [docs/raster_options.md](docs/raster_options.md) for their design
- [t.in.era5](https://github.com/YannChemin/t.in.era5) and
  [r.in.dem](https://github.com/YannChemin/r.in.dem) for the
  fetch-forcing/fetch-DEM-directly options

## Install

```
g.extension extension=r.hydro.hbv url=https://github.com/YannChemin/r.hydro.hbv
```

## License

Public domain — see [LICENSE](LICENSE) (Unlicense).

## Background

This addon builds on an existing standalone (non-GRASS) OpenMP HBV
implementation — see [HBV](https://github.com/YannChemin/HBV) for that
original codebase, its own history/README, and the Karkheh/DICRIM
reference materials (`muthuwatta.pdf`).
