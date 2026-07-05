# r.hydro.hbv

## NAME

**r.hydro.hbv** - HBV lumped-basin rainfall-runoff model with
Monte-Carlo/GLUE calibration, OpenMP parallelized.

## SYNOPSIS

**r.hydro.hbv**\
**r.hydro.hbv --help**\
**r.hydro.hbv** [**-l**] [**dataset**=*string*]
[**precipitation**=*name*] [**temperature**=*name*]
[**evapotranspiration**=*name*] [**eta_observed**=*name*]
[**discharge_observed**=*name*] [**parameters**=*name*]
[**basin_ids**=*name*] **output**=*name* [**n_basins**=*integer*]
[**n_calib_steps**=*integer*] [**n_days**=*integer*]
[**n_years**=*integer*] [**warmup**=*integer*]
[**n_realizations**=*integer*] [physical parameter options] [**--verbose**]
[**--quiet**]

### Flags

**-l**
&nbsp;&nbsp;&nbsp;&nbsp;Legacy mode: do not reset performance accumulators between
Monte-Carlo realizations (bit-for-bit reproduces the original 2005
Karkheh code, including its unreset-statistics quirk).

## DESCRIPTION

*r.hydro.hbv* implements the HBV (Hydrologiska Byrans
Vattenbalansavdelning) conceptual rainfall-runoff model over one or more
lumped sub-basins, run **ntot** times with randomly sampled parameter
sets (Monte-Carlo/GLUE style uncertainty analysis), and reports
Nash-Sutcliffe efficiency, relative volume error and related performance
criteria per realization. The simulation loop over time steps and basins
is parallelized with OpenMP.

This module generalizes what were originally two hand-forked, hard-coded
C programs into a single, parametrized engine:

- the original calibration from Lal Muthuwatta's PhD thesis (ITC, 2005),
  8 sub-basins of the Karkheh basin (Iran);
- a Monte-Carlo flood-risk study (11 SHYREG sub-basins) built for the
  Plumergat (Brittany, France) DICRIM municipal risk-information
  document.

Both are bundled as ready-to-run demonstrations via `dataset=original`
and `dataset=dicrim`. Use `dataset=custom` (the default) to point the
module at your own lumped-basin CSV inputs.

## INPUTS

All time-series inputs are headerless CSV files, one row per day and one
column per basin/station:

- **precipitation**: observed precipitation (mm/d)
- **temperature**: observed mean temperature (deg C)
- **evapotranspiration**: potential evapotranspiration (mm/d)
- **eta_observed**: observed actual evapotranspiration used as a
  secondary calibration target (mm/d)
- **discharge_observed**: observed discharge (m3/s)

**parameters** is a CSV with 22 rows and one column per basin: rows 0-15
are paired low/high sampling bounds for the 8 free HBV parameters (fc,
beta, lp, alpha, kf, ks, perc, cflux); row 16 is basin area (km2); rows
17-18 are forest/field area fractions; rows 19-21 are elevation
differences (per 100 m) between basin and station for precipitation,
temperature and evapotranspiration.

**basin_ids** is an optional text file, one identifier per line, used to
name the output files (`Basin01`, `Basin02`, ... if omitted).

### GRASS-native basin delineation (instead of `parameters`/`basin_ids`)

In place of the "parameters" CSV and "basin_ids" text file, basin
physiography and identifiers can come straight from a GRASS vector map
instead of any file on disk:

- **elevation** + **outlets** (+ **id_column**, **threshold**,
  **snap_radius**, **landcover**, **forest_cats**, **basins_ffo**,
  **basins_ffi**, **precip_station_elevation**,
  **temp_station_elevation**, **et_station_elevation**,
  **parameters_template**, **hru_bands**, **basins**, **basins_vector**):
  r.hydro.hbv spawns *[r.hydro.hbv.basins](r.hydro.hbv.basins.md)* itself
  to delineate the basins raster/vector and write physiography plus (via
  **parameters_template**, still a small CSV of just the 8 HBV
  parameters' sampling bounds -- these can't be derived from a DEM)
  the parameter bounds directly into **basins_vector**'s attribute
  table, then reads that table as this run's parameters input.
  **hru_bands** (default 1) splits each basin into that many
  elevation-band HRUs, each with independent state and its own
  Monte-Carlo parameter draw -- see "SEMI-DISTRIBUTED (HRU) MODE" below.
- **basins_vector** alone (no elevation/outlets): reads an
  already-delineated basins vector -- e.g. one built by a previous
  *r.hydro.hbv.basins* or *r.hydro.hbv* run -- directly, skipping
  delineation.

Either way, `n_basins` is not needed: the basin count comes from the
number of populated rows in the vector's attribute table.

**dem_source** (+ optional **dem_area**) replaces **elevation**
entirely: given together with **outlets**, r.hydro.hbv spawns
*[r.in.dem](r.in.dem.md)* itself to fetch a real DEM (Copernicus
GLO-30/GLO-90, no account/API key needed, truly global land coverage)
for the current region before delineating -- so any study area on
Earth can be run against real terrain with no local elevation data
prepared in advance.

### GRASS-native time-series I/O (instead of CSV)

Each of **precipitation**/**temperature**/**evapotranspiration**/
**eta_observed**/**discharge_observed** has a `*_table=` alternative
(`precipitation_table=`, etc.; column names overridable via
`table_id_column=`/`table_date_column=`/`table_value_column=`, default
`station_id`/`date`/`value`): a long-format DB table, e.g. one built by
*[r.hydro.hbv.forcing](r.hydro.hbv.forcing.md)* from an ERA5-Land/MODIS
STRDS via zonal means, instead of a CSV file. `n_days` isn't needed in
that case either -- the day count comes from the tables' own date axis
(every `*_table` input must resolve to exactly the same set of dates).

**output_tables**=*prefix* additionally writes the discharge/ETa results
as `<prefix>_basinout`/`<prefix>_etout` GRASS DB tables (long format:
`station_id,day,qo,qr,po` / `station_id,day,eta,eto`), alongside (not
instead of) the CSV outputs below.

### Fetching precipitation/temperature/PET directly from ERA5

**era5_start**/**era5_end** (`YYYY-MM-DD`, + optional **era5_area**/
**era5_cache_dir**) replace **precipitation**/**temperature**/
**evapotranspiration** (and their `*_table` equivalents) entirely: given
together, r.hydro.hbv spawns
*[t.in.era5](t.in.era5.md)* once (fetching ERA5-Land, falling back to
ERA5, for all three variables) and
*[r.hydro.hbv.forcing](r.hydro.hbv.forcing.md)* three times (zonal
means against **basins**/**basins_vector**, both required in this
mode) to build the same long-format tables the `*_table=` options
expect -- no manual `t.in.era5`/`r.hydro.hbv.forcing` calls needed.
**eta_observed**/**discharge_observed** are calibration targets (real
gauge/remote-sensing observations), not meteorological forcing, so
ERA5 has no equivalent for them -- they must still be given as CSV or
`*_table=`, exactly as in every other mode.

## SEMI-DISTRIBUTED (HRU) MODE

With `hru_bands > 1`, each basin becomes several elevation-band HRUs
(computation units, `btot` in the code), each independently simulated,
but discharge is summed and ETa area-weighted-averaged back to the
**reporting station** level (`rtot`, one per basin, matching observed
data) before scoring and before writing `Basinout`/`ETout`. The
`Output`/`EToutput` parameter-diagnostic files, however, only report
each station's first (lowest-elevation-band) HRU as a representative
draw, not every HRU individually -- a limitation of
`hbv_performance.c`'s own per-station (not per-HRU) selection-counter
bookkeeping, left unchanged; see `docs/raster_options.md` for detail.
`hru_bands=1` (the default) is an exact no-op relative to the lumped,
one-unit-per-basin behavior this module always had.

## OUTPUTS

For every basin, four CSV files are written under **output**:

- `Output<NN>-<id>.csv`: selected ("behavioural") HBV parameter sets
  and their NS/NSH/NSL/RVE/RMAEH/RMAEL/REVE performance scores
- `Basinout<NN>-<id>.csv`: observed discharge, simulated discharge and
  corrected precipitation time series (space-separated)
- `ETout<NN>-<id>.csv`: simulated vs. observed actual evapotranspiration
  time series
- `EToutput<NN>-<id>.csv`: selected parameter sets with
  NSET/RVEET (evapotranspiration-based) scores

## EXAMPLES

Run the bundled original Karkheh thesis dataset:

```sh
r.hydro.hbv dataset=original output=/tmp/hbv_original -l
```

(the `-l` legacy flag reproduces the original 2005 code's behavior
exactly, including not resetting statistics accumulators between
realizations)

Run the bundled Plumergat DICRIM dataset:

```sh
r.hydro.hbv dataset=dicrim output=/tmp/hbv_dicrim
```

Run against custom lumped-basin data:

```sh
r.hydro.hbv dataset=custom \
  precipitation=precip.csv temperature=temp.csv \
  evapotranspiration=evap.csv eta_observed=etobs.csv \
  discharge_observed=dischargeobs.csv parameters=param_sto.csv \
  n_basins=4 n_calib_steps=1000 n_days=1000 n_years=3 warmup=100 \
  n_realizations=1000 output=/tmp/hbv_custom
```

Run against custom data with GRASS-native basin delineation instead of
a "parameters" CSV (r.hydro.hbv spawns *r.hydro.hbv.basins* itself):

```sh
g.region raster=dem
r.hydro.hbv dataset=custom \
  precipitation=precip.csv temperature=temp.csv \
  evapotranspiration=evap.csv eta_observed=etobs.csv \
  discharge_observed=dischargeobs.csv \
  elevation=dem outlets=outlets id_column=ID threshold=1000 \
  landcover=landcover forest_cats=1,2,3 \
  parameters_template=hbv_param_bounds_16row.csv \
  basins=basins basins_vector=basins_v \
  n_calib_steps=1000 n_days=1000 n_years=3 warmup=100 \
  n_realizations=1000 output=/tmp/hbv_custom
```

## NOTES ON EXTENDING TO RASTER INPUTS

The current engine operates entirely on lumped, basin-averaged time
series (no raster I/O). See `docs/raster_options.md` in the module
source tree for a discussion of options to evolve toward
semi-distributed or fully distributed (raster-cell-based) HBV, and a
recommended incremental path.

## REFERENCES

- Muthuwatta, L., 2005. *Calibration of a semi distributed hydrological
  model using discharge and remote sensing data.* PhD Thesis, ITC, The
  Netherlands.
- HBV hydrology model: <https://en.wikipedia.org/wiki/HBV_hydrology_model>

## SEE ALSO

*[r.watershed](r.watershed.md)*, *[r.water.outlet](r.water.outlet.md)*

## AUTHOR

Yann Chemin
