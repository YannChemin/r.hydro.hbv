# Toward raster-native HBV: options and a recommended path

## Status update

Since this document was first written, both **Option A** (GRASS-native
basin delineation, below) and **Option B** (semi-distributed elevation-
band HRUs, below) have been implemented, plus GRASS-native time-series
I/O beyond delineation:

- `r.hydro.hbv.basins` delineates basins via `r.watershed` +
  `r.stream.snap` + `r.stream.basins` (Option A) and, with
  `hru_bands > 1`, further splits each basin into elevation-band HRUs
  (Option B), writing physiography/parameter bounds directly into the
  basins vector's attribute table (`-v`-flagged `r.to.vect`, so vector
  `cat` matches the raster category the zonal stats are keyed by --
  see "A correctness fix" below).
- `r.hydro.hbv` reads that vector directly (`basins_input.c`), or a
  `basins_vector=` built by a previous run, instead of a "parameters"
  CSV + "basin_ids" text file -- `n_basins` isn't needed either way.
- `r.hydro.hbv.forcing` (new module) reduces a climate STRDS (ERA5-Land,
  MODIS, or any other daily raster series) to a per-basin long-format
  DB table via one `t.rast.univar zones=<basins>` call, usable directly
  as `r.hydro.hbv`'s `precipitation_table=`/`temperature_table=`/etc.
  (an alternative to the CSV options, not a replacement).
- `r.hydro.hbv` can also write its own discharge/ETa time-series results
  as GRASS DB tables (`output_tables=<prefix>`), in addition to the
  existing CSVs.
- A `t.in.era5` module (ERA5/ERA5-Land -> STRDS import) is noted as a
  natural next addition but not yet built -- see that module's own
  section below.

For HRUs specifically: `hbv_model.c` and `hbv_performance.c` were left
**completely untouched**, as originally planned. Each HRU gets its own
Monte-Carlo parameter draw and state (reusing `hbv_performance.c`'s
own pre-existing, previously-unexercised `rtot`-vs-`btot` distinction --
see `main.c`'s station-aggregation comments), and reuses its parent
basin's forcing time series corrected by its own `dep`/`det`/`dee`
(`hbv_model.c` already applied that correction per index, so no engine
change was needed there). One real limitation surfaced by this: because
`hbv_performance.c`'s selection-counter bookkeeping (`m[]`/`mtot[]`) is
only ever advanced per **station**, not per HRU, the `Output`/
`EToutput` parameter-diagnostic files only report each station's
first (lowest-elevation-band) HRU as a representative draw, not every
HRU individually -- the `Basinout`/`ETout` discharge/ETa time series
(the results that actually matter for calibration) are fully and
correctly aggregated across all of a station's HRUs, though. Fixing the
diagnostic-file limitation would require changes to
`hbv_performance.c` itself, which was out of scope for this pass.

Two real, unrelated bugs were also found and fixed while wiring this
up, both pre-existing (silently harmless only because `rtot` had always
equaled `btot` until HRUs made them differ):
`mfc`/`mbeta`/`mlp`/`malpha`/`mkf`/`mks`/`mperc`/`mcflux` were
allocated at `rtot` size in `main.c` despite `hbv_performance.c` writing
them with a `b<btot-1` loop (an out-of-bounds write once `btot > rtot`);
and the per-realization accumulator reset loop was iterating `b<btot`
over arrays sized `rtot`. Both are now sized/bounded correctly.

### A correctness fix: `r.to.vect` needs `-v`

While building the HRU attribute-writing path, testing on a
deliberately out-of-order synthetic raster (categories 9, 5, 17 in that
scan order) revealed that plain `r.to.vect` assigns vector `cat` in
polygon-discovery order, **not** by raster cell value -- e.g. raster
value 9 could become vector `cat=1`. Since `r.hydro.hbv.basins` keys its
zonal stats (`r.stats`/`r.univar`) by raster category and then `UPDATE
... WHERE cat=<raster_category>`, this could silently swap which
basin's physiography/bounds land on which polygon whenever discovery
order didn't match category order (only working by coincidence
otherwise). Fixed by adding `r.to.vect`'s `-v` flag ("use raster values
as categories instead of unique sequence"), which makes vector `cat`
equal the raster category directly. Caught and fixed before it ever
shipped a wrong result in the testsuites -- worth keeping in mind for
any future raster-to-vector zonal-stats pipeline in this project.

## Where things stand today

`r.hydro.hbv` operates entirely on **lumped, basin-averaged time
series** supplied as CSV: one value per basin per day for
precipitation, temperature, PET, observed discharge and observed ETa,
plus one scalar per basin for area/forest-fraction/elevation-difference
terms. There is no raster I/O in the module at all.

For the Plumergat DICRIM study, the upstream Python pipeline
(`~/Documents/LaTex/DICRIM/Plumergat/scripts/09_hbv_prep.py`) already
does real raster work to arrive at those lumped numbers: it delineates
the 11 SHYREG sub-basin catchments from the RGE ALTI 1m DEM
(EPSG:2154, `~/RSDATA/RGE_D056/.../RGEALTI_2-0_1M_ASC_LAMB93-IGN69_...`)
using **pysheds** (`fill_pits`/`fill_depressions`/`resolve_flats`/
`flowdir`/`catchment`/`polygonize`), falling back to a circular-polygon
approximation sized to match the SHYREG-reported area when pysheds or
the DEM isn't available. None of this uses GRASS.

The question is how much of that lumped-basin averaging to keep, and
how to introduce actual raster grids (DEM-derived catchments, gridded
meteorological forcing, land cover) without turning this into a
ground-up rewrite.

## Option A — GRASS-native basin delineation, engine unchanged

Replace pysheds with `r.watershed` + `r.water.outlet` (or
`r.stream.basins` given known outlet coordinates) run on the same RGE
ALTI 1m DEM, then derive the per-basin scalars the engine already
expects (`area`, `ffo`/`ffi`, elevation-difference terms) via
`r.univar -t` / `v.rast.stats` zonal statistics against the delineated
catchment polygons, instead of the circular-polygon fallback.

- **Effort**: low. No change to `hbv_model.c`/`hbv_performance.c`/
  `main.c` at all -- only the pre-processing step that produces
  `param_sto.csv` changes.
- **Value**: real hydrologic delineation (flow-accumulation based)
  everywhere, not just where pysheds succeeds; one less third-party
  Python dependency; consistent with how a GRASS addon's users would
  expect basin delineation to happen (inside GRASS, via `g.region` +
  standard raster tools) rather than via an external library.
- **Risk**: low. `r.watershed`/`r.stream.basins` on a 1m DEM over even
  a small commune is a well-trodden, inexpensive raster operation.

This is the recommended next step: it's a pure pre-processing swap,
fully decoupled from the Monte-Carlo engine, and immediately reusable
by both `dataset=dicrim` and any future custom basin.

## Option B — semi-distributed HRUs

Split each current lumped sub-basin into Hydrological Response Units
(HRUs) by elevation band (from the DEM, via `r.watershed` sub-basins
crossed with elevation reclassification) and land cover (from BD TOPO
or a GRASS land-cover raster), each HRU keeping its **own** soil
moisture/snowpack/groundwater state, with basin discharge = sum of HRU
discharges (optionally routed with a per-HRU lag before summation).

- **Effort**: moderate. The state arrays in `main.c`
  (`ssm`/`ssw`/`sgw`/`ssp`/`smw`/`sgwx`, all currently `[basin][time]`)
  would become `[HRU][time]`; `hbv_model.c`'s physics is unchanged
  per-HRU (it already only knows about a single "basin" index `b`).
  The parameter/physiography tables (`param_sto.csv`'s `ffo`/`ffi`/
  elevation terms) would need one row per HRU rather than per basin,
  and `hbv_report.c`'s aggregation would need a basin->HRU grouping
  map to sum discharge back up before writing `Basinout*.csv`.
- **Value**: captures elevation/land-cover heterogeneity (snowline,
  forest interception) within a basin without paying for full
  cell-by-cell simulation.
- **Risk**: moderate -- mostly bookkeeping (which HRU belongs to which
  basin), not new physics.

## Option C — fully distributed, cell-based

Run HBV per DEM grid cell: read gridded precipitation/temperature/PET
rasters row-by-row via `Rast_get_row`, write discharge/state rasters
via `Rast_put_row`, route cell discharge downstream along the
flow-direction raster from `r.watershed`.

- **Effort**: high, and it's a genuine algorithmic change, not just a
  larger loop bound. The current lumped code allocates full
  `[basin][ttot+1]` time series for every state variable because
  `btot` is small (8-11); for a raster this becomes `[nrows*ncols][ttot+1]`,
  which is not viable even at moderate resolution over a small commune
  -- a 1m DEM clipped to a few km2 is already millions of cells, times
  thousands of days, times 4 bytes/float, times ~6 state variables.
  Two changes are unavoidable if this is pursued:
  - **Coarsen resolution**: aggregate to 25-50m cells (or coarser
    elevation-band HRU cells, blurring into Option B) rather than
    keeping the native 1m DEM resolution for the hydrologic model
    itself (1m is appropriate for delineation, not for per-cell state
    simulation).
  - **Drop full-time-series state storage**: since only the final
    output time series (discharge, ETa) is needed, not per-run
    intermediate storage, state should be kept as a rolling 2-slot
    buffer (`t`, `t+1`) rather than the current `[cell][ttot+1]`
    allocation -- a real departure from how `main.c` is structured
    today, where every intermediate day of every state variable is
    kept in memory for the whole run (cheap when `btot` is 8-11, not
    when it's a raster's cell count).
  - Flow routing between cells (needed once cells aren't independent
    lumped units) adds an ordering dependency that the current
    `#pragma omp parallel for` over time steps does not have to deal
    with today -- see the note below on an existing correctness gap
    that this would make worse, not better.
- **Value**: the most physically complete answer, but the cost only
  pays off if within-basin spatial variability of forcing/land surface
  actually matters for the DICRIM use case (flood risk at the commune
  outlet) more than better basin delineation (Option A) or HRU
  splitting (Option B) already would.

## An existing correctness note worth fixing regardless of raster work

While verifying `r.hydro.hbv` against the historical Karkheh reference
outputs (see `testsuite/test_original.py`'s module docstring), it was
confirmed -- by re-running the original, unmodified `~/dev/HBV/hbv`
binary twice -- that the Monte-Carlo summary files (`Output*.csv`,
`EToutput*.csv`) are **not** run-to-run reproducible even in the
pristine 2005 code. The likely cause: `main.c`'s inner time-step loop
is parallelized with `#pragma omp parallel for ... for (t=0;t<ttot;t++)`,
but the model's state update at `t+1` depends on state computed at `t`
(`ssm[b][t+1]` depends on `ssm[b][t]`, etc.) -- a sequential dependency
that plain OpenMP `parallel for` over `t` does not respect. The
deterministic (last-realization) `Basinout*.csv`/`ETout*.csv` files
usually reproduce bit-for-bit run-to-run (confirmed repeatedly at both
small and full 45000-realization scale for the original Karkheh
dataset), but not reliably: a later full-scale rerun (after an unrelated
fix elsewhere in `hbv_performance.c`, see below) showed one basin's
`Basinout` differing by a few m3/s in a handful of rows out of 916,
consistent with the same race occasionally manifesting there too. Any of
Options A/B/C -- and Option C especially, since cell-to-cell routing
genuinely does need correct ordering -- should parallelize over the
**basin/HRU/cell** dimension (which has no sequential dependency across
basins) rather than the **time** dimension (which does), or make the
time loop sequential and parallelize the Monte-Carlo realization loop
instead.

A second, unrelated bug was found (and fixed, since it isn't merely a
reproducibility quirk but an out-of-bounds memory access) while running
the full-scale DICRIM dataset through this same code path:
`hbv_performance.c`'s annual discharge-deficit loop computed
`minVal = 1 + (t-1)*TY` for the discharge-deficit window, which is
negative at `t=0` (`1 - 365`), reading `qo`/`qr` at a large negative
column index. This happened to stay within readable (if wrong) heap
memory for the small 8-basin/916-day Karkheh arrays, but segfaulted on
the larger 11-basin/8400-day/45000-realization DICRIM arrays. Fixed by
clamping `minVal` to `>= 0`; this only affects `qdo`/`qdr`/`dqd`
(discharge-deficit diagnostics, not written to any of the report CSVs),
not the discharge/ET time series themselves.

## Noted for later: `r.in.landcover`

`r.hydro.hbv.basins`' `landcover=`/`forest_cats=` options (for the
`ffo`/`ffi` forest/field-fraction physiography terms) currently require
the caller to already have a local land-cover raster -- there's no
equivalent of `r.in.dem`/`t.in.era5` (real, no-account-needed global
data fetched directly for the current region) for land cover yet.

**ESA WorldCover 10m 2021 (v200)** is a real, no-auth candidate,
confirmed reachable the same way the Copernicus DEM bucket was: plain
HTTPS, no signed S3 request, public AWS Open Data bucket at
`esa-worldcover.s3.eu-central-1.amazonaws.com`, 3x3-degree COG tiles
named e.g. `ESA_WorldCover_10m_2021_v200_N33E048_Map.tif` (tile ID =
south-west corner, snapped to multiples of 3 degrees -- a different
grid from the Copernicus DEM's 1-degree tiles, so the tiling logic
needs its own implementation, not a reuse of `r.in.dem`'s). Class 10
("Tree cover") is the natural `forest_cats=` value. Categorical data,
so any resampling must be nearest-neighbour, never bilinear/cubic --
unlike a continuous DEM, averaging land-cover class codes produces
meaningless values.

A throwaway version of this fetch-and-warp logic (same overview-based
fast-warp technique as `r.in.dem`'s default path) was written and
verified working for the `iran_karkheh_fetch_landcover.py` demo script
in `$HOME/grassdata` (4 tiles, ~5s for the whole Karkheh region) but
not turned into a proper standalone addon -- that's the follow-up this
note is for, mirroring `r.in.dem`'s structure (`source=`/`area=`/
`cache_dir=`/native-vs-region-resolution flag) but for WorldCover's
3-degree grid and nearest-neighbour-only resampling.

## Noted for later: harden `hbv_model.c` against `lp=0` (and similar) edge cases

Discovered while building the Iran_Karkheh demo project (`$HOME/grassdata/
run_iran_karkheh.sh` + `iran_karkheh_plots.py`): `dataset=original`'s own
bundled parameter bounds (`data/original/param_sto.csv`) set the `lp`
parameter's lower bound to exactly `0.0` for several sub-basins, and
Monte-Carlo sampling can and does draw `lp = 0`. `hbv_model.c`'s actual-ET
term,

```c
eta[b][t] = MIN(etp[b][t], (etp[b][t]*ssm[b][t] / (lp[b][n]*fc[b][n])));
```

divides by `lp[b][n]*fc[b][n]`, which is `0` whenever `lp = 0` (regardless
of `fc`) -- producing `inf`/`nan` for that timestep. Since `ssm[b][t+1]`
is updated from `eta[b][t]` (`hbv_model.c`'s water-balance recursion),
a single `nan` day propagates forward and can corrupt the rest of that
realization's state, not just one row. Confirmed via `ETout01-Doab.csv`
from a real `dataset=original` run: 343 of 916 days came back `nan` in
the simulated-ETa column, and for several sub-basins (Pole_Chehr,
Doabe_M, Ghor_B, Holilan in one such run) *every* reported day was
`nan`, presumably because the specific realization selected for the
final report happened to have drawn `lp = 0`.

This wasn't fixed in that session (`hbv_model.c` was intentionally kept
untouched throughout the GRASS-integration and raster-native work
described above) -- `iran_karkheh_plots.py`'s `nash_sutcliffe()` was
instead made robust to it (masks `nan` days out of the score rather than
letting them propagate into the whole statistic), which is a reasonable
plotting-side mitigation but not a fix for the underlying engine
behavior. The real fix belongs in `hbv_model.c` itself -- e.g. clamping
`lp` away from exactly `0` at the point of sampling (`main.c`, where
`lp[b][n]` is drawn from its prior bounds), or guarding the division in
`hbv_model.c` directly (e.g. treat `lp*fc <= 0` as `eta = 0` or
`eta = etp`, whichever is the hydrologically correct boundary case --
worth checking against the original HBV-96 formulation, Lindström et al.
1997, rather than guessing). Either fix should be checked against
`test_original.py`/`test_dicrim.py` for behavior changes, since both
bundled datasets' parameter bounds allow `lp = 0` and any observed
`nan`-day counts in their existing reference outputs would change.
