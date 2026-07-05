/* Builds the (btot x 22) "parameters" table and the per-row ID array
 * directly from a GRASS basins vector map (as produced by
 * r.hydro.hbv.basins), instead of reading a "parameters" CSV and a
 * "basin_ids" text file from disk. If elevation/outlets are given,
 * r.hydro.hbv.basins is spawned first to (re)build basins/basins_vector;
 * otherwise basins_vector is assumed to already exist (e.g. from a
 * previous r.hydro.hbv.basins run) and is read as-is.
 *
 * *btot_out rows are computation units: one per basin normally, or one
 * per elevation-band HRU when the vector was built with hru_bands > 1
 * (several rows can then share the same basin_id). *rtot_out/
 * *station_ids_out are the distinct basin_id count/values (reporting
 * stations); *station_of_out[b] gives each row's 0-based index into
 * station_ids_out, in row order matching *f7p_out and *basin_ids_out. For
 * hru_bands == 1 (today's lumped case), rtot_out == btot_out and
 * station_of_out[b] == b for every b. */
void hbv_basins_from_vector(
    const char *elevation, const char *outlets, const char *id_column,
    const char *threshold, const char *snap_radius, const char *landcover,
    const char *forest_cats, const char *ffo_default,
    const char *ffi_default, const char *precip_station_elevation,
    const char *temp_station_elevation, const char *et_station_elevation,
    const char *parameters_template, const char *hru_bands,
    const char *basins_raster, const char *basins_vector, int *btot_out,
    float ***f7p_out, char ***basin_ids_out, int *rtot_out,
    int **station_of_out, char ***station_ids_out);
