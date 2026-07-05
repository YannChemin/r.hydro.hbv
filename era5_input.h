/* Spawns t.in.era5 once (precipitation + temperature +
 * potential_evaporation, in one call, one STRDS per variable) and then
 * r.hydro.hbv.forcing once per resulting STRDS (zonal mean against
 * basins_raster/basins_vector) to build the long-format
 * (station_id,date,value) tables table_read_pivot() already knows how
 * to read -- an alternative to precipitation_table=/temperature_table=/
 * evapotranspiration_table= that fetches real ERA5(-Land) reanalysis
 * data instead of requiring a pre-built table.
 *
 * eta_observed/discharge_observed are NOT covered by this: those are
 * observed calibration targets (gauge discharge, SEBS/remote-sensing
 * ETa), not meteorological forcing, and ERA5 reanalysis has no
 * equivalent -- they must still come from a CSV or a user-supplied
 * table. */
void hbv_era5_forcing(const char *era5_start, const char *era5_end,
		      const char *era5_area, const char *era5_cache_dir,
		      const char *basins_raster, const char *basins_vector,
		      char **precip_table_out, char **temp_table_out,
		      char **evap_table_out);
