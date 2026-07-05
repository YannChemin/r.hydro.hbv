/****************************************************************************
 *
 * MODULE:       r.hydro.hbv
 * AUTHOR(S):    Yann Chemin
 * PURPOSE:      HBV (Hydrologiska Byrans Vattenbalansavdelning) rainfall-
 *               runoff model, Monte-Carlo/GLUE calibration, OpenMP
 *               parallelized. Runs either the original Karkheh basin
 *               (Muthuwatta, PhD thesis ITC 2005) dataset, the Plumergat
 *               DICRIM flood-risk dataset, or arbitrary lumped-basin CSV
 *               inputs following the same layout.
 * COPYRIGHT:    (C) 2026 by Yann Chemin
 *               Released into the public domain -- see LICENSE (Unlicense).
 *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <sys/param.h> /*MIN() and MAX()*/
#include <sys/stat.h>

#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/spawn.h>

#include "arrays.h"
#include "readcsv.h"
#include "hbv.h"
#include "basins_input.h"
#include "table_input.h"
#include "era5_input.h"
#include "dem_input.h"

#define PGMNAME "r.hydro.hbv"

/* ------------------------------------------------------------------ */
/* small helpers                                                       */
/* ------------------------------------------------------------------ */

static char *join_path(const char *a, const char *b) {
	size_t n = strlen(a) + strlen(b) + 2;
	char *out = G_malloc(n);
	snprintf(out, n, "%s/%s", a, b);
	return out;
}

/* G_tempfile()'s own name has no extension; OGR's CSV driver needs an
 * actual .csv one to recognize the file, so append one. */
static char *tempfile_with_suffix(const char *suffix) {
	char *base = G_tempfile();
	size_t n = strlen(base) + strlen(suffix) + 1;
	char *out = G_malloc(n);
	snprintf(out, n, "%s%s", base, suffix);
	return out;
}

static char *bundled_data_path(const char *dataset, const char *fname) {
	const char *base = G_gisbase();
	char etcdir[GPATH_MAX];
	snprintf(etcdir, sizeof(etcdir), "%s/etc/%s/data/%s", base, PGMNAME,
		 dataset);
	return join_path(etcdir, fname);
}

/* returns option->answer, or (for dataset != custom) a bundled default
 * path, or NULL if neither is available. */
static char *resolve_path(struct Option *opt, const char *dataset,
			   const char *bundled_fname) {
	if (opt->answer)
		return G_store(opt->answer);
	if (dataset && strcmp(dataset, "custom") != 0)
		return bundled_data_path(dataset, bundled_fname);
	return NULL;
}

static int resolve_int(struct Option *opt, const char *dataset,
			int orig_default, int dicrim_default, int have_default) {
	if (opt->answer)
		return atoi(opt->answer);
	if (dataset && strcmp(dataset, "original") == 0)
		return orig_default;
	if (dataset && strcmp(dataset, "dicrim") == 0)
		return dicrim_default;
	if (have_default)
		return orig_default;
	G_fatal_error(_("Option <%s> requires a value when dataset=custom"),
		      opt->key);
	return 0; /* not reached */
}

static double resolve_double(struct Option *opt) {
	return atof(opt->answer);
}

/* Load one basin ID per line from a text file; falls back to
 * "Basin01".."BasinNN" if path is NULL. Caller owns the returned array. */
static char **load_basin_ids(const char *path, int btot) {
	char **ids = (char **)G_malloc(btot * sizeof(char *));
	int i;

	if (path) {
		FILE *f = fopen(path, "r");
		char buf[512];

		if (!f)
			G_fatal_error(_("Unable to open basin ID file <%s>"), path);
		for (i = 0; i < btot; i++) {
			if (!fgets(buf, sizeof(buf), f))
				G_fatal_error(
				    _("Basin ID file <%s> has fewer than %d lines"),
				    path, btot);
			buf[strcspn(buf, "\r\n")] = '\0';
			ids[i] = G_store(buf);
		}
		fclose(f);
	} else {
		for (i = 0; i < btot; i++) {
			char name[32];
			snprintf(name, sizeof(name), "Basin%02d", i + 1);
			ids[i] = G_store(name);
		}
	}
	return ids;
}

static void ensure_dir(const char *path) {
	if (mkdir(path, 0755) != 0 && errno != EEXIST)
		G_fatal_error(_("Unable to create output directory <%s>"), path);
}

/* Imports a long-format CSV (+ its paired .csvt sidecar, written by
 * hbv_report()) into a GRASS DB table named "<prefix><suffix>", then
 * removes both temp files. */
static void import_long_csv_as_table(const char *csv_path, const char *prefix,
				      const char *suffix) {
	char table_name[GNAME_MAX];
	char csvt_path[GPATH_MAX];
	const char *argv[8];
	int argc = 0;
	char input_arg[GPATH_MAX + 8], output_arg[GNAME_MAX + 8];

	snprintf(table_name, sizeof(table_name), "%s%s", prefix, suffix);

	argv[argc++] = "db.in.ogr";
	snprintf(input_arg, sizeof(input_arg), "input=%s", csv_path);
	argv[argc++] = input_arg;
	snprintf(output_arg, sizeof(output_arg), "output=%s", table_name);
	argv[argc++] = output_arg;
	if (G_get_overwrite())
		argv[argc++] = "--overwrite";
	argv[argc++] = NULL;

	if (G_vspawn_ex(argv[0], argv) != 0)
		G_fatal_error(_("Unable to import <%s> as GRASS DB table <%s>"),
			      csv_path, table_name);

	snprintf(csvt_path, sizeof(csvt_path), "%s.csvt", csv_path);
	remove(csv_path);
	remove(csvt_path);

	G_message(_("Wrote GRASS DB table <%s>"), table_name);
}

/* ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
	struct GModule *module;
	struct Option *opt_dataset, *opt_precip, *opt_temp, *opt_evap,
	    *opt_eta_obs, *opt_q_obs, *opt_params, *opt_basin_ids, *opt_output,
	    *opt_n_basins, *opt_n_calib_steps, *opt_n_days, *opt_n_years,
	    *opt_warmup, *opt_n_realizations, *opt_tcon, *opt_ecevpfo,
	    *opt_ecalt, *opt_pcalt, *opt_tcalt, *opt_sfcf, *opt_foscf,
	    *opt_rfcf, *opt_tt, *opt_tti, *opt_dttm, *opt_cfmax, *opt_focfmax,
	    *opt_cfr, *opt_whc, *opt_sgwmax;
	struct Option *opt_elevation, *opt_outlets, *opt_id_column,
	    *opt_threshold, *opt_snap_radius, *opt_landcover, *opt_forest_cats,
	    *opt_basins_ffo, *opt_basins_ffi, *opt_precip_station_elev,
	    *opt_temp_station_elev, *opt_et_station_elev,
	    *opt_parameters_template, *opt_hru_bands, *opt_basins,
	    *opt_basins_vector;
	struct Option *opt_precip_table, *opt_temp_table, *opt_evap_table,
	    *opt_eta_obs_table, *opt_q_obs_table, *opt_table_id_column,
	    *opt_table_date_column, *opt_table_value_column,
	    *opt_output_tables;
	struct Option *opt_era5_start, *opt_era5_end, *opt_era5_area,
	    *opt_era5_cache_dir;
	struct Option *opt_dem_source, *opt_dem_area;
	struct Flag *flag_legacy;

	int b, n, t;

	G_gisinit(argv[0]);

	module = G_define_module();
	G_add_keyword(_("raster"));
	G_add_keyword(_("hydrology"));
	G_add_keyword(_("model"));
	G_add_keyword(_("HBV"));
	module->description =
	    _("HBV lumped-basin rainfall-runoff model with Monte-Carlo/GLUE "
	      "calibration (OpenMP parallel). Bundles the original Karkheh "
	      "basin thesis dataset and the Plumergat DICRIM flood-risk "
	      "dataset as ready-to-run demonstrations.");

	opt_dataset = G_define_option();
	opt_dataset->key = "dataset";
	opt_dataset->type = TYPE_STRING;
	opt_dataset->required = NO;
	opt_dataset->options = "original,dicrim,custom";
	opt_dataset->answer = "custom";
	opt_dataset->description =
	    _("Bundled demonstration dataset to run (fills in any input "
	      "option left unset); use 'custom' to supply your own inputs");

	opt_precip = G_define_standard_option(G_OPT_F_INPUT);
	opt_precip->key = "precipitation";
	opt_precip->required = NO;
	opt_precip->description =
	    _("CSV of observed precipitation (mm/d), rows=days, cols=basins");

	opt_temp = G_define_standard_option(G_OPT_F_INPUT);
	opt_temp->key = "temperature";
	opt_temp->required = NO;
	opt_temp->description =
	    _("CSV of observed temperature (deg C), rows=days, cols=basins");

	opt_evap = G_define_standard_option(G_OPT_F_INPUT);
	opt_evap->key = "evapotranspiration";
	opt_evap->required = NO;
	opt_evap->description =
	    _("CSV of potential evapotranspiration (mm/d), rows=days, "
	      "cols=basins");

	opt_eta_obs = G_define_standard_option(G_OPT_F_INPUT);
	opt_eta_obs->key = "eta_observed";
	opt_eta_obs->required = NO;
	opt_eta_obs->description =
	    _("CSV of observed actual evapotranspiration used for "
	      "calibration (mm/d), rows=days, cols=basins");

	opt_q_obs = G_define_standard_option(G_OPT_F_INPUT);
	opt_q_obs->key = "discharge_observed";
	opt_q_obs->required = NO;
	opt_q_obs->description =
	    _("CSV of observed discharge (m3/s), rows=days, cols=stations");

	/* GRASS-native alternative to each CSV option above: a long-format
	 * (station_id, date, value) DB table -- e.g. one built by
	 * r.hydro.hbv.forcing from an ERA5/ERA5-Land STRDS, or any other
	 * table following the same schema. Mutually exclusive with the
	 * matching CSV option. */
	opt_precip_table = G_define_option();
	opt_precip_table->key = "precipitation_table";
	opt_precip_table->type = TYPE_STRING;
	opt_precip_table->key_desc = "name";
	opt_precip_table->required = NO;
	opt_precip_table->description =
	    _("Long-format (station_id,date,value) DB table of observed "
	      "precipitation, instead of 'precipitation'");

	opt_temp_table = G_define_option();
	opt_temp_table->key = "temperature_table";
	opt_temp_table->type = TYPE_STRING;
	opt_temp_table->key_desc = "name";
	opt_temp_table->required = NO;
	opt_temp_table->description =
	    _("Long-format (station_id,date,value) DB table of observed "
	      "temperature, instead of 'temperature'");

	opt_evap_table = G_define_option();
	opt_evap_table->key = "evapotranspiration_table";
	opt_evap_table->type = TYPE_STRING;
	opt_evap_table->key_desc = "name";
	opt_evap_table->required = NO;
	opt_evap_table->description =
	    _("Long-format (station_id,date,value) DB table of potential "
	      "evapotranspiration, instead of 'evapotranspiration'");

	opt_eta_obs_table = G_define_option();
	opt_eta_obs_table->key = "eta_observed_table";
	opt_eta_obs_table->type = TYPE_STRING;
	opt_eta_obs_table->key_desc = "name";
	opt_eta_obs_table->required = NO;
	opt_eta_obs_table->description =
	    _("Long-format (station_id,date,value) DB table of observed "
	      "actual evapotranspiration, instead of 'eta_observed'");

	opt_q_obs_table = G_define_option();
	opt_q_obs_table->key = "discharge_observed_table";
	opt_q_obs_table->type = TYPE_STRING;
	opt_q_obs_table->key_desc = "name";
	opt_q_obs_table->required = NO;
	opt_q_obs_table->description =
	    _("Long-format (station_id,date,value) DB table of observed "
	      "discharge, instead of 'discharge_observed'");

	opt_table_id_column = G_define_option();
	opt_table_id_column->key = "table_id_column";
	opt_table_id_column->type = TYPE_STRING;
	opt_table_id_column->required = NO;
	opt_table_id_column->answer = "station_id";
	opt_table_id_column->description =
	    _("Station-identifier column name shared by all "
	      "*_table options above");

	opt_table_date_column = G_define_option();
	opt_table_date_column->key = "table_date_column";
	opt_table_date_column->type = TYPE_STRING;
	opt_table_date_column->required = NO;
	opt_table_date_column->answer = "date";
	opt_table_date_column->description =
	    _("Date column name shared by all *_table options above");

	opt_table_value_column = G_define_option();
	opt_table_value_column->key = "table_value_column";
	opt_table_value_column->type = TYPE_STRING;
	opt_table_value_column->required = NO;
	opt_table_value_column->answer = "value";
	opt_table_value_column->description =
	    _("Value column name shared by all *_table options above");

	/* GRASS-native alternative to precipitation_table=/temperature_table=/
	 * evapotranspiration_table=: fetch real ERA5(-Land) reanalysis data
	 * directly via t.in.era5 (one call for all three variables) and
	 * reduce it to the same long-format tables via r.hydro.hbv.forcing
	 * -- eta_observed/discharge_observed are calibration targets (real
	 * observations), not meteorological forcing, so they are never
	 * covered by this and must still come from a CSV or table option. */
	opt_era5_start = G_define_option();
	opt_era5_start->key = "era5_start";
	opt_era5_start->type = TYPE_STRING;
	opt_era5_start->required = NO;
	opt_era5_start->description =
	    _("Start date (YYYY-MM-DD) to fetch ERA5(-Land) precipitation/"
	      "temperature/potential_evaporation via t.in.era5, instead of "
	      "precipitation/temperature/evapotranspiration or their "
	      "_table equivalents; requires era5_end, and basins/"
	      "basins_vector (elevation/outlets or basins_vector mode)");

	opt_era5_end = G_define_option();
	opt_era5_end->key = "era5_end";
	opt_era5_end->type = TYPE_STRING;
	opt_era5_end->required = NO;
	opt_era5_end->description =
	    _("End date (YYYY-MM-DD, inclusive) for era5_start");

	opt_era5_area = G_define_option();
	opt_era5_area->key = "era5_area";
	opt_era5_area->type = TYPE_STRING;
	opt_era5_area->required = NO;
	opt_era5_area->key_desc = "north,west,south,east";
	opt_era5_area->description =
	    _("Bounding box in WGS84 degrees for the ERA5(-Land) request; "
	      "default derived from the current region");

	opt_era5_cache_dir = G_define_standard_option(G_OPT_M_DIR);
	opt_era5_cache_dir->key = "era5_cache_dir";
	opt_era5_cache_dir->required = NO;
	opt_era5_cache_dir->description =
	    _("Directory to cache downloaded ERA5(-Land) NetCDF files in "
	      "(default a temporary, run-scoped directory)");

	opt_params = G_define_standard_option(G_OPT_F_INPUT);
	opt_params->key = "parameters";
	opt_params->required = NO;
	opt_params->description =
	    _("CSV of HBV parameter bounds and basin physiography (22 rows "
	      "x n_basins columns)");

	opt_basin_ids = G_define_standard_option(G_OPT_F_INPUT);
	opt_basin_ids->key = "basin_ids";
	opt_basin_ids->required = NO;
	opt_basin_ids->description =
	    _("Text file with one basin identifier per line, used to name "
	      "output files (defaults to Basin01, Basin02, ...); ignored if "
	      "elevation/outlets or basins_vector is given");

	/* GRASS-native basin delineation, in place of the "parameters" CSV:
	 * mirrors r.hydro.hbv.basins' own options (it is spawned internally
	 * when elevation+outlets are given) and reads the resulting basins
	 * vector's attribute table directly -- see basins_input.c. */
	opt_elevation = G_define_standard_option(G_OPT_R_ELEV);
	opt_elevation->key = "elevation";
	opt_elevation->required = NO;
	opt_elevation->description =
	    _("Input DEM for basin delineation (requires outlets; replaces "
	      "'parameters'/'basin_ids' CSV inputs with a GRASS-native "
	      "delineation via r.hydro.hbv.basins)");

	opt_outlets = G_define_standard_option(G_OPT_V_INPUT);
	opt_outlets->key = "outlets";
	opt_outlets->required = NO;
	opt_outlets->description =
	    _("Vector points map, one point per basin outlet/station "
	      "(requires elevation or dem_source)");

	/* GRASS-native alternative to elevation=: fetch a real DEM directly
	 * via r.in.dem (Copernicus GLO-30/GLO-90, no account/API key
	 * needed, truly global land coverage) instead of requiring the
	 * caller to already have a local elevation raster -- for any study
	 * area on Earth with the current region/project set up over it. */
	opt_dem_source = G_define_option();
	opt_dem_source->key = "dem_source";
	opt_dem_source->type = TYPE_STRING;
	opt_dem_source->required = NO;
	opt_dem_source->options = "copernicus_glo30,copernicus_glo90";
	opt_dem_source->description =
	    _("Fetch a DEM via r.in.dem instead of giving elevation= "
	      "directly (requires outlets, not elevation)");

	opt_dem_area = G_define_option();
	opt_dem_area->key = "dem_area";
	opt_dem_area->type = TYPE_STRING;
	opt_dem_area->required = NO;
	opt_dem_area->key_desc = "north,west,south,east";
	opt_dem_area->description =
	    _("Bounding box in WGS84 degrees for dem_source; default "
	      "derived from the current region");

	opt_id_column = G_define_standard_option(G_OPT_DB_COLUMN);
	opt_id_column->key = "id_column";
	opt_id_column->required = NO;
	opt_id_column->answer = "id";
	opt_id_column->description =
	    _("Attribute column in outlets holding the basin identifier");

	opt_threshold = G_define_option();
	opt_threshold->key = "threshold";
	opt_threshold->type = TYPE_INTEGER;
	opt_threshold->required = NO;
	opt_threshold->answer = "1000";
	opt_threshold->description =
	    _("Minimum flow accumulation (cells) for basin delineation");

	opt_snap_radius = G_define_option();
	opt_snap_radius->key = "snap_radius";
	opt_snap_radius->type = TYPE_INTEGER;
	opt_snap_radius->required = NO;
	opt_snap_radius->answer = "3";
	opt_snap_radius->description =
	    _("Maximum distance (cells) to snap outlet points onto the "
	      "derived stream network");

	opt_landcover = G_define_standard_option(G_OPT_R_INPUT);
	opt_landcover->key = "landcover";
	opt_landcover->required = NO;
	opt_landcover->description =
	    _("Land-cover raster used to derive forest/field fraction per "
	      "basin during delineation; omit to use basins_ffo/basins_ffi");

	opt_forest_cats = G_define_option();
	opt_forest_cats->key = "forest_cats";
	opt_forest_cats->type = TYPE_INTEGER;
	opt_forest_cats->multiple = YES;
	opt_forest_cats->required = NO;
	opt_forest_cats->description =
	    _("Land-cover category values considered \"forest\" (required "
	      "if landcover is given)");

	opt_basins_ffo = G_define_option();
	opt_basins_ffo->key = "basins_ffo";
	opt_basins_ffo->type = TYPE_DOUBLE;
	opt_basins_ffo->required = NO;
	opt_basins_ffo->description =
	    _("Uniform forest fraction for every basin during delineation "
	      "when landcover is not given");

	opt_basins_ffi = G_define_option();
	opt_basins_ffi->key = "basins_ffi";
	opt_basins_ffi->type = TYPE_DOUBLE;
	opt_basins_ffi->required = NO;
	opt_basins_ffi->description =
	    _("Uniform field fraction for every basin during delineation "
	      "when landcover is not given");

	opt_precip_station_elev = G_define_option();
	opt_precip_station_elev->key = "precip_station_elevation";
	opt_precip_station_elev->type = TYPE_DOUBLE;
	opt_precip_station_elev->required = NO;
	opt_precip_station_elev->description =
	    _("Elevation of the precipitation station/grid reference point "
	      "during delineation; omit to skip the correction term");

	opt_temp_station_elev = G_define_option();
	opt_temp_station_elev->key = "temp_station_elevation";
	opt_temp_station_elev->type = TYPE_DOUBLE;
	opt_temp_station_elev->required = NO;
	opt_temp_station_elev->description =
	    _("Elevation of the temperature station/grid reference point "
	      "during delineation; omit to skip the correction term");

	opt_et_station_elev = G_define_option();
	opt_et_station_elev->key = "et_station_elevation";
	opt_et_station_elev->type = TYPE_DOUBLE;
	opt_et_station_elev->required = NO;
	opt_et_station_elev->description =
	    _("Elevation of the evapotranspiration station/grid reference "
	      "point during delineation; omit to skip the correction term");

	opt_parameters_template = G_define_standard_option(G_OPT_F_INPUT);
	opt_parameters_template->key = "parameters_template";
	opt_parameters_template->required = NO;
	opt_parameters_template->description =
	    _("CSV of HBV parameter sampling bounds (16 rows: "
	      "fc,beta,lp,alpha,kf,ks,perc,cflux low/high), 1 column "
	      "(broadcast) or one per delineated basin; required when "
	      "elevation/outlets are given");

	opt_hru_bands = G_define_option();
	opt_hru_bands->key = "hru_bands";
	opt_hru_bands->type = TYPE_INTEGER;
	opt_hru_bands->required = NO;
	opt_hru_bands->answer = "1";
	opt_hru_bands->description =
	    _("Number of elevation-band HRUs per basin during delineation "
	      "(1 = one lumped unit per basin, no splitting); only used "
	      "with elevation/outlets");

	opt_basins = G_define_standard_option(G_OPT_R_OUTPUT);
	opt_basins->key = "basins";
	opt_basins->required = NO;
	opt_basins->description =
	    _("Name for output delineated basins raster map (required with "
	      "elevation/outlets)");

	/* deliberately not G_OPT_V_OUTPUT/G_OPT_V_INPUT: this name is used
	 * both ways depending on mode (created fresh when elevation/outlets
	 * are given, read as-is otherwise), so neither the "must not exist"
	 * nor the "must exist" automatic checks apply. */
	opt_basins_vector = G_define_option();
	opt_basins_vector->key = "basins_vector";
	opt_basins_vector->type = TYPE_STRING;
	opt_basins_vector->key_desc = "name";
	opt_basins_vector->required = NO;
	opt_basins_vector->description =
	    _("Name for the basins vector map: created when elevation/"
	      "outlets are given, or read as-is (from a previous "
	      "r.hydro.hbv.basins run) otherwise -- either way, its "
	      "attribute table is used as this run's parameters input");

	opt_output = G_define_standard_option(G_OPT_M_DIR);
	opt_output->key = "output";
	opt_output->required = YES;
	opt_output->answer = ".";
	opt_output->description =
	    _("Output directory for Output*/Basinout*/ETout*/EToutput* CSVs");

	opt_output_tables = G_define_option();
	opt_output_tables->key = "output_tables";
	opt_output_tables->type = TYPE_STRING;
	opt_output_tables->key_desc = "name";
	opt_output_tables->required = NO;
	opt_output_tables->description =
	    _("Name prefix for optional long-format (station_id,day,...) "
	      "GRASS DB tables <prefix>_basinout/<prefix>_etout, written in "
	      "addition to (not instead of) the CSVs above");

	opt_n_basins = G_define_option();
	opt_n_basins->key = "n_basins";
	opt_n_basins->type = TYPE_INTEGER;
	opt_n_basins->required = NO;
	opt_n_basins->description = _("Number of sub-basins / stations");

	opt_n_calib_steps = G_define_option();
	opt_n_calib_steps->key = "n_calib_steps";
	opt_n_calib_steps->type = TYPE_INTEGER;
	opt_n_calib_steps->required = NO;
	opt_n_calib_steps->description =
	    _("Number of daily time steps used in the calibration loop");

	opt_n_days = G_define_option();
	opt_n_days->key = "n_days";
	opt_n_days->type = TYPE_INTEGER;
	opt_n_days->required = NO;
	opt_n_days->description = _("Number of daily rows in the input CSVs");

	opt_n_years = G_define_option();
	opt_n_years->key = "n_years";
	opt_n_years->type = TYPE_INTEGER;
	opt_n_years->required = NO;
	opt_n_years->description = _("Number of years spanned by the data");

	opt_warmup = G_define_option();
	opt_warmup->key = "warmup";
	opt_warmup->type = TYPE_INTEGER;
	opt_warmup->required = NO;
	opt_warmup->description =
	    _("Warm-up period (time steps) excluded from performance stats");

	opt_n_realizations = G_define_option();
	opt_n_realizations->key = "n_realizations";
	opt_n_realizations->type = TYPE_INTEGER;
	opt_n_realizations->required = NO;
	opt_n_realizations->description =
	    _("Number of Monte-Carlo parameter realizations to run");

#define PHYSPARAM(varname, optkey, defval, desc)                              \
	opt_##varname = G_define_option();                                    \
	opt_##varname->key = optkey;                                          \
	opt_##varname->type = TYPE_DOUBLE;                                    \
	opt_##varname->required = NO;                                         \
	opt_##varname->answer = defval;                                       \
	opt_##varname->description = _(desc)

	PHYSPARAM(tcon, "tcon", "86.4", "Unit conversion mm/d -> m3/s factor");
	PHYSPARAM(ecevpfo, "ecevpfo", "1.15", "Forest PET correction factor");
	PHYSPARAM(ecalt, "ecalt", "0.1", "PET elevation correction factor");
	PHYSPARAM(pcalt, "pcalt", "0.2", "Precipitation elevation correction factor");
	PHYSPARAM(tcalt, "tcalt", "0.6", "Temperature elevation correction factor");
	PHYSPARAM(sfcf, "sfcf", "1.0", "Snowfall correction factor");
	PHYSPARAM(foscf, "foscf", "1.0", "Forest snow correction factor");
	PHYSPARAM(rfcf, "rfcf", "1.0", "Rainfall correction factor");
	PHYSPARAM(tt, "tt", "-0.5", "Threshold temperature (deg C)");
	PHYSPARAM(tti, "tti", "2.0", "Temperature interval snow/rain (deg C)");
	PHYSPARAM(dttm, "dttm", "0.54391", "Melt threshold temperature offset");
	PHYSPARAM(cfmax, "cfmax", "3.5", "Degree-day snowmelt factor (mm/d/degC)");
	PHYSPARAM(focfmax, "focfmax", "0.6", "Forest degree-day factor ratio");
	PHYSPARAM(cfr, "cfr", "0.05", "Refreezing coefficient");
	PHYSPARAM(whc, "whc", "0.1", "Water holding capacity of snowpack");
	PHYSPARAM(sgwmax, "sgwmax", "5.0", "Maximum groundwater box storage (mm)");
#undef PHYSPARAM

	flag_legacy = G_define_flag();
	flag_legacy->key = 'l';
	flag_legacy->description =
	    _("Legacy mode: do not reset performance accumulators between "
	      "Monte-Carlo realizations (bit-for-bit reproduces the "
	      "original 2005 Karkheh code, including its unreset-statistics "
	      "quirk)");

	if (G_parser(argc, argv))
		exit(EXIT_FAILURE);

	const char *dataset = opt_dataset->answer;

	if (opt_precip->answer && opt_precip_table->answer)
		G_fatal_error(_("Give either precipitation or "
				"precipitation_table, not both"));
	if (opt_temp->answer && opt_temp_table->answer)
		G_fatal_error(_("Give either temperature or temperature_table, "
				"not both"));
	if (opt_evap->answer && opt_evap_table->answer)
		G_fatal_error(_("Give either evapotranspiration or "
				"evapotranspiration_table, not both"));
	if (opt_eta_obs->answer && opt_eta_obs_table->answer)
		G_fatal_error(_("Give either eta_observed or "
				"eta_observed_table, not both"));
	if (opt_q_obs->answer && opt_q_obs_table->answer)
		G_fatal_error(_("Give either discharge_observed or "
				"discharge_observed_table, not both"));

	if (opt_dem_source->answer) {
		if (opt_elevation->answer)
			G_fatal_error(
			    _("Give either elevation or dem_source, not both"));
		if (!opt_outlets->answer)
			G_fatal_error(
			    _("outlets is required together with dem_source"));

		const char *dem_output = "hbv_dem_source";
		hbv_dem_fetch(opt_dem_source->answer, opt_dem_area->answer,
			      dem_output);
		opt_elevation->answer = G_store(dem_output);
	}

	/* GRASS-native basin delineation replaces the "parameters"/
	 * "basin_ids" CSV inputs whenever elevation+outlets (delineate
	 * fresh, spawning r.hydro.hbv.basins) or basins_vector alone (reuse
	 * an already-delineated vector) are given. */
	if ((opt_elevation->answer == NULL) != (opt_outlets->answer == NULL))
		G_fatal_error(_("elevation and outlets must be given together"));
	if (opt_elevation->answer &&
	    (!opt_basins->answer || !opt_basins_vector->answer))
		G_fatal_error(_("basins and basins_vector are required "
				"together with elevation/outlets"));
	int use_vector_basins =
	    (opt_elevation->answer != NULL) || (opt_basins_vector->answer != NULL);

	int btot;
	float **f7p = NULL;
	char **basin_ids = NULL;
	/* rtot/station_of/station_ids: reporting-station grouping. In CSV
	 * mode (no vector), and in vector mode with hru_bands == 1
	 * (r.hydro.hbv.basins' default), rtot == btot and station_of[b]
	 * == b for every b -- one basin is its own station, today's only
	 * behavior. With hru_bands > 1, several btot rows (HRUs) can share
	 * one station (basin). */
	int rtot;
	int *station_of = NULL;
	char **station_ids = NULL;

	if (use_vector_basins)
		hbv_basins_from_vector(
		    opt_elevation->answer, opt_outlets->answer,
		    opt_id_column->answer, opt_threshold->answer,
		    opt_snap_radius->answer, opt_landcover->answer,
		    opt_forest_cats->answer, opt_basins_ffo->answer,
		    opt_basins_ffi->answer, opt_precip_station_elev->answer,
		    opt_temp_station_elev->answer, opt_et_station_elev->answer,
		    opt_parameters_template->answer, opt_hru_bands->answer,
		    opt_basins->answer, opt_basins_vector->answer, &btot, &f7p,
		    &basin_ids, &rtot, &station_of, &station_ids);
	else {
		btot = resolve_int(opt_n_basins, dataset, 8, 11, 0);
		rtot = btot;
	}

	if ((opt_era5_start->answer == NULL) != (opt_era5_end->answer == NULL))
		G_fatal_error(_("era5_start and era5_end must be given together"));
	if (opt_era5_start->answer) {
		if (opt_precip->answer || opt_precip_table->answer ||
		    opt_temp->answer || opt_temp_table->answer ||
		    opt_evap->answer || opt_evap_table->answer)
			G_fatal_error(
			    _("era5_start/era5_end fetch precipitation/"
			      "temperature/evapotranspiration directly -- "
			      "give none of precipitation[_table]/"
			      "temperature[_table]/evapotranspiration[_table] "
			      "alongside them"));
		if (!opt_basins->answer || !opt_basins_vector->answer)
			G_fatal_error(
			    _("basins and basins_vector are required "
			      "together with era5_start/era5_end (needed "
			      "for r.hydro.hbv.forcing's zonal means)"));

		hbv_era5_forcing(opt_era5_start->answer, opt_era5_end->answer,
				 opt_era5_area->answer,
				 opt_era5_cache_dir->answer, opt_basins->answer,
				 opt_basins_vector->answer,
				 &opt_precip_table->answer,
				 &opt_temp_table->answer,
				 &opt_evap_table->answer);
	}

	int any_table_used = opt_precip_table->answer || opt_temp_table->answer ||
			      opt_evap_table->answer ||
			      opt_eta_obs_table->answer || opt_q_obs_table->answer;

	int ntot = resolve_int(opt_n_realizations, dataset, 45000, 45000, 1);
	int ttot = resolve_int(opt_n_calib_steps, dataset, 916, 8400, 0);
	/* dtot (day count) comes from the *_table inputs' own canonical
	 * date axis when any are used, so n_days isn't needed/resolved in
	 * that case -- see the table_read_pivot() block below. */
	int dtot = any_table_used ? 0 : resolve_int(opt_n_days, dataset, 1910, 8400, 0);
	int ytot = resolve_int(opt_n_years, dataset, 5, 22, 0);
	int tstart = resolve_int(opt_warmup, dataset, 274, 365, 0);
	int th = 3;   /* temporal range moving average high flow */
	int tl = 15;  /* temporal range moving average low flow */
	int TY = 365; /* days per year */
	int TS = 86400; /* seconds per day */

	float tcon = resolve_double(opt_tcon);
	float ecevpfo = resolve_double(opt_ecevpfo);
	float ecalt = resolve_double(opt_ecalt);
	float pcalt = resolve_double(opt_pcalt);
	float tcalt = resolve_double(opt_tcalt);
	float sfcf = resolve_double(opt_sfcf);
	float foscf = resolve_double(opt_foscf);
	float rfcf = resolve_double(opt_rfcf);
	float tt = resolve_double(opt_tt);
	float tti = resolve_double(opt_tti);
	float dttm = resolve_double(opt_dttm);
	float cfmax = resolve_double(opt_cfmax);
	float focfmax = resolve_double(opt_focfmax);
	float cfr = resolve_double(opt_cfr);
	float whc = resolve_double(opt_whc);
	float sgwmax = resolve_double(opt_sgwmax);

	/* optimization criteria constants, unchanged between the two
	 * bundled datasets */
	float rqh = 5.;
	float rql = 0.5;
	float QDB = 100.;
	float KL = 1.30456; /* lower is 10 years and upper is 100 years */
	float KU = 3.13668; /* lower is 10 years and upper is 100 years */

	char *precip_path = opt_precip_table->answer
				? NULL
				: resolve_path(opt_precip, dataset, "precip.csv");
	char *temp_path = opt_temp_table->answer
			      ? NULL
			      : resolve_path(opt_temp, dataset, "temp.csv");
	char *evap_path = opt_evap_table->answer
			      ? NULL
			      : resolve_path(opt_evap, dataset, "evap.csv");
	char *eta_obs_path =
	    opt_eta_obs_table->answer
		? NULL
		: resolve_path(opt_eta_obs, dataset, "etobs.csv");
	char *q_obs_path =
	    opt_q_obs_table->answer
		? NULL
		: resolve_path(opt_q_obs, dataset, "dischargeobs.csv");

	if ((!precip_path && !opt_precip_table->answer) ||
	    (!temp_path && !opt_temp_table->answer) ||
	    (!evap_path && !opt_evap_table->answer) ||
	    (!eta_obs_path && !opt_eta_obs_table->answer) ||
	    (!q_obs_path && !opt_q_obs_table->answer))
		G_fatal_error(
		    _("precipitation/temperature/evapotranspiration/"
		      "eta_observed/discharge_observed (or their _table "
		      "equivalents) are required when dataset=custom"));

	if (!use_vector_basins) {
		/* legacy path: basin physiography/parameter bounds and basin
		 * IDs come from CSV/text files instead of a basins vector. */
		char *params_path =
		    resolve_path(opt_params, dataset, "param_sto.csv");
		char *basin_ids_path =
		    opt_basin_ids->answer
			? G_store(opt_basin_ids->answer)
			: (strcmp(dataset, "custom") != 0
			       ? bundled_data_path(dataset, "basin_ids.txt")
			       : NULL);

		if (!params_path)
			G_fatal_error(_("parameters is required when "
					"dataset=custom and elevation/"
					"outlets/basins_vector are not "
					"given"));

		basin_ids = load_basin_ids(basin_ids_path, btot);
		f7p = readcsv(params_path, btot, 22);

		/* no HRU grouping in CSV mode: every row is its own station */
		station_of = (int *)G_malloc(btot * sizeof(int));
		station_ids = basin_ids;
		for (b = 0; b < btot; b++)
			station_of[b] = b;
	}

	ensure_dir(opt_output->answer);

	/* storage soil moisture,surface water,ground water */
	float **ssm = af2d(btot, ttot + 1);
	float **ssw = af2d(btot, ttot + 1);
	float **sgw = af2d(btot, ttot + 1);
	/* storage snow pack,melting water */
	float **ssp = af2d(btot, ttot + 1);
	float **smw = af2d(btot, ttot + 1);
	float **sgwx = af2d(btot, ttot + 1);
	/* precipitation corrected,temperature corrected,snowfall,rainfall */
	float **p = af2d(btot, ttot);
	float **ta = af2d(btot, ttot);
	float **s = af2d(btot, ttot);
	float **r = af2d(btot, ttot);
	/* snow melt,refreezing rate,infiltration,PET corrected */
	float **sm = af2d(btot, ttot);
	float **sr = af2d(btot, ttot);
	float **inp = af2d(btot, ttot);
	float **etp = af2d(btot, ttot);
	/* direct discharge,indirect discharge,slow flow,fast flow */
	float **qd = af2d(btot, ttot);
	float **qin = af2d(btot, ttot);
	float **qs = af2d(btot, ttot);
	float **qf = af2d(btot, ttot);
	/* total discharge,capillary rise,actual evapotranspiration */
	float **qt = af2d(btot, ttot);
	float **qc = af2d(btot, ttot);
	float **eta = af2d(btot, ttot);

	float **fc = af2d(btot, ntot);
	float **beta = af2d(btot, ntot);
	float **lp = af2d(btot, ntot);
	float **alpha = af2d(btot, ntot);
	float **kf = af2d(btot, ntot);
	float **ks = af2d(btot, ntot);
	float **perc = af2d(btot, ntot);
	float **cflux = af2d(btot, ntot);

	/* river discharge gathered from qt after each hbv run (HRU-level) */
	float **qr = af2d(btot, ttot);
	/* qr/eta aggregated per reporting station (summed discharge,
	 * area-weighted-averaged ETa across an HRU group) -- this is what
	 * hbv_performance() actually scores against qo/eto, both of which
	 * are inherently station-level (gauge/reference observations). For
	 * hru_bands == 1 this is an identity copy of qr/eta. */
	float **qr_station = af2d(rtot, ttot);
	float **eta_station = af2d(rtot, ttot);

	float *m = af1d(rtot);
	float *mtot = af1d(rtot);
	float *qotot = af1d(rtot);
	float *yotot = af1d(rtot);
	float *etotot = af1d(rtot);

	/* m[] is cast to an array index in hbv_performance() and persists
	 * (unreset) across every Monte-Carlo realization, so it must start
	 * at exactly 0 -- af1d() only malloc()s, it never zeroes. Relying on
	 * malloc'd memory happening to already be zero (true for a
	 * process's first, small allocations on most systems, which is why
	 * this went unnoticed) segfaulted here once the allocation
	 * footprint grew large enough (11 basins x 8400 days x 45000
	 * realizations) that this assumption no longer held. */
	for (b = 0; b < rtot; b++) {
		m[b] = 0.;
		mtot[b] = 0.;
	}

	float **NS = af2d(rtot, ntot);
	float **NSET = af2d(rtot, ntot);
	float **NSH = af2d(rtot, ntot);
	float **NSL = af2d(rtot, ntot);
	float **RVE = af2d(rtot, ntot);
	float **RVEET = af2d(rtot, ntot);
	float **RMAEH = af2d(rtot, ntot);
	float **RMAEL = af2d(rtot, ntot);
	float **REVE = af2d(rtot, ntot);
	float *dqd = af1d(ntot);

	float **MNS = af2d(rtot, ntot);
	float **MNSET = af2d(rtot, ntot);
	float **MNSH = af2d(rtot, ntot);
	float **MNSL = af2d(rtot, ntot);
	float **MRVE = af2d(rtot, ntot);
	float **MRVEET = af2d(rtot, ntot);
	float **MRMAEH = af2d(rtot, ntot);
	float **MRMAEL = af2d(rtot, ntot);
	float **MREVE = af2d(rtot, ntot);
	float *Mdqd = af1d(ntot);

	float *difq = af1d(rtot);
	float *dq = af1d(rtot);
	float *dqh = af1d(rtot);
	float *dql = af1d(rtot);
	float *mqo = af1d(rtot);
	float *meto = af1d(rtot);
	float *difet = af1d(rtot);
	float *deta = af1d(rtot);
	float *veto = af1d(rtot);
	float *mqoh = af1d(rtot);
	float *mqol = af1d(rtot);
	float *vqo = af1d(rtot);
	float *vqoh = af1d(rtot);
	float *vqol = af1d(rtot);
	float **qrh = af2d(rtot, ttot);
	float **qoh = af2d(rtot, ttot);
	float **qrl = af2d(rtot, ttot);
	float **qol = af2d(rtot, ttot);
	float *dqha = af1d(rtot);
	float *dqla = af1d(rtot);
	float *qdo = af1d(ytot);
	float *qdr = af1d(ytot);
	float *qodef = af1d(ntot);
	float *qrdef = af1d(ntot);
	/* selected-realization parameter draws are per computation unit
	 * (HRU), matching hbv_performance()'s `b<btot-1` storage loop --
	 * not rtot, which is what NS/RVE-family scores below use. */
	float **mfc = af2d(btot, ntot);
	float **mbeta = af2d(btot, ntot);
	float **mlp = af2d(btot, ntot);
	float **malpha = af2d(btot, ntot);
	float **mkf = af2d(btot, ntot);
	float **mks = af2d(btot, ntot);
	float **mperc = af2d(btot, ntot);
	float **mcflux = af2d(btot, ntot);
	float *tyear = af1d(ttot);
	float **qom = af2d(rtot, ytot);
	float **qrm = af2d(rtot, ytot);
	float *mqom = af1d(rtot);
	float *mqrm = af1d(rtot);
	float *sqom = af1d(rtot);
	float *sqmr = af1d(rtot);

	/* observed precipitation,potential ET,discharge,temperature,ETa.
	 * precipitation/temperature/PET are station-level (one column per
	 * reporting station, matching the input CSV/table regardless of
	 * HRU splitting) and get expanded to HRU width just below;
	 * discharge/ETa observations are inherently station-level
	 * (gauge/reference measurements) and stay that way throughout. */
	char **table_dates = NULL;
	int table_n_dates = 0;
	const char *tic = opt_table_id_column->answer;
	const char *tdc = opt_table_date_column->answer;
	const char *tvc = opt_table_value_column->answer;

	/* every *_table read happens before any readcsv() call below: dtot
	 * is only known (from table_n_dates) once at least one *_table has
	 * been read, and any_table_used means readcsv()'s dtot argument
	 * would otherwise still be the placeholder 0 it was given above --
	 * reading a real CSV with col=0 corrupts af2d()'s allocation and
	 * segfaults. This can happen with a mix of *_table and plain CSV
	 * inputs (e.g. era5_start/era5_end for precipitation/temperature/
	 * evapotranspiration, plain CSV for eta_observed/
	 * discharge_observed). */
	float **po_station =
	    opt_precip_table->answer
		? table_read_pivot(opt_precip_table->answer, tic, tdc, tvc,
				    station_ids, rtot, &table_dates,
				    &table_n_dates)
		: NULL;
	float **etpo_station =
	    opt_evap_table->answer
		? table_read_pivot(opt_evap_table->answer, tic, tdc, tvc,
				    station_ids, rtot, &table_dates,
				    &table_n_dates)
		: NULL;
	float **qo = opt_q_obs_table->answer
			 ? table_read_pivot(opt_q_obs_table->answer, tic, tdc,
					     tvc, station_ids, rtot,
					     &table_dates, &table_n_dates)
			 : NULL;
	float **tm_station =
	    opt_temp_table->answer
		? table_read_pivot(opt_temp_table->answer, tic, tdc, tvc,
				    station_ids, rtot, &table_dates,
				    &table_n_dates)
		: NULL;
	float **eto = opt_eta_obs_table->answer
			  ? table_read_pivot(opt_eta_obs_table->answer, tic,
					      tdc, tvc, station_ids, rtot,
					      &table_dates, &table_n_dates)
			  : NULL;

	if (any_table_used)
		dtot = table_n_dates;

	if (!po_station)
		po_station = readcsv(precip_path, rtot, dtot);
	if (!etpo_station)
		etpo_station = readcsv(evap_path, rtot, dtot);
	if (!qo)
		qo = readcsv(q_obs_path, rtot, dtot);
	if (!tm_station)
		tm_station = readcsv(temp_path, rtot, dtot);
	if (!eto)
		eto = readcsv(eta_obs_path, rtot, dtot);

	/* expand station-level forcing to HRU width -- each HRU reuses its
	 * parent station's time series (hbv_model.c already applies each
	 * HRU's own dep/det/dee elevation correction on top of it, so no
	 * separate per-HRU forcing data is needed). For hru_bands == 1 this
	 * is an identity copy (station_of[b] == b, rtot == btot). */
	float **po = af2d(btot, dtot);
	float **etpo = af2d(btot, dtot);
	float **tm = af2d(btot, dtot);
	for (b = 0; b < btot; b++) {
		for (t = 0; t < dtot; t++) {
			po[b][t] = po_station[station_of[b]][t];
			etpo[b][t] = etpo_station[station_of[b]][t];
			tm[b][t] = tm_station[station_of[b]][t];
		}
	}
	/* f7p (basin physiography + HBV parameter sampling bounds) was
	 * already populated above, either from a basins vector's attribute
	 * table (hbv_basins_from_vector) or from the "parameters" CSV. */

	float *lfc = af1d(btot);
	float *hfc = af1d(btot);
	float *lbeta = af1d(btot);
	float *hbeta = af1d(btot);
	float *llp = af1d(btot);
	float *hlp = af1d(btot);
	float *lalpha = af1d(btot);
	float *halpha = af1d(btot);
	float *lkf = af1d(btot);
	float *hkf = af1d(btot);
	float *lks = af1d(btot);
	float *hks = af1d(btot);
	float *lperc = af1d(btot);
	float *hperc = af1d(btot);
	float *lcflux = af1d(btot);
	float *hcflux = af1d(btot);

	float *area = af1d(btot);
	float *ffo = af1d(btot);
	float *ffi = af1d(btot);
	float *dep = af1d(btot);
	float *det = af1d(btot);
	float *dee = af1d(btot);

	for (b = 0; b < btot; b++) {
		lfc[b] = f7p[b][0];
		hfc[b] = f7p[b][1];
		lbeta[b] = f7p[b][2];
		hbeta[b] = f7p[b][3];
		llp[b] = f7p[b][4];
		hlp[b] = f7p[b][5];
		lalpha[b] = f7p[b][6];
		halpha[b] = f7p[b][7];
		lkf[b] = f7p[b][8];
		hkf[b] = f7p[b][9];
		lks[b] = f7p[b][10];
		hks[b] = f7p[b][11];
		lperc[b] = f7p[b][12];
		hperc[b] = f7p[b][13];
		lcflux[b] = f7p[b][14];
		hcflux[b] = f7p[b][15];
		area[b] = f7p[b][16];
		ffo[b] = f7p[b][17];
		ffi[b] = f7p[b][18];
		dep[b] = f7p[b][19];
		det[b] = f7p[b][20];
		dee[b] = f7p[b][21];
	}

	/* per-station total HRU area, for area-weighting eta_station below
	 * (invariant across realizations, computed once) */
	float *station_area = af1d(rtot);
	for (b = 0; b < rtot; b++)
		station_area[b] = 0.;
	for (b = 0; b < btot; b++)
		station_area[station_of[b]] += area[b];

	for (n = 0; n < ntot; n++) {
		G_percent(n, ntot, 1);
		for (b = 0; b < btot; b++) {
			fc[b][n] = lfc[b] + (hfc[b] - lfc[b]) * (rand() % 1000) / 1000;
			beta[b][n] =
			    lbeta[b] + (hbeta[b] - lbeta[b]) * (rand() % 1000) / 1000;
			lp[b][n] = llp[b] + (hlp[b] - llp[b]) * (rand() % 1000) / 1000;
			alpha[b][n] = lalpha[b] +
				      (halpha[b] - lalpha[b]) * (rand() % 1000) / 1000;
			kf[b][n] = lkf[b] + (hkf[b] - lkf[b]) * (rand() % 1000) / 1000;
			ks[b][n] = lks[b] + (hks[b] - lks[b]) * (rand() % 1000) / 1000;
			perc[b][n] =
			    lperc[b] + (hperc[b] - lperc[b]) * (rand() % 1000) / 1000;
			cflux[b][n] = lcflux[b] +
				      (hcflux[b] - lcflux[b]) * (rand() % 1000) / 1000;
			ssp[b][0] = 0.0;
			smw[b][0] = 0.0;
			ssm[b][0] = MIN(15.0, fc[b][n]);
			ssw[b][0] = 15.0;
			sgw[b][0] = 15.0;
			ssp[b][1] = ssp[b][0];
			smw[b][1] = smw[b][0];
			ssm[b][1] = ssm[b][0];
			ssw[b][1] = ssw[b][0];
			sgw[b][1] = sgw[b][0];
			qf[b][0] = kf[b][n] * pow(ssw[b][0], 1 + alpha[b][n]);
			qs[b][0] = ks[b][n] * sgw[b][0];
			qt[b][0] = (qf[b][0] + qs[b][0]) * area[b] / tcon;
		}

#pragma omp parallel for default(shared) private(t, b)
		for (t = 0; t < ttot; t++) {
			for (b = 0; b < btot; b++) {
				hbv_model(n, t, b, p, po, pcalt, dep, ta, tm, tcalt, det,
					  tt, tti, s, ffo, foscf, ffi, sfcf, r, rfcf, sm,
					  cfmax, dttm, ssp, cfr, focfmax, smw, inp, sr,
					  whc, qd, ssm, fc, qin, beta, etp, etpo, ecevpfo,
					  ecalt, dee, lp, qc, cflux, qf, kf, ssw, alpha,
					  qs, ks, sgw, qt, area, tcon, sgwx, sgwmax, perc,
					  eta);
				qr[b][t] = qt[b][t];
			}
		}

		if (!flag_legacy->answer) {
			/* reset per-realization accumulators; the original 2005
			 * code never did this, silently accumulating stats
			 * across every Monte-Carlo run (see -l flag). These are
			 * all rtot-sized (station-level scoring inputs), not
			 * btot-sized. */
			for (b = 0; b < rtot; b++) {
				difq[b] = 0.;
				dq[b] = 0.;
				dqh[b] = 0.;
				dql[b] = 0.;
				mqo[b] = 0.;
				qotot[b] = 0.;
				meto[b] = 0.;
				etotot[b] = 0.;
				mqoh[b] = 0.;
				mqol[b] = 0.;
				vqo[b] = 0.;
				vqoh[b] = 0.;
				vqol[b] = 0.;
				dqha[b] = 0.;
				dqla[b] = 0.;
				yotot[b] = 0.;
				mqom[b] = 0.;
				mqrm[b] = 0.;
				sqom[b] = 0.;
				sqmr[b] = 0.;
				difet[b] = 0.;
				deta[b] = 0.;
				veto[b] = 0.;
				for (t = 0; t < ytot; t++) {
					qom[b][t] = 0.;
					qrm[b][t] = 0.;
				}
				for (t = 0; t < ttot; t++) {
					qrh[b][t] = 0.;
					qoh[b][t] = 0.;
					qrl[b][t] = 0.;
					qol[b][t] = 0.;
				}
			}
			for (t = 0; t < ytot; t++) {
				qdo[t] = 0.;
				qdr[t] = 0.;
			}
		}

		if (!flag_legacy->answer) {
			/* clip non-finite/negative discharge before scoring: the
			 * physics can produce inf/NaN (or small negative values)
			 * for extreme sampled parameter combinations, which
			 * would otherwise poison NS/RVE stats. Skipped in legacy
			 * mode, which never applied this clipping. */
			for (b = 0; b < btot; b++) {
				for (t = 0; t < ttot; t++) {
					if (!isfinite(qr[b][t]) || qr[b][t] < 0.)
						qr[b][t] = 0.;
				}
			}
		}

		/* aggregate HRU-level qr/eta to the station level that
		 * hbv_performance() scores against qo/eto -- summed discharge
		 * (each HRU's own area-scaled contribution), area-weighted
		 * mean ETa. Identity copy when hru_bands == 1. */
		for (b = 0; b < rtot; b++) {
			for (t = 0; t < ttot; t++) {
				qr_station[b][t] = 0.;
				eta_station[b][t] = 0.;
			}
		}
		for (b = 0; b < btot; b++) {
			int st = station_of[b];
			float w = area[b] / station_area[st];

			for (t = 0; t < ttot; t++) {
				qr_station[st][t] += qr[b][t];
				eta_station[st][t] += eta[b][t] * w;
			}
		}

#pragma omp barrier
		hbv_performance(
		    n, rtot, ttot, ytot, btot, ntot, tstart, qo, mqo, qotot, eto,
		    meto, etotot, difq, qr_station, dq, dqh, dql, vqo, vqoh, vqol,
		    difet, eta_station, deta, veto, th, qrh, qoh, mqoh, tl, qrl,
		    qol, rql, rqh, dqha,
		    dqla, mqol, tyear, qom, qrm, yotot, mqom, mqrm, sqom, sqmr, NS,
		    NSET, NSH, NSL, RVE, RVEET, RMAEH, RMAEL, REVE, KL, KU, qodef,
		    qrdef, TY, qdo, qdr, QDB, TS, fc, mfc, mbeta, beta, mlp, lp,
		    malpha, alpha, mkf, kf, mks, ks, mperc, perc, mcflux, cflux, MNS,
		    MNSET, MNSH, MNSL, MRVE, MRVEET, MRMAEH, MRMAEL, MREVE, dqd,
		    Mdqd, mtot, m);
	}
	G_percent(ntot, ntot, 1);

	/* the original 2005 report writer used a single hardcoded mtot[1]
	 * as the "selected realizations" loop bound for every basin
	 * (rather than each basin's own mtot[b]) -- reproduce that quirk
	 * verbatim under -l instead of the corrected per-basin bound. */
	float *mtot_report = mtot;
	if (flag_legacy->answer) {
		mtot_report = af1d(rtot);
		for (b = 0; b < rtot; b++)
			mtot_report[b] = mtot[1];
	}

	char *table_basinout_csv = NULL, *table_etout_csv = NULL;

	if (opt_output_tables->answer) {
		/* OGR's CSV driver needs an actual .csv extension to
		 * recognize the file -- G_tempfile()'s own name has none. */
		table_basinout_csv = tempfile_with_suffix(".csv");
		table_etout_csv = tempfile_with_suffix(".csv");
	}

	/* one report row per reporting station (rtot), not per computation
	 * unit (btot): qr_station/eta_station/po_station are the
	 * station-level aggregates from above, and mfc/mbeta/.../MNS/...
	 * (allocated at btot/rtot respectively, btot >= rtot always) are
	 * read here only up to index rtot-1 -- for hru_bands > 1 this
	 * reports each station's first (lowest elevation band) HRU's
	 * parameter draw as that station's representative one, since
	 * hbv_performance()'s selection bookkeeping (m[]/mtot[]) is itself
	 * only ever advanced per station, not per HRU (see basins_input.c/
	 * hbv_performance.c's own btot-1-vs-rtot-1 loop split). */
	hbv_report(opt_output->answer, station_ids, rtot, mtot_report, ttot,
		   qo, qr_station, po_station, eta_station, eto, mfc, mbeta,
		   mlp, malpha, mkf, mks, mperc, mcflux, MNSET, MRVEET, MNS,
		   MNSH, MNSL, MRVE, MRMAEH, MRMAEL, MREVE, table_basinout_csv,
		   table_etout_csv);

	if (opt_output_tables->answer) {
		import_long_csv_as_table(table_basinout_csv,
					  opt_output_tables->answer,
					  "_basinout");
		import_long_csv_as_table(table_etout_csv,
					  opt_output_tables->answer, "_etout");
	}

	exit(EXIT_SUCCESS);
}
