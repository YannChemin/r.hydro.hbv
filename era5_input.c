#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/spawn.h>

#include "era5_input.h"

static void spawn_or_die(const char **argv, const char *what) {
	int ret = G_vspawn_ex(argv[0], argv);
	if (ret != 0)
		G_fatal_error(_("%s failed (exit code %d)"), what, ret);
}

void hbv_era5_forcing(const char *era5_start, const char *era5_end,
		      const char *era5_area, const char *era5_cache_dir,
		      const char *basins_raster, const char *basins_vector,
		      char **precip_table_out, char **temp_table_out,
		      char **evap_table_out) {
	const char *prefix = "hbv_era5";
	const char *argv[16];
	int argc;
	char start_arg[64], end_arg[64], area_arg[128], cache_arg[GPATH_MAX + 16];
	char output_prefix_arg[64];

	argc = 0;
	argv[argc++] = "t.in.era5";
	argv[argc++] = "variables=precipitation,temperature,"
		       "potential_evaporation";
	snprintf(start_arg, sizeof(start_arg), "start=%s", era5_start);
	argv[argc++] = start_arg;
	snprintf(end_arg, sizeof(end_arg), "end=%s", era5_end);
	argv[argc++] = end_arg;
	if (era5_area) {
		snprintf(area_arg, sizeof(area_arg), "area=%s", era5_area);
		argv[argc++] = area_arg;
	}
	snprintf(output_prefix_arg, sizeof(output_prefix_arg),
		 "output_prefix=%s", prefix);
	argv[argc++] = output_prefix_arg;
	if (era5_cache_dir) {
		snprintf(cache_arg, sizeof(cache_arg), "cache_dir=%s",
			 era5_cache_dir);
		argv[argc++] = cache_arg;
	}
	if (G_get_overwrite())
		argv[argc++] = "--overwrite";
	argv[argc++] = NULL;

	G_message(_("Fetching ERA5(-Land) precipitation/temperature/"
		    "potential_evaporation via t.in.era5 (this calls out to "
		    "the Copernicus Climate Data Store and may take a "
		    "while)..."));
	spawn_or_die(argv, "t.in.era5");

	{
		static const char *var_keys[3] = {"precipitation", "temperature",
						   "potential_evaporation"};
		char **table_outs[3];
		int i;

		table_outs[0] = precip_table_out;
		table_outs[1] = temp_table_out;
		table_outs[2] = evap_table_out;

		for (i = 0; i < 3; i++) {
			char strds_arg[128], table_name[128], output_table_arg[160];
			const char *fargv[10];
			int fargc = 0;

			snprintf(strds_arg, sizeof(strds_arg), "strds=%s_%s",
				 prefix, var_keys[i]);
			snprintf(table_name, sizeof(table_name), "%s_%s_table",
				 prefix, var_keys[i]);
			snprintf(output_table_arg, sizeof(output_table_arg),
				 "output_table=%s", table_name);

			fargv[fargc++] = "r.hydro.hbv.forcing";
			fargv[fargc++] = strds_arg;
			{
				char basins_arg[GNAME_MAX + 16],
				    basins_vector_arg[GNAME_MAX + 16];

				snprintf(basins_arg, sizeof(basins_arg),
					 "basins=%s", basins_raster);
				fargv[fargc++] = basins_arg;
				snprintf(basins_vector_arg,
					 sizeof(basins_vector_arg),
					 "basins_vector=%s", basins_vector);
				fargv[fargc++] = basins_vector_arg;
				fargv[fargc++] = output_table_arg;
				if (G_get_overwrite())
					fargv[fargc++] = "--overwrite";
				fargv[fargc++] = NULL;

				spawn_or_die(fargv, "r.hydro.hbv.forcing");
			}

			*table_outs[i] = G_store(table_name);
		}
	}
}
