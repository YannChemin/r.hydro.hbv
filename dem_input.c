#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/spawn.h>

#include "dem_input.h"

void hbv_dem_fetch(const char *dem_source, const char *dem_area,
		    const char *dem_output) {
	const char *argv[8];
	int argc = 0;
	char source_arg[32], area_arg[128], output_arg[GNAME_MAX + 16];

	snprintf(source_arg, sizeof(source_arg), "source=%s", dem_source);
	snprintf(output_arg, sizeof(output_arg), "output=%s", dem_output);

	argv[argc++] = "r.in.dem";
	argv[argc++] = source_arg;
	argv[argc++] = output_arg;
	if (dem_area) {
		snprintf(area_arg, sizeof(area_arg), "area=%s", dem_area);
		argv[argc++] = area_arg;
	}
	if (G_get_overwrite())
		argv[argc++] = "--overwrite";
	argv[argc++] = NULL;

	G_message(_("Fetching a DEM via r.in.dem (this calls out to the "
		    "public Copernicus DEM AWS Open Data bucket)..."));
	int ret = G_vspawn_ex(argv[0], argv);
	if (ret != 0)
		G_fatal_error(_("r.in.dem failed (exit code %d)"), ret);
}
