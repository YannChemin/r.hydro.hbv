/* Spawns r.in.dem to fetch a real, globally-available DEM (Copernicus
 * GLO-30/GLO-90, no account/API key needed) covering the current
 * region (or dem_area, if given) into raster dem_output -- an
 * alternative to requiring the caller to already have a local
 * elevation raster before delineating basins. */
void hbv_dem_fetch(const char *dem_source, const char *dem_area,
		    const char *dem_output);
