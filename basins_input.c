#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/spawn.h>

#include "arrays.h"
#include "csvline.h"
#include "basins_input.h"

#define N_PARAM_COLS 22

static const char *PARAM_COLUMNS[N_PARAM_COLS] = {
    "fc_lo",  "fc_hi",  "beta_lo",  "beta_hi",  "lp_lo",    "lp_hi",
    "alpha_lo", "alpha_hi", "kf_lo",  "kf_hi",  "ks_lo",    "ks_hi",
    "perc_lo",  "perc_hi",  "cflux_lo", "cflux_hi", "area_km2", "ffo",
    "ffi",      "dep",      "det",      "dee",
};

static void spawn_or_die(const char **argv, const char *what) {
    int ret = G_vspawn_ex(argv[0], argv);
    if (ret != 0)
        G_fatal_error(_("%s failed (exit code %d)"), what, ret);
}

static void delineate_with_basins_module(
    const char *elevation, const char *outlets, const char *id_column,
    const char *threshold, const char *snap_radius, const char *landcover,
    const char *forest_cats, const char *ffo_default,
    const char *ffi_default, const char *precip_station_elevation,
    const char *temp_station_elevation, const char *et_station_elevation,
    const char *parameters_template, const char *hru_bands,
    const char *basins_raster, const char *basins_vector) {
    const char *argv[40];
    int argc = 0;
    char elevation_arg[GNAME_MAX + 16], outlets_arg[GNAME_MAX + 16];
    char id_column_arg[GNAME_MAX + 16], threshold_arg[64],
        snap_radius_arg[64];
    char landcover_arg[GNAME_MAX + 16], forest_cats_arg[512];
    char ffo_arg[64], ffi_arg[64];
    char precip_elev_arg[64], temp_elev_arg[64], et_elev_arg[64];
    char parameters_template_arg[GPATH_MAX + 32], hru_bands_arg[64];
    char basins_arg[GNAME_MAX + 16], basins_vector_arg[GNAME_MAX + 16];

    if (!parameters_template)
        G_fatal_error(_("parameters_template is required when delineating "
                        "basins from elevation/outlets (HBV parameter "
                        "sampling bounds cannot be derived from a DEM)"));

    argv[argc++] = "r.hydro.hbv.basins";

    snprintf(elevation_arg, sizeof(elevation_arg), "elevation=%s", elevation);
    argv[argc++] = elevation_arg;
    snprintf(outlets_arg, sizeof(outlets_arg), "outlets=%s", outlets);
    argv[argc++] = outlets_arg;
    if (id_column) {
        snprintf(id_column_arg, sizeof(id_column_arg), "id_column=%s",
                 id_column);
        argv[argc++] = id_column_arg;
    }
    snprintf(threshold_arg, sizeof(threshold_arg), "threshold=%s", threshold);
    argv[argc++] = threshold_arg;
    if (snap_radius) {
        snprintf(snap_radius_arg, sizeof(snap_radius_arg), "snap_radius=%s",
                 snap_radius);
        argv[argc++] = snap_radius_arg;
    }
    if (landcover) {
        snprintf(landcover_arg, sizeof(landcover_arg), "landcover=%s",
                 landcover);
        argv[argc++] = landcover_arg;
    }
    if (forest_cats) {
        snprintf(forest_cats_arg, sizeof(forest_cats_arg), "forest_cats=%s",
                 forest_cats);
        argv[argc++] = forest_cats_arg;
    }
    if (ffo_default) {
        snprintf(ffo_arg, sizeof(ffo_arg), "ffo=%s", ffo_default);
        argv[argc++] = ffo_arg;
    }
    if (ffi_default) {
        snprintf(ffi_arg, sizeof(ffi_arg), "ffi=%s", ffi_default);
        argv[argc++] = ffi_arg;
    }
    if (precip_station_elevation) {
        snprintf(precip_elev_arg, sizeof(precip_elev_arg),
                 "precip_station_elevation=%s", precip_station_elevation);
        argv[argc++] = precip_elev_arg;
    }
    if (temp_station_elevation) {
        snprintf(temp_elev_arg, sizeof(temp_elev_arg),
                 "temp_station_elevation=%s", temp_station_elevation);
        argv[argc++] = temp_elev_arg;
    }
    if (et_station_elevation) {
        snprintf(et_elev_arg, sizeof(et_elev_arg),
                 "et_station_elevation=%s", et_station_elevation);
        argv[argc++] = et_elev_arg;
    }
    snprintf(parameters_template_arg, sizeof(parameters_template_arg),
             "parameters_template=%s", parameters_template);
    argv[argc++] = parameters_template_arg;
    if (hru_bands) {
        snprintf(hru_bands_arg, sizeof(hru_bands_arg), "hru_bands=%s",
                 hru_bands);
        argv[argc++] = hru_bands_arg;
    }
    snprintf(basins_arg, sizeof(basins_arg), "basins=%s", basins_raster);
    argv[argc++] = basins_arg;
    snprintf(basins_vector_arg, sizeof(basins_vector_arg),
             "basins_vector=%s", basins_vector);
    argv[argc++] = basins_vector_arg;
    if (G_get_overwrite())
        argv[argc++] = "--overwrite";
    argv[argc++] = NULL;

    G_message(_("Delineating basins via r.hydro.hbv.basins..."));
    spawn_or_die(argv, "r.hydro.hbv.basins");
}

static void read_basins_table(const char *basins_vector, int *btot_out,
                               float ***f7p_out, char ***basin_ids_out,
                               int *rtot_out, int **station_of_out,
                               char ***station_ids_out) {
    char columns[1024];
    const char *argv[16];
    int argc = 0;
    char map_arg[GNAME_MAX + 16], columns_arg[1024 + 16];
    char *tmpfile;
    FILE *f;
    char buf[4096];
    int n_rows = 0;
    int cats[4096];
    char *ids_tmp[4096];
    float vals_tmp[4096][N_PARAM_COLS];
    int i, j, k;

    snprintf(columns, sizeof(columns), "cat,basin_id");
    for (i = 0; i < N_PARAM_COLS; i++) {
        strcat(columns, ",");
        strcat(columns, PARAM_COLUMNS[i]);
    }

    tmpfile = G_tempfile();

    argv[argc++] = "v.db.select";
    snprintf(map_arg, sizeof(map_arg), "map=%s", basins_vector);
    argv[argc++] = map_arg;
    snprintf(columns_arg, sizeof(columns_arg), "columns=%s", columns);
    argv[argc++] = columns_arg;
    argv[argc++] = "format=csv";
    argv[argc++] = SF_REDIRECT_FILE;
    argv[argc++] = SF_STDOUT;
    argv[argc++] = SF_MODE_OUT;
    argv[argc++] = tmpfile;
    argv[argc++] = NULL;

    if (G_vspawn_ex(argv[0], argv) != 0) {
        remove(tmpfile);
        G_fatal_error(_("Unable to read attribute table of basins vector "
                        "<%s> (missing basin_id/physiography/parameter-"
                        "bound columns? see r.hydro.hbv.basins)"),
                      basins_vector);
    }

    f = fopen(tmpfile, "r");
    if (!f)
        G_fatal_error(_("Unable to open temporary file <%s>"), tmpfile);

    if (!G_getl(buf, sizeof(buf), f)) /* header */
        G_fatal_error(_("Empty attribute table for basins vector <%s>"),
                      basins_vector);

    while (G_getl(buf, sizeof(buf), f)) {
        char *fields[2 + N_PARAM_COLS];
        int n_fields = csv_split_line(buf, fields, 2 + N_PARAM_COLS);

        if (n_fields != 2 + N_PARAM_COLS)
            G_fatal_error(
                _("Malformed attribute row in basins vector <%s>: "
                  "expected %d fields, got %d"),
                basins_vector, 2 + N_PARAM_COLS, n_fields);
        if (n_rows >= 4096)
            G_fatal_error(_("Too many basins (max 4096 supported)"));

        cats[n_rows] = atoi(fields[0]);
        ids_tmp[n_rows] = G_store(fields[1]);
        for (j = 0; j < N_PARAM_COLS; j++)
            vals_tmp[n_rows][j] = (float)atof(fields[2 + j]);
        n_rows++;
    }
    fclose(f);
    remove(tmpfile);

    if (n_rows == 0)
        G_fatal_error(_("Basins vector <%s> has no populated basin rows"),
                      basins_vector);

    /* sort by ascending category (insertion sort -- basin counts are
     * always small: tens, not thousands) */
    for (i = 1; i < n_rows; i++) {
        int cat_i = cats[i];
        char *id_i = ids_tmp[i];
        float row_i[N_PARAM_COLS];

        memcpy(row_i, vals_tmp[i], sizeof(row_i));
        k = i - 1;
        while (k >= 0 && cats[k] > cat_i) {
            cats[k + 1] = cats[k];
            ids_tmp[k + 1] = ids_tmp[k];
            memcpy(vals_tmp[k + 1], vals_tmp[k], sizeof(row_i));
            k--;
        }
        cats[k + 1] = cat_i;
        ids_tmp[k + 1] = id_i;
        memcpy(vals_tmp[k + 1], row_i, sizeof(row_i));
    }

    *btot_out = n_rows;
    *f7p_out = af2d(n_rows, N_PARAM_COLS);
    *basin_ids_out = (char **)G_malloc(n_rows * sizeof(char *));
    for (i = 0; i < n_rows; i++) {
        (*basin_ids_out)[i] = ids_tmp[i];
        for (j = 0; j < N_PARAM_COLS; j++)
            (*f7p_out)[i][j] = vals_tmp[i][j];
    }

    /* group rows sharing the same basin_id into one reporting station
     * (several HRU rows of one basin all resolve to the same station);
     * for hru_bands == 1 (one row per basin) this yields rtot == btot
     * and station_of[b] == b for every b, i.e. today's behavior. */
    {
        char *station_ids_tmp[4096];
        int n_stations = 0;
        int *station_of = (int *)G_malloc(n_rows * sizeof(int));

        for (i = 0; i < n_rows; i++) {
            int found = -1;

            for (k = 0; k < n_stations; k++)
                if (strcmp(station_ids_tmp[k], ids_tmp[i]) == 0) {
                    found = k;
                    break;
                }
            if (found < 0) {
                found = n_stations;
                station_ids_tmp[n_stations++] = ids_tmp[i];
            }
            station_of[i] = found;
        }

        *rtot_out = n_stations;
        *station_of_out = station_of;
        *station_ids_out = (char **)G_malloc(n_stations * sizeof(char *));
        for (k = 0; k < n_stations; k++)
            (*station_ids_out)[k] = station_ids_tmp[k];
    }
}

void hbv_basins_from_vector(
    const char *elevation, const char *outlets, const char *id_column,
    const char *threshold, const char *snap_radius, const char *landcover,
    const char *forest_cats, const char *ffo_default,
    const char *ffi_default, const char *precip_station_elevation,
    const char *temp_station_elevation, const char *et_station_elevation,
    const char *parameters_template, const char *hru_bands,
    const char *basins_raster, const char *basins_vector, int *btot_out,
    float ***f7p_out, char ***basin_ids_out, int *rtot_out,
    int **station_of_out, char ***station_ids_out) {
    if (elevation && outlets)
        delineate_with_basins_module(
            elevation, outlets, id_column, threshold, snap_radius,
            landcover, forest_cats, ffo_default, ffi_default,
            precip_station_elevation, temp_station_elevation,
            et_station_elevation, parameters_template, hru_bands,
            basins_raster, basins_vector);

    read_basins_table(basins_vector, btot_out, f7p_out, basin_ids_out,
                       rtot_out, station_of_out, station_ids_out);
}
