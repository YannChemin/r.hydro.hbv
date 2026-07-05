#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/spawn.h>

#include "arrays.h"
#include "csvline.h"
#include "table_input.h"

struct raw_row {
	char *id;
	char *date;
	float value;
};

static int strcmp_cb(const void *a, const void *b) {
	return strcmp(*(char *const *)a, *(char *const *)b);
}

static int strcmp_ptr_cb(const void *a, const void *b) {
	const char *const *pa = a;
	const char *const *pb = b;
	return strcmp(*pa, *pb);
}

/* Spawns db.select and returns the parsed (unpivoted) rows; *n_rows_out
 * is set to the row count. Caller owns the returned array and its
 * per-row id/date strings (G_store()'d). */
static struct raw_row *spawn_and_read(const char *table, const char *id_column,
				      const char *date_column,
				      const char *value_column,
				      int *n_rows_out) {
	char sql[2048];
	char *tmpfile;
	const char *argv[16];
	int argc = 0;
	char sql_arg[2048 + 8], output_arg[GPATH_MAX + 16];
	FILE *f;
	char buf[1024];
	struct raw_row *rows = NULL;
	int n_rows = 0, cap = 0;

	snprintf(sql, sizeof(sql), "SELECT %s,%s,%s FROM %s ORDER BY %s,%s",
		 id_column, date_column, value_column, table, date_column,
		 id_column);

	tmpfile = G_tempfile();

	argv[argc++] = "db.select";
	snprintf(sql_arg, sizeof(sql_arg), "sql=%s", sql);
	argv[argc++] = sql_arg;
	argv[argc++] = "format=csv";
	snprintf(output_arg, sizeof(output_arg), "output=%s", tmpfile);
	argv[argc++] = output_arg;
	argv[argc++] = NULL;

	if (G_vspawn_ex(argv[0], argv) != 0) {
		remove(tmpfile);
		G_fatal_error(_("Unable to read table <%s> (missing columns "
				"%s/%s/%s?)"),
			      table, id_column, date_column, value_column);
	}

	f = fopen(tmpfile, "r");
	if (!f)
		G_fatal_error(_("Unable to open temporary file <%s>"), tmpfile);

	if (!G_getl(buf, sizeof(buf), f)) /* header */
		G_fatal_error(_("Table <%s> is empty"), table);

	while (G_getl(buf, sizeof(buf), f)) {
		char *fields[3];

		if (csv_split_line(buf, fields, 3) != 3)
			G_fatal_error(_("Malformed row in table <%s>"), table);

		if (n_rows >= cap) {
			cap = cap ? cap * 2 : 4096;
			rows = G_realloc(rows, cap * sizeof(struct raw_row));
		}
		rows[n_rows].id = G_store(fields[0]);
		rows[n_rows].date = G_store(fields[1]);
		rows[n_rows].value = (float)atof(fields[2]);
		n_rows++;
	}
	fclose(f);
	remove(tmpfile);

	*n_rows_out = n_rows;
	return rows;
}

static char **distinct_sorted(struct raw_row *rows, int n_rows,
			      int is_date, int *n_out) {
	char **all = G_malloc(n_rows * sizeof(char *));
	char **uniq;
	int i, n_uniq;

	for (i = 0; i < n_rows; i++)
		all[i] = is_date ? rows[i].date : rows[i].id;

	qsort(all, n_rows, sizeof(char *), strcmp_cb);

	uniq = G_malloc(n_rows * sizeof(char *));
	n_uniq = 0;
	for (i = 0; i < n_rows; i++) {
		if (i == 0 || strcmp(all[i], all[i - 1]) != 0)
			uniq[n_uniq++] = G_store(all[i]);
	}
	G_free(all);

	*n_out = n_uniq;
	return uniq;
}

float **table_read_pivot(const char *table, const char *id_column,
			  const char *date_column, const char *value_column,
			  char **station_ids, int n_stations, char ***dates,
			  int *n_dates) {
	struct raw_row *rows;
	int n_rows, i, j;
	char **table_dates;
	int n_table_dates;
	float **result;
	char *filled;

	rows = spawn_and_read(table, id_column, date_column, value_column,
			       &n_rows);

	table_dates = distinct_sorted(rows, n_rows, 1, &n_table_dates);

	if (*dates == NULL) {
		*dates = table_dates;
		*n_dates = n_table_dates;
	} else {
		if (n_table_dates != *n_dates)
			G_fatal_error(
			    _("Table <%s> has %d distinct dates, expected "
			      "%d (does not match the date range of the "
			      "first time-series table read) -- all "
			      "time-series tables must share exactly the "
			      "same set of dates"),
			    table, n_table_dates, *n_dates);
		for (i = 0; i < n_table_dates; i++)
			if (strcmp(table_dates[i], (*dates)[i]) != 0)
				G_fatal_error(
				    _("Table <%s>'s dates do not match the "
				      "first time-series table read (first "
				      "mismatch: <%s> vs <%s>)"),
				    table, table_dates[i], (*dates)[i]);
	}

	result = af2d(n_stations, *n_dates);
	filled = G_calloc(n_stations * (size_t)*n_dates, sizeof(char));

	for (i = 0; i < n_rows; i++) {
		int station_idx = -1, date_idx;
		char **found;

		for (j = 0; j < n_stations; j++)
			if (strcmp(station_ids[j], rows[i].id) == 0) {
				station_idx = j;
				break;
			}
		if (station_idx < 0)
			G_fatal_error(
			    _("Table <%s> contains station id <%s>, which "
			      "is not one of this run's basin IDs"),
			    table, rows[i].id);

		found = bsearch(&rows[i].date, *dates, *n_dates,
				 sizeof(char *), strcmp_ptr_cb);
		if (!found)
			G_fatal_error(_("Internal error: date <%s> not "
					"found in canonical date axis"),
				      rows[i].date);
		date_idx = (int)(found - *dates);

		result[station_idx][date_idx] = rows[i].value;
		filled[station_idx * (size_t)*n_dates + date_idx] = 1;
	}

	for (i = 0; i < n_stations; i++)
		for (j = 0; j < *n_dates; j++)
			if (!filled[i * (size_t)*n_dates + j])
				G_fatal_error(
				    _("Table <%s> is missing a value for "
				      "station <%s> on date <%s>"),
				    table, station_ids[i], (*dates)[j]);

	G_free(filled);
	return result;
}
