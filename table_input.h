/* Reads a long-format (id,date,value) GRASS DB table -- schema
 * overridable via id_column/date_column/value_column, default
 * "station_id"/"date"/"value" -- and pivots it into a
 * [n_stations][n_dates] float array, station order matching the given
 * station_ids array.
 *
 * *dates and *n_dates are an in-out canonical date axis: on the first call
 * (with *dates == NULL) they are populated from this table's own
 * distinct, sorted date values; every subsequent call (temperature
 * after precipitation, etc.) must resolve to exactly that same date
 * set or this calls G_fatal_error() -- kept strict rather than
 * interpolating/filling gaps, so a misaligned input fails loudly
 * instead of silently corrupting the run. Every station in
 * station_ids must have exactly one value per canonical date; any
 * missing cell, or any station_id in the table not present in
 * station_ids, is also a fatal error. Caller owns *dates once set
 * (each string G_store()'d). */
float **table_read_pivot(const char *table, const char *id_column,
			  const char *date_column, const char *value_column,
			  char **station_ids, int n_stations, char ***dates,
			  int *n_dates);
