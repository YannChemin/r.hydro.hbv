/* splits one CSV line (as produced by GRASS modules' format=csv output:
 * double-quoted text fields, doubled quotes as escape) into fields, in
 * place. Returns the number of fields found (capped at max_fields). */
int csv_split_line(char *line, char **fields, int max_fields);
