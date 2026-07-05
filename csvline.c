#include "csvline.h"

int csv_split_line(char *line, char **fields, int max_fields) {
	int n = 0;
	char *p = line;

	while (*p && n < max_fields) {
		char *out = p;
		char *start = p;

		if (*p == '"') {
			p++;
			start = p;
			out = p;
			while (*p) {
				if (*p == '"' && *(p + 1) == '"') {
					*out++ = '"';
					p += 2;
				} else if (*p == '"') {
					p++;
					break;
				} else {
					*out++ = *p++;
				}
			}
			*out = '\0';
		} else {
			while (*p && *p != ',')
				p++;
		}
		fields[n++] = start;
		if (*p == ',')
			p++;
		else
			break;
	}
	return n;
}
