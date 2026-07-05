#include <stdio.h>
#include <stdlib.h>
#include "arrays.h"

/* http://cboard.cprogramming.com/c-programming/47105-how-read-csv-file.html
 * NOTE: row is the true column count of one line in the file; the
 * previous ARRAYSIZE(*array) == sizeof(float*) trick only worked by
 * coincidence when row happened to equal 8 (pointer size on 64-bit).
 */
float ** readcsv(char *filename, int row, int col){
	int x,y;
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		fprintf(stderr, "readcsv: cannot open '%s'\n", filename);
		exit(EXIT_FAILURE);
	}
	float **array	= af2d(col,row);//array[row][col] Transposed for reading
	size_t i, j;
	char buffer[BUFSIZ], *ptr;
	/* array only has col rows allocated (af2d(col,row) above): the
	 * outer loop must stop there even if the file itself has more
	 * lines (e.g. a full multi-year CSV read with a shorter dtot/col
	 * to match a shorter co-input, such as an ERA5-fetched date range
	 * narrower than a bundled historical observation file) -- reading
	 * every fgets()-able line regardless of col wrote past the
	 * allocation and segfaulted. */
	for ( i = 0; i < (size_t)col && fgets(buffer, sizeof buffer, file); ++i ){
		for ( j = 0, ptr = buffer; j < (size_t)row; ++j, ++ptr ){
			array[i][j] = (float)strtof(ptr, &ptr);
		}
	}
	fclose(file);
	float **arrayout= af2d(row,col);//array[row][col]
	for (x=0;x<row;x++)
		for (y=0;y<col;y++)
			arrayout[x][y]=array[y][x];
	return &arrayout[0];
}
