#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arrays.h"

static FILE *open_report_file(const char *outdir, const char *prefix, int i,
			       const char *basin_id) {
	char path[4096];
	snprintf(path, sizeof(path), "%s/%s%02d-%s.csv", outdir, prefix,
		 i + 1, basin_id);
	FILE *f = fopen(path, "w");
	if (f == NULL) {
		fprintf(stderr, "hbv_report: cannot open '%s'\n", path);
		exit(EXIT_FAILURE);
	}
	return f;
}

static FILE *open_long_csv(const char *path, const char *header,
			    const char *csvt) {
	char csvt_path[4096];
	const char *ext;
	size_t stem_len;
	FILE *f, *ft;

	if (!path)
		return NULL;
	f = fopen(path, "w");
	if (f == NULL) {
		fprintf(stderr, "hbv_report: cannot open '%s'\n", path);
		exit(EXIT_FAILURE);
	}
	fprintf(f, "%s\n", header);

	/* OGR's CSV driver expects the sidecar type file to share the data
	 * file's basename with a .csvt extension (e.g. foo.csv ->
	 * foo.csvt), not foo.csv.csvt. */
	ext = strrchr(path, '.');
	stem_len = (ext && strcmp(ext, ".csv") == 0) ? (size_t)(ext - path)
						      : strlen(path);
	snprintf(csvt_path, sizeof(csvt_path), "%.*s.csvt", (int)stem_len,
		 path);
	ft = fopen(csvt_path, "w");
	if (ft == NULL) {
		fprintf(stderr, "hbv_report: cannot open '%s'\n", csvt_path);
		exit(EXIT_FAILURE);
	}
	fprintf(ft, "%s\n", csvt);
	fclose(ft);

	return f;
}

/* Generalized report writer: one basin loop instead of hardcoded per-basin
 * fopen() calls, output directory and basin IDs supplied by the caller.
 * table_basinout_csv/table_etout_csv, if non-NULL, additionally get a
 * long-format (station_id,day,...) CSV (+ a paired .csvt sidecar giving
 * OGR explicit column types) written across all basins, for the caller
 * to import into GRASS DB tables with db.in.ogr -- an alternative,
 * GRASS-native view of the same Basinout/ETout data. */
int hbv_report(const char *outdir, char **basin_ids, int btot, float *mtot,
	       int ttot, float **qo, float **qr, float **po, float **eta,
	       float **eto, float **mfc, float **mbeta, float **mlp,
	       float **malpha, float **mkf, float **mks, float **mperc,
	       float **mcflux, float **MNSET, float **MRVEET, float **MNS,
	       float **MNSH, float **MNSL, float **MRVE, float **MRMAEH,
	       float **MRMAEL, float **MREVE, const char *table_basinout_csv,
	       const char *table_etout_csv) {
	int i, t, q;
	FILE *tb = open_long_csv(table_basinout_csv, "station_id,day,qo,qr,po",
				  "String,Integer,Real,Real,Real");
	FILE *te = open_long_csv(table_etout_csv, "station_id,day,eta,eto",
				  "String,Integer,Real,Real");

	for (i = 0; i < btot; i++) {
		FILE *f1 = open_report_file(outdir, "Output", i, basin_ids[i]);
		FILE *f2 = open_report_file(outdir, "Basinout", i, basin_ids[i]);
		FILE *f3 = open_report_file(outdir, "ETout", i, basin_ids[i]);
		FILE *f4 = open_report_file(outdir, "EToutput", i, basin_ids[i]);

		for (t = 0; t < ttot; t++) {
			fprintf(f2, "%f %f %f \n", qo[i][t], qr[i][t], po[i][t]);
			fprintf(f3, "%f %f \n", eta[i][t], eto[i][t]);
			if (tb)
				fprintf(tb, "%s,%d,%f,%f,%f\n", basin_ids[i], t,
					qo[i][t], qr[i][t], po[i][t]);
			if (te)
				fprintf(te, "%s,%d,%f,%f\n", basin_ids[i], t,
					eta[i][t], eto[i][t]);
		}
		fprintf(f1, "fc(mm),beta(-),lp(-),alpha(-),kf(1/d),ks(1/d),perc(mm/d),cflux(mm/d),NS(-),NSH(-),NSL(-),RVE(pc),RMAEH(-),RMAEL(-),REVE(pc)\n");
		fprintf(f4, "fc(mm),beta(-),lp(-),alpha(-),kf(1/d),ks(1/d),perc(mm/d),cflux(mm/d),NSET(-),RVEET(pc)\n");
		for (q = 0; q < (int)mtot[i]; q++) {
			fprintf(f1, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
				mfc[i][q], mbeta[i][q], mlp[i][q], malpha[i][q],
				mkf[i][q], mks[i][q], mperc[i][q], mcflux[i][q],
				MNS[i][q], MNSH[i][q], MNSL[i][q], MRVE[i][q],
				MRMAEH[i][q], MRMAEL[i][q], MREVE[i][q]);
			fprintf(f4, "%f %f %f %f %f %f %f %f %f %f\n",
				mfc[i][q], mbeta[i][q], mlp[i][q], malpha[i][q],
				mkf[i][q], mks[i][q], mperc[i][q], mcflux[i][q],
				MNSET[i][q], MRVEET[i][q]);
		}

		fclose(f1);
		fclose(f2);
		fclose(f3);
		fclose(f4);
	}
	if (tb)
		fclose(tb);
	if (te)
		fclose(te);
	return (1);
}
