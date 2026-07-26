/*
OpenFEM post processor (ofe_post)

ofe.out を読み、回路パラメータを CSV と SPICE サブサーキットに書き出す。
ソルバー本体のグローバルには依存しない (単体で動く)。

出力:
	rlc.csv       行列 (C, L, G, R) の一覧
	ofe_circuit.sp  等価回路 (.SUBCKT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define FN_out   "ofe.out"
#define FN_csv   "rlc.csv"
#define FN_spice "ofe_circuit.sp"
#define OUT_MAGIC "OFEOUT02"

typedef struct {
	int32_t np;
	int32_t haveC, haveL, haveR, nsection, haveM, haveS;
	int32_t tline;
	double  length, volt, freq;
	char    title[256];
	double  *c, *l, *g, *r, *m, *s;
	const double *lind;			// 等価回路に使う L (静磁場があればそちら)
} result_t;


static int readout(const char *fname, result_t *res)
{
	char magic[8];

	FILE *fp = fopen(fname, "rb");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fname);
		return 1;
	}

	int ierr = 0;
	int32_t flag[6];
	double dval[3];
	if ((fread(magic, 1, 8, fp) != 8) || memcmp(magic, OUT_MAGIC, 8)) {
		printf("*** %s is not an %s output file.\n", fname, "OpenFEM");
		ierr = 1;
	}
	if (!ierr) {
		if ((fread(&res->np, sizeof(int32_t), 1, fp) != 1) ||
		    (fread(flag, sizeof(int32_t), 6, fp) != 6) ||
		    (fread(&res->tline, sizeof(int32_t), 1, fp) != 1) ||
		    (fread(dval, sizeof(double), 3, fp) != 3) ||
		    (fread(res->title, 1, sizeof(res->title), fp) != sizeof(res->title))) {
			printf("*** %s is broken.\n", fname);
			ierr = 1;
		}
	}
	if (!ierr) {
		res->haveC   = flag[0];
		res->haveL   = flag[1];
		res->haveR   = flag[2];
		res->nsection = ((flag[3] > 0) ? flag[3] : 1);
		res->haveM   = flag[4];
		res->haveS   = flag[5];
		res->length  = dval[0];
		res->volt    = dval[1];
		res->freq    = dval[2];
		res->title[sizeof(res->title) - 1] = '\0';

		const size_t nn = (size_t)res->np * res->np;
		if ((res->np < 1) || (nn > 1000000)) {
			printf("*** invalid port count in %s.\n", fname);
			ierr = 1;
		}
		else {
			res->c = (double *)malloc(nn * sizeof(double));
			res->l = (double *)malloc(nn * sizeof(double));
			res->g = (double *)malloc(nn * sizeof(double));
			res->r = (double *)malloc(nn * sizeof(double));
			res->m = (double *)malloc(nn * sizeof(double));
			res->s = (double *)malloc(nn * sizeof(double));
			if ((fread(res->c, sizeof(double), nn, fp) != nn) ||
			    (fread(res->l, sizeof(double), nn, fp) != nn) ||
			    (fread(res->g, sizeof(double), nn, fp) != nn) ||
			    (fread(res->r, sizeof(double), nn, fp) != nn) ||
			    (fread(res->m, sizeof(double), nn, fp) != nn) ||
			    (fread(res->s, sizeof(double), nn, fp) != nn)) {
				printf("*** %s is broken (matrix).\n", fname);
				ierr = 1;
			}
			// 高周波モデルを既定とし、外部インダクタンス (TEM) を優先して使う。
			// 静磁場の L_dc は analysis から L を外したときに使われる
			res->lind = (res->haveL ? res->l : res->m);
		}
	}

	fclose(fp);

	return ierr;
}


static void write_matrix_csv(FILE *fp, const char *name, const char *unit,
	const double *m, int np)
{
	fprintf(fp, "%s,%s\n", name, unit);
	fprintf(fp, "i\\j");
	for (int j = 0; j < np; j++) {
		fprintf(fp, ",%d", j + 1);
	}
	fprintf(fp, "\n");
	for (int i = 0; i < np; i++) {
		fprintf(fp, "%d", i + 1);
		for (int j = 0; j < np; j++) {
			fprintf(fp, ",%.8e", m[(i * np) + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}


static void write_csv(const result_t *res)
{
	FILE *fp = fopen(FN_csv, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", FN_csv);
		return;
	}

	const int pul = (res->tline != 0);
	const int np = res->np;

	fprintf(fp, "title,%s\n", res->title);
	fprintf(fp, "ports,%d\n", np);
	fprintf(fp, "per unit length,%s\n", (pul ? "yes" : "no"));
	if (pul) fprintf(fp, "line length (equivalent circuit) [m],%.8e\n", res->length);
	fprintf(fp, "\n");

	if (res->haveC) {
		write_matrix_csv(fp, "C", (pul ? "F/m" : "F"), res->c, np);
		fprintf(fp, "SPICE capacitance,%s\n", (pul ? "F/m" : "F"));
		for (int i = 0; i < np; i++) {
			double sum = 0;
			for (int j = 0; j < np; j++) sum += res->c[(i * np) + j];
			fprintf(fp, "C(%d;gnd),%.8e\n", i + 1, sum);
		}
		for (int i = 0; i < np; i++) {
			for (int j = i + 1; j < np; j++) {
				fprintf(fp, "C(%d;%d),%.8e\n", i + 1, j + 1, -res->c[(i * np) + j]);
			}
		}
		fprintf(fp, "\n");
	}
	if (res->haveL) write_matrix_csv(fp, "L", "H/m", res->l, np);
	if (res->haveM) write_matrix_csv(fp, "Ldc", "H/m", res->m, np);
	if (res->haveR) {
		if (res->freq > 0) fprintf(fp, "frequency [Hz],%.8e\n\n", res->freq);
		write_matrix_csv(fp, "G", (pul ? "S/m" : "S"), res->g, np);
		write_matrix_csv(fp, "R", (pul ? "ohm*m" : "ohm"), res->r, np);
	}
	if (res->haveS) write_matrix_csv(fp, "Rs", "ohm/m", res->s, np);

	if (res->haveC && (res->haveL || res->haveM) && (np == 1)
	 && (res->c[0] > 0) && (res->lind[0] > 0)) {
		const double c = 2.99792458e8;
		fprintf(fp, "Z0 [ohm],%.8e\n", sqrt(res->lind[0] / res->c[0]));
		fprintf(fp, "eps_eff,%.8e\n", c * c * res->lind[0] * res->c[0]);
		fprintf(fp, "delay [s/m],%.8e\n", sqrt(res->lind[0] * res->c[0]));
	}

	fclose(fp);
	printf("output : %s\n", FN_csv);
}


// 等価回路 (SPICE サブサーキット)
static void write_spice(const result_t *res)
{
	FILE *fp = fopen(FN_spice, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", FN_spice);
		return;
	}

	const int np = res->np;
	const int pul = (res->tline != 0);
	const int ns = (pul ? res->nsection : 1);
	// 単位長あたりの値を段あたりの値に直す係数
	const double len = (pul ? (res->length / ns) : 1.0);

	fprintf(fp, "* %s : equivalent circuit generated by OpenFEM\n", res->title);
	fprintf(fp, "* ports = %d, %s\n", np, (pul ? "distributed (RLGC ladder)" : "lumped"));
	if (pul && (res->haveL || res->haveM)) {
		fprintf(fp, "* L = %s\n",
			(res->haveL ? "TEM (external only)" : "magnetostatic (incl. internal)"));
	}
	if (res->haveS) fprintf(fp, "* series R is the DC conductor resistance (no skin effect)\n");
	if (res->freq > 0) fprintf(fp, "* shunt G evaluated at %.6e [Hz]\n", res->freq);
	if (pul) fprintf(fp, "* line length = %.6e [m], sections = %d\n", res->length, ns);
	fprintf(fp, "*\n");

	fprintf(fp, ".SUBCKT OFE_CIRCUIT");
	for (int i = 0; i < np; i++) {
		fprintf(fp, " IN%d", i + 1);
	}
	if (pul) {
		for (int i = 0; i < np; i++) {
			fprintf(fp, " OUT%d", i + 1);
		}
	}
	fprintf(fp, " GND\n");

	if (!pul) {
		// 集中定数 : ポート - GND 間、ポート間の C と R
		for (int i = 0; i < np; i++) {
			double sum = 0;
			for (int j = 0; j < np; j++) sum += res->c[(i * np) + j];
			if (res->haveC && (sum > 0)) {
				fprintf(fp, "C%d_G IN%d GND %.6e\n", i + 1, i + 1, sum);
			}
		}
		for (int i = 0; i < np; i++) {
			for (int j = i + 1; j < np; j++) {
				const double cm = -res->c[(i * np) + j];
				if (res->haveC && (cm > 0)) {
					fprintf(fp, "C%d_%d IN%d IN%d %.6e\n", i + 1, j + 1, i + 1, j + 1, cm);
				}
			}
		}
		if (res->haveR) {
			for (int i = 0; i < np; i++) {
				double sum = 0;
				for (int j = 0; j < np; j++) sum += res->g[(i * np) + j];
				if (sum > 0) {
					fprintf(fp, "R%d_G IN%d GND %.6e\n", i + 1, i + 1, 1 / sum);
				}
			}
			for (int i = 0; i < np; i++) {
				for (int j = i + 1; j < np; j++) {
					const double gm = -res->g[(i * np) + j];
					if (gm > 0) {
						fprintf(fp, "R%d_%d IN%d IN%d %.6e\n", i + 1, j + 1, i + 1, j + 1, 1 / gm);
					}
				}
			}
		}
	}
	else {
		// 分布定数 : ns 段の RLGC 梯子 (直列 L、並列 C/G)
		for (int s = 0; s < ns; s++) {
			fprintf(fp, "* --- section %d ---\n", s + 1);
			for (int i = 0; i < np; i++) {
				// 節点名 : N<port>_<section>  (0 段目の入口は IN<port>)
				char na[64], nb[64];
				if (s == 0) sprintf(na, "IN%d", i + 1);
				else        sprintf(na, "N%d_%d", i + 1, s);
				if (s == ns - 1) sprintf(nb, "OUT%d", i + 1);
				else             sprintf(nb, "N%d_%d", i + 1, s + 1);

				// 直列 R (導体の DC 抵抗) -> 直列 L の順に並べる
				char nmid[64];
				if (res->haveS && (res->s[(i * np) + i] > 0)) {
					sprintf(nmid, "NR%d_%d", i + 1, s + 1);
					fprintf(fp, "R%d_%d %s %s %.6e\n", i + 1, s + 1, na, nmid,
						res->s[(i * np) + i] * len);
				}
				else {
					snprintf(nmid, sizeof(nmid), "%s", na);
				}

				if ((res->haveL || res->haveM) && (res->lind[(i * np) + i] > 0)) {
					fprintf(fp, "L%d_%d %s %s %.6e\n", i + 1, s + 1, nmid, nb,
						res->lind[(i * np) + i] * len);
				}
				else if (strcmp(nmid, nb)) {
					fprintf(fp, "V%d_%d %s %s 0\n", i + 1, s + 1, nmid, nb);
				}

				double sum = 0;
				for (int j = 0; j < np; j++) sum += res->c[(i * np) + j];
				if (res->haveC && (sum > 0)) {
					fprintf(fp, "C%d_%d %s GND %.6e\n", i + 1, s + 1, nb, sum * len);
				}
				if (res->haveR) {
					double gsum = 0;
					for (int j = 0; j < np; j++) gsum += res->g[(i * np) + j];
					if (gsum > 0) {
						// 並列コンダクタンス G は抵抗 1/G で表す
						fprintf(fp, "RG%d_%d %s GND %.6e\n", i + 1, s + 1, nb, 1 / (gsum * len));
					}
				}
			}
			// 線路間の結合
			for (int i = 0; i < np; i++) {
				for (int j = i + 1; j < np; j++) {
					const double cm = -res->c[(i * np) + j];
					if (res->haveC && (cm > 0)) {
						char nb1[64], nb2[64];
						if (s == ns - 1) { sprintf(nb1, "OUT%d", i + 1); sprintf(nb2, "OUT%d", j + 1); }
						else             { sprintf(nb1, "N%d_%d", i + 1, s + 1); sprintf(nb2, "N%d_%d", j + 1, s + 1); }
						fprintf(fp, "CM%d%d_%d %s %s %.6e\n", i + 1, j + 1, s + 1, nb1, nb2, cm * len);
					}
					if (res->haveL || res->haveM) {
						const double lii = res->lind[(i * np) + i];
						const double ljj = res->lind[(j * np) + j];
						const double lij = res->lind[(i * np) + j];
						if ((lii > 0) && (ljj > 0)) {
							const double kc = lij / sqrt(lii * ljj);
							if (fabs(kc) > 1e-12) {
								fprintf(fp, "K%d%d_%d L%d_%d L%d_%d %.6f\n",
									i + 1, j + 1, s + 1, i + 1, s + 1, j + 1, s + 1, kc);
							}
						}
					}
				}
			}
		}
	}

	fprintf(fp, ".ENDS OFE_CIRCUIT\n");

	fclose(fp);
	printf("output : %s\n", FN_spice);
}


int main(int argc, char *argv[])
{
	result_t res;
	memset(&res, 0, sizeof(res));

	const char *fname = FN_out;
	for (int n = 1; n < argc; n++) {
		if (!strcmp(argv[n], "-n")) {
			n++;			// スレッド指定は互換のため受け付けて無視する
		}
		else if (strstr(argv[n], ".out") != NULL) {
			fname = argv[n];
		}
	}

	if (readout(fname, &res)) return 1;

	write_csv(&res);
	write_spice(&res);

	free(res.c);
	free(res.l);
	free(res.g);
	free(res.r);
	free(res.m);
	free(res.s);

	return 0;
}
