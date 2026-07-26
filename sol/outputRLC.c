/*
outputRLC.c

抽出した回路パラメータのログ出力と、ポスト処理用バイナリ (ofe.out) の出力。
*/

#include "fem.h"
#include "fem_prototype.h"

#define OUT_MAGIC "OFEOUT01"

static void print_matrix(FILE *fp, const char *name, const char *unit, const double *m, int np)
{
	fprintf(fp, "\n%s [%s]:\n", name, unit);
	for (int i = 0; i < np; i++) {
		fprintf(fp, " ");
		for (int j = 0; j < np; j++) {
			fprintf(fp, " %14.6e", m[(i * np) + j]);
		}
		fprintf(fp, "\n");
	}
}


void outputRLC(FILE *fp)
{
	const int np = NPort;
	const int pul = (Tline && (TlineLength > 0));		// 単位長あたりか

	fprintf(fp, "\n=== circuit parameters ===\n");
	if (pul) {
		fprintf(fp, "per unit length along %c (model length = %.6e [m])\n", Tline, TlineLength);
		fprintf(fp, "equivalent circuit line length = %.6e [m], sections = %d\n", LineLength, NSection);
	}
	else {
		fprintf(fp, "lumped (total) values\n");
	}

	if (HaveC) {
		print_matrix(fp, "Capacitance matrix C", (pul ? "F/m" : "F"), Cmat, np);

		// 等価回路 (SPICE) 形式の容量
		fprintf(fp, "\nSPICE capacitances [%s]:\n", (pul ? "F/m" : "F"));
		for (int i = 0; i < np; i++) {
			double sum = 0;
			for (int j = 0; j < np; j++) {
				sum += Cmat[(i * np) + j];
			}
			fprintf(fp, "  C(%d,gnd) = %14.6e\n", i + 1, sum);
		}
		for (int i = 0; i < np; i++) {
			for (int j = i + 1; j < np; j++) {
				fprintf(fp, "  C(%d,%d)   = %14.6e\n", i + 1, j + 1, -Cmat[(i * np) + j]);
			}
		}
	}

	if (HaveL) {
		print_matrix(fp, "Inductance matrix L (TEM)", "H/m", Lmat, np);
	}

	if (HaveR) {
		print_matrix(fp, "Conductance matrix G", (pul ? "S/m" : "S"), Gmat, np);
		print_matrix(fp, "Resistance matrix R = inv(G)", (pul ? "ohm*m" : "ohm"), Rmat, np);
	}

	// 単一ポートの伝送線路定数
	if (pul && HaveC && HaveL && (np == 1)) {
		const double c = Cmat[0];
		const double l = Lmat[0];
		if ((c > 0) && (l > 0)) {
			const double z0 = sqrt(l / c);
			const double vp = 1 / sqrt(l * c);
			const double eeff = (C0 * C0) * l * c;
			fprintf(fp, "\nTransmission line (TEM):\n");
			fprintf(fp, "  Z0        = %14.6e [ohm]\n", z0);
			fprintf(fp, "  eps_eff   = %14.6e\n", eeff);
			fprintf(fp, "  v_p       = %14.6e [m/s] (%.4f c)\n", vp, vp / C0);
			fprintf(fp, "  delay     = %14.6e [s/m]\n", 1 / vp);
		}
	}

	fflush(fp);
}


void writeout(FILE *fp)
{
	const int32_t np = NPort;
	char title[256];

	memset(title, 0, sizeof(title));
	const size_t tlen = strlen(Title);
	memcpy(title, Title, ((tlen < sizeof(title) - 1) ? tlen : sizeof(title) - 1));

	const int32_t flag[4] = {HaveC, HaveL, HaveR, NSection};
	const int32_t tline = (int32_t)Tline;
	const double dval[2] = {LineLength, Volt};

	fwrite(OUT_MAGIC, 1, 8, fp);
	fwrite(&np, sizeof(int32_t), 1, fp);
	fwrite(flag, sizeof(int32_t), 4, fp);
	fwrite(&tline, sizeof(int32_t), 1, fp);
	fwrite(dval, sizeof(double), 2, fp);
	fwrite(title, 1, sizeof(title), fp);

	const size_t nn = (size_t)np * np;
	fwrite(Cmat, sizeof(double), nn, fp);
	fwrite(Lmat, sizeof(double), nn, fp);
	fwrite(Gmat, sizeof(double), nn, fp);
	fwrite(Rmat, sizeof(double), nn, fp);
}
