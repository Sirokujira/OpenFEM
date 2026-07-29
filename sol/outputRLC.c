/*
outputRLC.c

抽出した回路パラメータのログ出力と、ポスト処理用バイナリ (ofe.out) の出力。
*/

#include "fem.h"
#include "fem_prototype.h"

#define OUT_MAGIC "OFEOUT04"

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
		print_matrix(fp, "Inductance matrix L (TEM, external only)", "H/m", Lmat, np);
	}

	if (HavePfe) {
		// 鉄損は B の非線形関数なので重ね合わせが効かない。対角 [k][k] が
		// 「ポート k だけを current で励振したときの損失」で、非対角は定義できない
		fprintf(fp, "\n%s [%s]:\n", "Iron loss Pfe (Bertotti, per driven port)",
			(pul ? "W/m" : "W"));
		for (int i = 0; i < np; i++) {
			fprintf(fp, "  port %d = %14.6e\n", i + 1, Pfemat[(i * np) + i]);
		}
	}

	if (HaveM) {
		print_matrix(fp, "Inductance matrix L (magnetostatic, incl. internal)", "H/m", Mmat, np);
	}

	if (HaveR) {
		if (Freq > 0) {
			fprintf(fp, "\n(shunt loss evaluated at %.6e [Hz]; tan-delta enters as "
				"sigma_d = omega*eps0*epsr*tand)\n", Freq);
		}
		print_matrix(fp, "Shunt conductance matrix G", (pul ? "S/m" : "S"), Gmat, np);
		print_matrix(fp, "Shunt resistance matrix R = inv(G)", (pul ? "ohm*m" : "ohm"), Rmat, np);
	}

	if (HaveS) {
		print_matrix(fp, "Series resistance matrix Rs (DC, conductor loss)", "ohm/m", Smat, np);
	}

	if (HaveF) {
		fprintf(fp, "\n(eddy current solution at %.6e [Hz]; includes skin and "
			"proximity effect)\n", Freq);
		print_matrix(fp, "Series resistance matrix R(f)", (pul ? "ohm/m" : "ohm"), Rfmat, np);
		print_matrix(fp, "Series inductance matrix L(f)", (pul ? "H/m" : "H"), Lfmat, np);
		if (HaveS && (np >= 1) && (Smat[0] > 0)) {
			fprintf(fp, "  R(f)/R(dc) = %.4f (port 1)\n", Rfmat[0] / Smat[0]);
		}
	}

	// 単一ポートの伝送線路定数
	// Z0 は高周波量なので、外部インダクタンス (TEM) があればそちらを使う
	// (静磁場の L_dc は表皮効果の無い低周波モデルの値)
	if (pul && HaveC && (HaveL || HaveM) && (np == 1)) {
		const double c = Cmat[0];
		const double l = (HaveL ? Lmat[0] : Mmat[0]);
		if ((c > 0) && (l > 0)) {
			const double z0 = sqrt(l / c);
			const double vp = 1 / sqrt(l * c);
			const double eeff = (C0 * C0) * l * c;
			fprintf(fp, "\nTransmission line (L = %s):\n",
				(HaveL ? "TEM, external" : "magnetostatic, incl. internal"));
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

	const int32_t flag[8] = {HaveC, HaveL, HaveR, NSection, HaveM, HaveS, HaveF, HavePfe};
	const int32_t tline = (int32_t)Tline;
	const double dval[3] = {LineLength, Volt, Freq};

	fwrite(OUT_MAGIC, 1, 8, fp);
	fwrite(&np, sizeof(int32_t), 1, fp);
	fwrite(flag, sizeof(int32_t), 8, fp);
	fwrite(&tline, sizeof(int32_t), 1, fp);
	fwrite(dval, sizeof(double), 3, fp);
	fwrite(title, 1, sizeof(title), fp);

	const size_t nn = (size_t)np * np;
	fwrite(Cmat, sizeof(double), nn, fp);
	fwrite(Lmat, sizeof(double), nn, fp);
	fwrite(Gmat, sizeof(double), nn, fp);
	fwrite(Rmat, sizeof(double), nn, fp);
	fwrite(Mmat, sizeof(double), nn, fp);
	fwrite(Smat, sizeof(double), nn, fp);
	fwrite(Rfmat, sizeof(double), nn, fp);
	fwrite(Lfmat, sizeof(double), nn, fp);
	fwrite(Pfemat, sizeof(double), nn, fp);
}
