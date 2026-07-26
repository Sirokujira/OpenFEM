/*
input_data.c

.ofe 入力ファイルの読み込み。

書式は OpenFDTD の .ofd に準拠する (1 行 1 キー、"key = value..."、
先頭行がプログラム名 + バージョン、"end" で終端、'#' 以降はコメント)。
mesh / material / geometry の書式は .ofd からそのまま流用できる。
*/

#include "fem.h"
#include "fem_prototype.h"

#define MAXTOKEN (1000)
#define ARRAY_INC (1000)

// 空白区切りでトークンに分解する ("str" は破壊される)
static int tokenize(char *str, const char *sep, char *token[], int maxtoken)
{
	int n = 0;

	char *p = strtok(str, sep);
	while ((p != NULL) && (n < maxtoken)) {
		token[n++] = p;
		p = strtok(NULL, sep);
	}

	return n;
}

// 領域境界と分割数から節点座標を作る
static double *mesh_expand(int nreg, const double *pos, const int *ndiv, int *ncell)
{
	int n = 0;
	for (int i = 0; i < nreg; i++) {
		n += ndiv[i];
	}
	if (n <= 0) return NULL;

	double *p = (double *)malloc((n + 1) * sizeof(double));
	int m = 0;
	p[0] = pos[0];
	for (int i = 0; i < nreg; i++) {
		const double d = (pos[i + 1] - pos[i]) / ndiv[i];
		for (int j = 1; j <= ndiv[i]; j++) {
			p[++m] = pos[i] + (d * j);
		}
	}
	*ncell = n;

	return p;
}


int input_data(FILE *fp)
{
	int    nline = 0;
	int    nxr = 0, nyr = 0, nzr = 0;
	int    *dxr = NULL, *dyr = NULL, *dzr = NULL;
	double *xr = NULL, *yr = NULL, *zr = NULL;
	char   prog[BUFSIZ];
	char   strline[BUFSIZ], strkey[BUFSIZ], strsave[BUFSIZ];
	char   *token[MAXTOKEN];
	const char sep[] = " \t\n\r";
	const char errfmt1[] = "*** too many %s data\n";
	const char errfmt2[] = "*** invalid %s data\n";
	const char errfmt3[] = "*** invalid %s data #%d\n";

	// 既定値 (キー省略時の動作)

	strcpy(Title, "");
	strcpy(prog, "");

	NMaterial = 2;		// 0 : 真空, 1 : PEC (予約)
	Material = (material_t *)malloc(NMaterial * sizeof(material_t));
	for (int m = 0; m < NMaterial; m++) {
		Material[m].epsr  = 1;
		Material[m].sigma = 0;
		Material[m].mur   = 1;
		Material[m].tand  = 0;
		Material[m].debye = 0;
		Material[m].einf  = 1;
		Material[m].deps  = 0;
		Material[m].tau   = 0;
	}

	NGeometry = 0;
	Geometry = NULL;

	NConductor = 0;
	Conductor = NULL;
	NPort = 0;

	Analysis = ANALYSIS_C;
	Tline = 0;
	LineLength = 0;
	Volt = 1;
	Freq = 0;
	Curr = 1;
	NSection = 1;
	for (int p = 0; p < MAXPORT; p++) {
		CondSigma[p] = 0;
	}

	Solver.maxiter = 10000;
	Solver.nout = 100;
	Solver.converg = 1e-9;

	// 読み込み

	while (fgets(strline, sizeof(strline), fp) != NULL) {
		// コメント
		char *pcomment = strchr(strline, '#');
		if (pcomment != NULL) *pcomment = '\0';

		// 空行
		if (strspn(strline, " \t\n\r") == strlen(strline)) continue;

		// 終端
		if (!strncmp(strline, "end", 3)) break;

		strcpy(strsave, strline);

		const int ntoken = tokenize(strline, sep, token, MAXTOKEN);
		if (ntoken < 1) continue;

		// ヘッダー以外は "key = value..." の形であること
		if ((nline > 0) && ((ntoken < 3) || strcmp(token[1], "="))) continue;

		strcpy(strkey, token[0]);
		const int nval = ntoken - 2;		// 値の個数

		if      (nline == 0) {
			strcpy(prog, strkey);
			if (strcmp(prog, PROGRAM)) {
				printf("*** not %s data\n", PROGRAM);
				return 1;
			}
			nline++;
		}
		else if (!strcmp(strkey, "title")) {
			char *peq = strchr(strsave, '=');
			if (peq != NULL) {
				strcpy(Title, peq + 1);
				// 前後の空白を落とす
				char *p = Title;
				while ((*p == ' ') || (*p == '\t')) p++;
				memmove(Title, p, strlen(p) + 1);
				size_t len = strlen(Title);
				while ((len > 0) && (strchr(" \t\n\r", Title[len - 1]) != NULL)) {
					Title[--len] = '\0';
				}
			}
		}
		else if (!strcmp(strkey, "xmesh") || !strcmp(strkey, "ymesh") || !strcmp(strkey, "zmesh")) {
			if ((nval < 3) || (nval % 2 == 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int nreg = (nval - 1) / 2;
			double *pos = (double *)malloc((nreg + 1) * sizeof(double));
			int   *ndiv = (int *)malloc(nreg * sizeof(int));
			pos[0] = atof(token[2]);
			for (int i = 0; i < nreg; i++) {
				ndiv[i] = atoi(token[(2 * i) + 3]);
				pos[i + 1] = atof(token[(2 * i) + 4]);
				if ((ndiv[i] <= 0) || (pos[i + 1] <= pos[i])) {
					printf(errfmt2, strkey);
					free(pos);
					free(ndiv);
					return 1;
				}
			}
			if      (strkey[0] == 'x') { nxr = nreg; xr = pos; dxr = ndiv; }
			else if (strkey[0] == 'y') { nyr = nreg; yr = pos; dyr = ndiv; }
			else                       { nzr = nreg; zr = pos; dzr = ndiv; }
		}
		else if (!strcmp(strkey, "material")) {
			// material = <epsr> <sigma>
			// (.ofd 互換 : material = [type] <epsr> <esgm> <amur> <msgm> も受け付け、
			//  透磁率は準静的解析では使わないため無視する)
			double epsr = 1, sigma = 0;
			if      (nval == 2) {
				epsr  = atof(token[2]);
				sigma = atof(token[3]);
			}
			else if (nval == 4) {
				epsr  = atof(token[2]);
				sigma = atof(token[3]);
			}
			else if (nval == 5) {
				epsr  = atof(token[3]);
				sigma = atof(token[4]);
			}
			else {
				printf(errfmt3, strkey, NMaterial - 1);
				return 1;
			}
			if ((epsr <= 0) || (sigma < 0)) {
				printf(errfmt3, strkey, NMaterial - 1);
				return 1;
			}
			if (NMaterial % ARRAY_INC == 2) {
				Material = (material_t *)realloc(Material, (NMaterial + ARRAY_INC) * sizeof(material_t));
			}
			if (NMaterial >= UCHAR_MAX) {
				printf(errfmt1, strkey);
				return 1;
			}
			Material[NMaterial].epsr  = epsr;
			Material[NMaterial].sigma = sigma;
			Material[NMaterial].mur   = 1;
			Material[NMaterial].tand  = 0;
			Material[NMaterial].debye = 0;
			Material[NMaterial].einf  = 1;
			Material[NMaterial].deps  = 0;
			Material[NMaterial].tau   = 0;
			NMaterial++;
		}
		else if (!strcmp(strkey, "geometry")) {
			// geometry = <material_id> <shape> <g0> ... <g5>
			if (nval < 8) {
				printf(errfmt3, strkey, NGeometry + 1);
				return 1;
			}
			if (NGeometry % ARRAY_INC == 0) {
				Geometry = (geometry_t *)realloc(Geometry, (NGeometry + ARRAY_INC) * sizeof(geometry_t));
			}
			Geometry[NGeometry].m     = atoi(token[2]);
			Geometry[NGeometry].shape = atoi(token[3]);
			for (int n = 0; n < 6; n++) {
				Geometry[NGeometry].g[n] = atof(token[4 + n]);
			}
			Geometry[NGeometry].g[6] = Geometry[NGeometry].g[7] = 0;
			NGeometry++;
		}
		else if (!strcmp(strkey, "conductor")) {
			// conductor = <id> <shape> <g0> ... <g5>   (id = 0 : 基準導体)
			if (nval < 8) {
				printf(errfmt3, strkey, NConductor + 1);
				return 1;
			}
			// id = -1 は「導体を取り消す」指定 (先に置いた導体をくり抜く)
			const int id = atoi(token[2]);
			if ((id < -1) || (id >= MAXPORT)) {
				printf(errfmt3, strkey, NConductor + 1);
				return 1;
			}
			if (NConductor % ARRAY_INC == 0) {
				Conductor = (conductor_t *)realloc(Conductor, (NConductor + ARRAY_INC) * sizeof(conductor_t));
			}
			Conductor[NConductor].id    = id;
			Conductor[NConductor].shape = atoi(token[3]);
			for (int n = 0; n < 6; n++) {
				Conductor[NConductor].g[n] = atof(token[4 + n]);
			}
			Conductor[NConductor].g[6] = Conductor[NConductor].g[7] = 0;
			if (id > NPort) NPort = id;
			NConductor++;
		}
		else if (!strcmp(strkey, "mur") || !strcmp(strkey, "tand")) {
			// mur  = <material_id> <mur>    : 比透磁率 (静磁場解析)
			// tand = <material_id> <tand>   : 誘電正接 (frequency 指定時の誘電損)
			if (nval < 2) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int mid = atoi(token[2]);
			const double val = atof(token[3]);
			if ((mid < 0) || (mid >= NMaterial) || (val < 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			if (strkey[0] == 'm') {
				if (val <= 0) {
					printf(errfmt2, strkey);
					return 1;
				}
				Material[mid].mur = val;
			}
			else {
				Material[mid].tand = val;
			}
		}
		else if (!strcmp(strkey, "debye")) {
			// debye = <material_id> <eps_inf> <delta_eps> <tau [s]>
			// εr(ω) = eps_inf + delta_eps / (1 + jωτ)。frequency 指定時に
			// epsr (実部) と tand (= εr''/εr') に展開する
			if (nval < 4) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int mid = atoi(token[2]);
			if ((mid < 2) || (mid >= NMaterial)) {
				printf(errfmt2, strkey);
				return 1;
			}
			const double einf = atof(token[3]);
			const double deps = atof(token[4]);
			const double tau  = atof(token[5]);
			if ((einf <= 0) || (deps < 0) || (tau < 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			Material[mid].debye = 1;
			Material[mid].einf  = einf;
			Material[mid].deps  = deps;
			Material[mid].tau   = tau;
		}
		else if (!strcmp(strkey, "conductorsigma")) {
			// conductorsigma = <conductor_id> <sigma>  : 導体の DC 直列抵抗の計算に使う
			if (nval < 2) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int cid = atoi(token[2]);
			const double val = atof(token[3]);
			if ((cid < 0) || (cid >= MAXPORT) || (val <= 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			CondSigma[cid] = val;
		}
		else if (!strcmp(strkey, "frequency")) {
			// frequency = <f [Hz]>  : tanδ による並列コンダクタンスの計算に使う
			Freq = atof(token[2]);
			if (Freq < 0) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "current")) {
			// current = <I [A]>  : 静磁場解析の励振電流
			Curr = atof(token[2]);
			if (Curr == 0) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "analysis")) {
			// analysis = C L R (部分集合)
			Analysis = 0;
			for (int n = 0; n < nval; n++) {
				const int c = toupper((int)token[2 + n][0]);
				if      (c == 'C') Analysis |= ANALYSIS_C;
				else if (c == 'L') Analysis |= ANALYSIS_L;
				else if (c == 'R') Analysis |= ANALYSIS_R;
				else if (c == 'M') Analysis |= ANALYSIS_M;
				else if (c == 'F') Analysis |= ANALYSIS_F;
				else {
					printf(errfmt2, strkey);
					return 1;
				}
			}
			// L は真空静電界を使うので C の求解と同じ道具立てが要る
			if (Analysis == 0) Analysis = ANALYSIS_C;
		}
		else if (!strcmp(strkey, "tline")) {
			const int c = toupper((int)token[2][0]);
			if ((c != 'X') && (c != 'Y') && (c != 'Z')) {
				printf(errfmt2, strkey);
				return 1;
			}
			Tline = (char)c;
		}
		else if (!strcmp(strkey, "linelength")) {
			// 等価回路 (ofe_post) を作るときの線路長 [m]。省略時は解析領域長
			LineLength = atof(token[2]);
			if (LineLength <= 0) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "voltage")) {
			Volt = atof(token[2]);
			if (Volt == 0) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "nsection")) {
			NSection = atoi(token[2]);
			if (NSection < 1) NSection = 1;
		}
		else if (!strcmp(strkey, "solver")) {
			if (nval < 3) {
				printf(errfmt2, strkey);
				return 1;
			}
			Solver.maxiter = atoi(token[2]);
			Solver.nout    = atoi(token[3]);
			Solver.converg = atof(token[4]);
			if (Solver.maxiter < 1) Solver.maxiter = 1;
			if (Solver.nout    < 1) Solver.nout = 1;
			if (Solver.converg <= 0) Solver.converg = 1e-9;
		}
		else {
			// 未知のキーは無視する (ポスト処理用キーの共存を許す)
			;
		}
	}

	// Debye 分散材料 を frequency の値で展開する
	// (frequency 未指定なら静的な εr(0) = einf + deps、損失なしとして扱う)

	for (int m = 0; m < NMaterial; m++) {
		if (!Material[m].debye) continue;
		const double wt = 2 * PI * Freq * Material[m].tau;
		const double den = 1 + (wt * wt);
		const double er = Material[m].einf + (Material[m].deps / den);
		const double ei = Material[m].deps * wt / den;
		Material[m].epsr = er;
		Material[m].tand = ei / er;
	}

	// 格子を展開する

	if ((nxr == 0) || (nyr == 0) || (nzr == 0)) {
		printf("%s\n", "*** no mesh data");
		return 1;
	}
	Xn = mesh_expand(nxr, xr, dxr, &Nx);
	Yn = mesh_expand(nyr, yr, dyr, &Ny);
	Zn = mesh_expand(nzr, zr, dzr, &Nz);
	free(xr); free(yr); free(zr);
	free(dxr); free(dyr); free(dzr);
	if ((Xn == NULL) || (Yn == NULL) || (Zn == NULL)) {
		printf("%s\n", "*** invalid mesh data");
		return 1;
	}

	// 整合性の確認

	for (int n = 0; n < NGeometry; n++) {
		if ((Geometry[n].m < 0) || (Geometry[n].m >= NMaterial)) {
			printf("*** invalid geometry material id #%d\n", n + 1);
			return 1;
		}
	}

	if (NConductor < 1) {
		printf("%s\n", "*** no conductor data");
		return 1;
	}
	if (NPort < 1) {
		printf("%s\n", "*** no port conductor (conductor id >= 1) data");
		return 1;
	}
	// 基準導体 (id = 0) が要る
	int have_ref = 0;
	for (int n = 0; n < NConductor; n++) {
		if (Conductor[n].id == 0) have_ref = 1;
	}
	// (実際に節点が残っているかは setup() で確認する)
	if (!have_ref) {
		printf("%s\n", "*** no reference conductor (conductor id = 0) data");
		return 1;
	}
	// ポート番号は 1..NPort が連続していること
	for (int p = 1; p <= NPort; p++) {
		int found = 0;
		for (int n = 0; n < NConductor; n++) {
			if (Conductor[n].id == p) found = 1;
		}
		if (!found) {
			printf("*** conductor id %d is missing\n", p);
			return 1;
		}
	}

	// L 解析は TEM 仮定 (単位長あたり) でのみ意味を持つ
	if ((Analysis & ANALYSIS_L) && !Tline) {
		printf("%s\n", "*** analysis L requires the tline key (TEM per-unit-length)");
		return 1;
	}
	// M 解析 (静磁場) は断面 2 次元の定式化なので伝送線路軸が要る
	if ((Analysis & ANALYSIS_M) && !Tline) {
		printf("%s\n", "*** analysis M requires the tline key (2D magnetostatic)");
		return 1;
	}
	// F 解析 (渦電流) も断面 2 次元。周波数と導体の導電率が要る
	if (Analysis & ANALYSIS_F) {
		if (!Tline) {
			printf("%s\n", "*** analysis F requires the tline key (2D eddy current)");
			return 1;
		}
		if (Freq <= 0) {
			printf("%s\n", "*** analysis F requires the frequency key");
			return 1;
		}
	}

	return 0;
}
