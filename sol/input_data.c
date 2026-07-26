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
		Material[m].npole = 0;
		Material[m].einf  = 1;
		Material[m].bhaniso = 0;
		Material[m].ja.on = 0;
		for (int d = 0; d < 3; d++) Material[m].nbh[d] = 0;
		for (int d = 0; d < 6; d++) {
			Material[m].eps6[d] = 0;	// 0 : anisoeps 未指定 (最後に epsr で埋める)
			Material[m].mu6[d]  = 0;
		}
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

	NSweep = 0;
	JaSub = 20;

	NlMaxiter = 50;
	NlTol = 1e-5;
	NlRelax = 1.0;

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
			Material[NMaterial].npole = 0;
			Material[NMaterial].einf  = 1;
			Material[NMaterial].bhaniso = 0;
			Material[NMaterial].ja.on = 0;
			for (int d = 0; d < 3; d++) Material[NMaterial].nbh[d] = 0;
			for (int d = 0; d < 6; d++) {
				Material[NMaterial].eps6[d] = 0;
				Material[NMaterial].mu6[d]  = 0;
			}
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
		else if (!strcmp(strkey, "debye") || !strcmp(strkey, "lorentz")) {
			// debye   = <material_id> <eps_inf> (<deps> <tau>)...
			// lorentz = <material_id> <eps_inf> (<deps> <f0> <delta>)...
			//   εr(ω) = eps_inf + Σ Debye Δε/(1+jωτ)
			//                   + Σ Lorentz Δε ω0^2/(ω0^2-ω^2+jωδ)
			// frequency の値で epsr (実部) と tand (= εr''/εr') に展開する。
			// 同じ材料に複数行書くと極が追加される (eps_inf は最後の指定が効く)。
			const int lor = !strcmp(strkey, "lorentz");
			const int nper = (lor ? 3 : 2);
			if ((nval < 2 + nper) || (((nval - 2) % nper) != 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int mid = atoi(token[2]);
			if ((mid < 2) || (mid >= NMaterial)) {
				printf(errfmt2, strkey);
				return 1;
			}
			material_t *mt = &Material[mid];
			const double einf = atof(token[3]);
			if (einf <= 0) {
				printf(errfmt2, strkey);
				return 1;
			}
			mt->einf = einf;
			const int npair = (nval - 2) / nper;
			for (int q = 0; q < npair; q++) {
				if (mt->npole >= MAXPOLE) {
					printf(errfmt1, strkey);
					return 1;
				}
				const char **tk = (const char **)&token[4 + (q * nper)];
				pole_t *pl = &mt->pole[mt->npole];
				pl->type = (lor ? 2 : 1);
				pl->a = atof(tk[0]);
				pl->b = atof(tk[1]);
				pl->c = (lor ? atof(tk[2]) : 0);
				if ((pl->a < 0) || (pl->b <= 0) || (lor && (pl->c < 0))) {
					printf(errfmt2, strkey);
					return 1;
				}
				mt->npole++;
			}
		}
		else if (!strcmp(strkey, "anisoeps") || !strcmp(strkey, "anisomur")) {
			// anisoeps = <material_id> <exx> <eyy> <ezz> [<exy> <eyz> <ezx>]
			// anisomur = <material_id> <mxx> <myy> <mzz> [<mxy> <myz> <mzx>]
			//   3 値なら対角テンソル、6 値なら非対角を含む対称テンソル
			//   (主軸が格子軸と一致しない材料は回転させた 6 成分を与える)
			if ((nval != 4) && (nval != 7)) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int mid = atoi(token[2]);
			if ((mid < 0) || (mid >= NMaterial)) {
				printf(errfmt2, strkey);
				return 1;
			}
			double v[6] = {0, 0, 0, 0, 0, 0};
			for (int d = 0; d < nval - 1; d++) {
				v[d] = atof(token[3 + d]);
			}
			// 対角成分は正、テンソルは正定値であること
			if ((v[0] <= 0) || (v[1] <= 0) || (v[2] <= 0)) {
				printf(errfmt2, strkey);
				return 1;
			}
			double chk[6];
			if (tensor6_inverse(v, chk)) {
				printf("*** %s : the tensor is singular\n", strkey);
				return 1;
			}
			double *dst = ((strkey[5] == 'e') ? Material[mid].eps6 : Material[mid].mu6);
			for (int d = 0; d < 6; d++) dst[d] = v[d];
		}
		else if (!strcmp(strkey, "bh")) {
			// bh = <material_id> <H [A/m]> <B [T]>            等方性 (|B| に対する曲線)
			// bh = <material_id> <X|Y|Z> <H [A/m]> <B [T]>    軸毎 (直交異方性)
			// B は狭義単調増加、B > 0。原点は書かない (B < B1 は初期透磁率で扱う)
			if (nval < 3) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int mid = atoi(token[2]);
			if ((mid < 2) || (mid >= NMaterial)) {
				printf(errfmt2, strkey);
				return 1;
			}
			material_t *mt = &Material[mid];

			// 軸指定の有無を判定する (数字で始まらないトークンなら軸)
			int ax0 = 0, ax1 = 2, iv = 3;
			const int c0 = toupper((int)token[3][0]);
			if ((c0 == 'X') || (c0 == 'Y') || (c0 == 'Z')) {
				if (nval < 4) {
					printf(errfmt2, strkey);
					return 1;
				}
				ax0 = ax1 = ((c0 == 'X') ? 0 : (c0 == 'Y') ? 1 : 2);
				iv = 4;
				mt->bhaniso = 1;
			}

			const double hh = atof(token[iv]);
			const double bb = atof(token[iv + 1]);
			if ((hh <= 0) || (bb <= 0)) {
				printf("*** bh : H and B must be positive (origin is implicit)\n");
				return 1;
			}
			for (int ax = ax0; ax <= ax1; ax++) {
				const int nb = mt->nbh[ax];
				if (nb >= MAXBH) {
					printf(errfmt1, strkey);
					return 1;
				}
				if ((nb > 0) && ((bb <= mt->bh_b[ax][nb - 1]) || (hh <= mt->bh_h[ax][nb - 1]))) {
					printf("*** bh : H and B must increase monotonically (material %d)\n", mid);
					return 1;
				}
				mt->bh_h[ax][nb] = hh;
				mt->bh_b[ax][nb] = bb;
				mt->nbh[ax]++;
			}
		}
		else if (!strcmp(strkey, "ja")) {
			// ja = <material_id> <Ms> <a> <alpha> <k> <c>
			//   Jiles-Atherton ヒステリシスモデル。currentsweep と併用する
			if (nval < 6) {
				printf(errfmt2, strkey);
				return 1;
			}
			const int mid = atoi(token[2]);
			if ((mid < 2) || (mid >= NMaterial)) {
				printf(errfmt2, strkey);
				return 1;
			}
			ja_t *ja = &Material[mid].ja;
			ja->ms    = atof(token[3]);
			ja->a     = atof(token[4]);
			ja->alpha = atof(token[5]);
			ja->k     = atof(token[6]);
			ja->c     = atof(token[7]);
			if ((ja->ms <= 0) || (ja->a <= 0) || (ja->k <= 0)
			 || (ja->alpha < 0) || (ja->c < 0) || (ja->c >= 1)) {
				printf(errfmt2, strkey);
				return 1;
			}
			ja->on = 1;
		}
		else if (!strcmp(strkey, "currentsweep")) {
			// currentsweep = <I1> <I2> ...  : 履歴に沿って順に解く
			if (nval < 1) {
				printf(errfmt2, strkey);
				return 1;
			}
			if (nval > MAXSWEEP) {
				printf(errfmt1, strkey);
				return 1;
			}
			NSweep = nval;
			for (int q = 0; q < nval; q++) {
				Sweep[q] = atof(token[2 + q]);
			}
		}
		else if (!strcmp(strkey, "jasub")) {
			// jasub = <n>  : 1 ステップあたりの J-A 部分積分数
			JaSub = atoi(token[2]);
			if (JaSub < 1) JaSub = 1;
		}
		else if (!strcmp(strkey, "nlsolver")) {
			// nlsolver = <maxiter> <tol> <relax>  : 非線形 (B-H) 反復の設定
			if (nval < 3) {
				printf(errfmt2, strkey);
				return 1;
			}
			NlMaxiter = atoi(token[2]);
			NlTol     = atof(token[3]);
			NlRelax   = atof(token[4]);
			if (NlMaxiter < 1) NlMaxiter = 1;
			if (NlTol <= 0) NlTol = 1e-5;
			if ((NlRelax <= 0) || (NlRelax > 1)) NlRelax = 1.0;
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

	// 分散材料を frequency の値で展開する
	// (frequency 未指定なら静的な εr(0) = ε∞ + ΣΔε、損失なしとして扱う)

	const double omega_d = 2 * PI * Freq;
	for (int m = 0; m < NMaterial; m++) {
		material_t *mt = &Material[m];
		if (mt->npole <= 0) continue;
		double er = mt->einf, ei = 0;
		for (int q = 0; q < mt->npole; q++) {
			const pole_t *pl = &mt->pole[q];
			if (pl->type == 1) {
				// Debye : Δε / (1 + jωτ)
				const double wt = omega_d * pl->b;
				const double den = 1 + (wt * wt);
				er += pl->a / den;
				ei += pl->a * wt / den;
			}
			else {
				// Lorentz : Δε ω0^2 / (ω0^2 - ω^2 + jωδ)
				const double w0 = 2 * PI * pl->b;
				const double dl = 2 * PI * pl->c;
				const double d1 = (w0 * w0) - (omega_d * omega_d);
				const double d2 = omega_d * dl;
				const double den = (d1 * d1) + (d2 * d2);
				if (den <= 0) continue;
				er += pl->a * w0 * w0 * d1 / den;
				ei += pl->a * w0 * w0 * d2 / den;
			}
		}
		if (er <= 0) {
			printf("*** dispersive material %d has a non-positive epsr at %.4e Hz\n", m, Freq);
			return 1;
		}
		mt->epsr = er;
		mt->tand = ei / er;
	}

	// 異方性が指定されていない材料は等方性 (epsr / mur) で埋める

	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].eps6[0] <= 0) {
			Material[m].eps6[0] = Material[m].eps6[1] = Material[m].eps6[2] = Material[m].epsr;
			Material[m].eps6[3] = Material[m].eps6[4] = Material[m].eps6[5] = 0;
		}
		if (Material[m].mu6[0] <= 0) {
			Material[m].mu6[0] = Material[m].mu6[1] = Material[m].mu6[2] = Material[m].mur;
			Material[m].mu6[3] = Material[m].mu6[4] = Material[m].mu6[5] = 0;
		}
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
	// 非線形磁性体 (B-H) は重ね合わせが成り立たないので単一ポートに限る
	int nonlinear = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].nbh[0] > 0) nonlinear = 1;
		// 軸毎に与えるときは 3 軸とも必要
		if (Material[m].bhaniso) {
			for (int d = 0; d < 3; d++) {
				if (Material[m].nbh[d] <= 0) {
					printf("*** bh : material %d needs a curve for all three axes "
						"when axes are specified\n", m);
					return 1;
				}
			}
		}
	}
	// ヒステリシス (Jiles-Atherton)
	int hysteresis = 0;
	for (int m = 0; m < NMaterial; m++) {
		if (Material[m].ja.on) hysteresis = 1;
		if (Material[m].ja.on && (Material[m].nbh[0] > 0)) {
			printf("*** material %d cannot have both bh and ja\n", m);
			return 1;
		}
	}
	if (hysteresis) {
		if (!(Analysis & ANALYSIS_M)) {
			printf("%s\n", "*** the ja key needs analysis M (magnetostatic)");
			return 1;
		}
		if (Analysis & ANALYSIS_F) {
			printf("%s\n", "*** the ja key cannot be combined with analysis F");
			return 1;
		}
		if (NPort != 1) {
			printf("%s\n", "*** a hysteretic (ja) model requires exactly one port");
			return 1;
		}
		if (NSweep < 1) {
			// 掃引が無いときは current の 1 点だけを処女曲線として解く
			NSweep = 1;
			Sweep[0] = Curr;
		}
	}
	else if (NSweep > 0) {
		printf("%s\n", "*** currentsweep needs a hysteretic (ja) material");
		return 1;
	}

	if (nonlinear) {
		if (!(Analysis & ANALYSIS_M)) {
			printf("%s\n", "*** the bh key needs analysis M (magnetostatic)");
			return 1;
		}
		if (Analysis & ANALYSIS_F) {
			printf("%s\n", "*** the bh key cannot be combined with analysis F "
				"(time-harmonic eddy current assumes a linear material)");
			return 1;
		}
		if (NPort != 1) {
			printf("%s\n", "*** a nonlinear (bh) model requires exactly one port "
				"(superposition does not hold)");
			return 1;
		}
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
