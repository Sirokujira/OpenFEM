// fem.h
//
// OpenFEM : 準静的 FEM による回路パラメータ (R/L/C) 抽出ソルバー
// 共有ヘッダ (格子・材料・導体・行列のグローバル)
//
// OpenFDTD (https://github.com/Sirokujira/OpenFDTD) の構成・入力書式に
// 準拠する (mesh/material/geometry の書式は共通)。

#ifndef _FEM_H_
#define _FEM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "OpenFEM"
#define VERSION_MAJOR (1)
#define VERSION_MINOR (0)
#define VERSION_BUILD (0)

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

#define FN_log   "ofe.log"
#define FN_out   "ofe.out"
#define FN_csv   "rlc.csv"
#define FN_spice "ofe_circuit.sp"

// 数学・物理定数 (OpenFDTD の ofd.h と同じ定義)
#define PI     (4.0 * atan(1.0))
#define C0     (2.99792458e8)
#define MU0    (4 * PI * 1e-7)
#define EPS0   (1 / (C0 * C0 * MU0))
#define ETA0   (C0 * MU0)
#define EPS    (1e-6)

// 解析種別ビット (analysis キー)
#define ANALYSIS_C (1 << 0)		// 静電界 -> 容量行列 C
#define ANALYSIS_L (1 << 1)		// 真空静電界 -> インダクタンス行列 L (TEM 近似)
#define ANALYSIS_R (1 << 2)		// 定常電流界 -> 抵抗/コンダクタンス行列 R, G
#define ANALYSIS_M (1 << 3)		// 静磁場 -> DC インダクタンス行列 (内部インダクタンス込み)
#define ANALYSIS_F (1 << 4)		// 時間調和渦電流 -> 表皮効果を含む R(f), L(f)

#define MAXPORT (16)			// 導体 (基準導体 + ポート) の最大数
#define MAXBH   (64)			// B-H 曲線の点数の最大数
#define MAXPOLE (8)				// 分散材料の極の最大数
#define MAXSWEEP (64)			// 電流掃引の点数の最大数

// Jiles-Atherton ヒステリシスモデルのパラメータ
//   M_an(He) = Ms (coth(He/a) - a/He),  He = H + α M
//   dM/dH = [ (1-c) δ_M (M_an-M)/(δk - α(M_an-M)) + c dM_an/dHe ]
//         / [ 1 - α c dM_an/dHe ]
typedef struct {
	int    on;					// 1 : このモデルを使う
	double ms;					// 飽和磁化 Ms [A/m]
	double a;					// 形状パラメータ a [A/m]
	double alpha;				// 磁区間結合 α
	double k;					// ピン止め k [A/m]
	double c;					// 可逆成分の割合 c
} ja_t;

// 分散材料の極
//   type 1 (Debye)   : Δε / (1 + jωτ)                     a=Δε, b=τ [s]
//   type 2 (Lorentz) : Δε ω0^2 / (ω0^2 - ω^2 + jωδ)       a=Δε, b=f0 [Hz], c=δ [Hz]
typedef struct {
	int    type;
	double a, b, c;
} pole_t;

// 材料 (id=0 : 真空, id=1 : PEC 予約, id>=2 : ユーザー定義)
typedef struct {
	double epsr;				// 比誘電率
	double sigma;				// 導電率 [S/m]
	double mur;					// 比透磁率 (静磁場解析で使う、既定 1)
	double tand;				// 誘電正接 tanδ (frequency 指定時に G へ寄与、既定 0)
	int    npole;				// 分散極の数 (0 : 分散なし)
	double einf;				// ε∞
	pole_t pole[MAXPOLE];
	// 異方性テンソル (対称、成分順 xx, yy, zz, xy, yz, zx)。既定は等方性
	double eps6[6];				// 比誘電率テンソル
	double mu6[6];				// 比透磁率テンソル
	// B-H 曲線 (軸毎)。bhaniso = 0 なら軸 0 の曲線を |B| に対して等方的に使う
	int    bhaniso;				// 1 : 軸毎に別の曲線 (直交異方性)
	int    nbh[3];				// 各軸の点数 (0 : 線形、mu6 を使う)
	double bh_h[3][MAXBH];		// H [A/m]
	double bh_b[3][MAXBH];		// B [T] (狭義単調増加、B > 0)
	ja_t   ja;					// ヒステリシス (Jiles-Atherton)
} material_t;

// 形状 (shape/g は OpenFDTD の geometry と共通)
typedef struct {
	int    m;					// 材料番号
	int    shape;				// 形状番号 (1:直方体, 2:楕円体, 11/12/13:円柱)
	double g[8];
} geometry_t;

// 導体 (電極)。id=0 が基準導体 (グランド)、id>=1 がポート
typedef struct {
	int    id;					// 導体番号
	int    shape;
	double g[8];
} conductor_t;

typedef struct {
	int    maxiter;				// 最大反復回数
	int    nout;				// 収束履歴の出力間隔
	double converg;				// 収束判定値 (相対残差)
} solver_t;

// CRS 疎行列 (対称だが上下三角とも保持する)
typedef struct {
	int64_t n;					// 次元 (節点数)
	int64_t nnz;				// 非ゼロ要素数
	int64_t *rowptr;			// [n+1]
	int32_t *col;				// [nnz]
	double  *val;				// [nnz]
} crs_t;

// ---- グローバル ----

EXTERN char Title[BUFSIZ];

// 格子 (節点座標)。セル数 Nx x Ny x Nz、節点数 (Nx+1) x (Ny+1) x (Nz+1)
EXTERN int Nx, Ny, Nz;
EXTERN double *Xn, *Yn, *Zn;

// セル中心座標 (形状判定用)
EXTERN double *Xc, *Yc, *Zc;

EXTERN int NMaterial;
EXTERN material_t *Material;

EXTERN int NGeometry;
EXTERN geometry_t *Geometry;

EXTERN int NConductor;			// conductor 行の数
EXTERN conductor_t *Conductor;

EXTERN int NPort;				// ポート数 (= 導体番号の最大値、基準導体 0 を除く)

EXTERN solver_t Solver;

EXTERN int Analysis;			// ANALYSIS_* のビット OR
EXTERN char Tline;				// 'X'/'Y'/'Z' : 単位長あたりで出力、0 : 絶対値
EXTERN double LineLength;		// 等価回路を作るときの線路長 [m] (既定 : 解析領域長)
EXTERN double Volt;				// 励振電圧 [V]
EXTERN double Freq;				// 周波数 [Hz] (tanδ による誘電損の計算に使う、0 = 無効)
EXTERN double Curr;				// 静磁場解析の励振電流 [A]
EXTERN int NSection;			// SPICE 等価回路の梯子段数

// 電流掃引 (ヒステリシス解析)。履歴に沿って順に解く
EXTERN int NSweep;
EXTERN double Sweep[MAXSWEEP];

// Jiles-Atherton の状態 (Gauss 点毎、添字 = cell * 8 + g)
EXTERN double *JaB, *JaH, *JaM;		// 収束済みの |B|, H, M
EXTERN double *JaBn, *JaHn, *JaMn;	// 反復中の |B|, H, M
EXTERN double *JaD, *JaDn;		// 磁化方向の基準ベクトル (符号付き B を取るため)
EXTERN int JaSub;				// 1 ステップあたりの部分積分数

// 非線形 (B-H) 反復の設定
EXTERN int NlMaxiter;			// 最大反復回数
EXTERN double NlTol;			// 収束判定 (ν の最大相対変化)
EXTERN double NlRelax;			// 緩和係数 (0 < w <= 1)

// セル毎の材料番号
EXTERN unsigned char *CellMaterial;

// 節点毎の導体番号 (-1 : 導体でない、0.. : 導体番号)
EXTERN signed char *NodeConductor;

// セル毎の導体番号 (断面積・電流密度の計算に使う)
EXTERN signed char *CellConductor;

// 導体毎の量 ([0] = 基準導体、[1..NPort] = ポート)
EXTERN double CondSigma[MAXPORT];	// 導電率 [S/m] (0 = 未指定)
EXTERN double CondArea[MAXPORT];	// 断面積 [m^2] (tline 指定時)

// 抽出結果 (行列は [NPort][NPort]、単位は Tline 指定時 F/m, H/m, S/m)
EXTERN double *Cmat;			// 短絡容量行列 (Maxwell 行列)
EXTERN double *Lmat;			// インダクタンス行列 (TEM 近似)
EXTERN double *Gmat;			// コンダクタンス行列
EXTERN double *Rmat;			// 並列抵抗行列 (= inv(G))
EXTERN double *Mmat;			// DC インダクタンス行列 (静磁場、内部インダクタンス込み)
EXTERN double *Smat;			// 直列抵抗行列 (導体の DC 抵抗) [ohm/m]
EXTERN double *Rfmat;			// 直列抵抗行列 R(f) (渦電流、表皮効果込み) [ohm/m]
EXTERN double *Lfmat;			// 直列インダクタンス行列 L(f) (渦電流) [H/m]
EXTERN int HaveC, HaveL, HaveR, HaveM, HaveS, HaveF;	// 各行列が計算済みか
EXTERN double TlineLength;		// 伝送線路長 [m] (Tline 指定時)

#ifdef __cplusplus
}
#endif

#endif		// _FEM_H_
