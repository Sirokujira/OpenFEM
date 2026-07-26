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

#define MAXPORT (16)			// 導体 (基準導体 + ポート) の最大数

// 材料 (id=0 : 真空, id=1 : PEC 予約, id>=2 : ユーザー定義)
typedef struct {
	double epsr;				// 比誘電率
	double sigma;				// 導電率 [S/m]
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
EXTERN int NSection;			// SPICE 等価回路の梯子段数

// セル毎の材料番号
EXTERN unsigned char *CellMaterial;

// 節点毎の導体番号 (-1 : 導体でない、0.. : 導体番号)
EXTERN signed char *NodeConductor;

// 抽出結果 (行列は [NPort][NPort]、単位は Tline 指定時 F/m, H/m, S/m)
EXTERN double *Cmat;			// 短絡容量行列 (Maxwell 行列)
EXTERN double *Lmat;			// インダクタンス行列
EXTERN double *Gmat;			// コンダクタンス行列
EXTERN double *Rmat;			// 抵抗行列 (= inv(G))
EXTERN int HaveC, HaveL, HaveR;	// 各行列が計算済みか
EXTERN double TlineLength;		// 伝送線路長 [m] (Tline 指定時)

#ifdef __cplusplus
}
#endif

#endif		// _FEM_H_
