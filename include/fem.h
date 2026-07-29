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
#define FN_field "ofe_field.vtk"
#define FN_sweep "ofe_sweep.csv"

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
#define ANALYSIS_E (1 << 5)		// 辺要素 (Nedelec) の自己検証 (3 次元渦電流の基盤)
#define ANALYSIS_A (1 << 6)		// 3 次元渦電流 (A-φ、辺要素) -> R(f), L(f)
#define ANALYSIS_P (1 << 7)		// 節点要素 (P1/P2) の自己検証 (多項式再現の恒等式)

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

// 分散材料の極 (時間因子 e^{jωt}、ε = ε' - jε'')
//   type 1 (Debye)     : Δε / (1 + jωτ)                   a=Δε,  b=τ [s]
//   type 2 (Lorentz)   : Δε ω0^2 / (ω0^2 - ω^2 + jωδ)     a=Δε,  b=f0 [Hz], c=δ [Hz]
//   type 3 (Drude)     : -ωp^2 / (ω^2 - jωΓ)              a=fp [Hz], b=Γ/2π [Hz]
//   type 4 (Cole-Cole) : Δε / (1 + (jωτ)^(1-α))           a=Δε,  b=τ [s],  c=α (0<=α<1)
//     α = 0 の Cole-Cole は Debye に厳密に一致する (検証で使う)
//   type 5 (Havriliak-Negami) : Δε / (1 + (jωτ)^(1-α))^β
//                                                          a=Δε, b=τ, c=α, d=β (0<β<=1)
//     β = 1 で Cole-Cole、α = 0 で Cole-Davidson、両方で Debye に厳密一致する。
//     この 3 つの極限がそのまま検証の恒等式になる
typedef struct {
	int    type;
	double a, b, c, d;
} pole_t;

// 材料 (id=0 : 真空, id=1 : PEC 予約, id>=2 : ユーザー定義)
typedef struct {
	double epsr;				// 比誘電率
	double sigma;				// 導電率 [S/m]
	double mur;					// 比透磁率 (静磁場解析で使う、既定 1)
	double tand;				// 誘電正接 tanδ (frequency 指定時に G へ寄与、既定 0)
	// 鉄損 (Bertotti の損失分離)。bertotti キーで与える。0 なら評価しない
	//   P = kh f B^alpha  +  (pi^2 sigma d^2 / 6) f^2 B^2  +  ke (f B)^1.5   [W/m^3]
	//        ヒステリシス            古典渦電流 (積層厚 d)        異常 (過剰)
	double be_kh, be_alpha, be_ke, be_d;
	int    npole;				// 分散極の数 (0 : 分散なし)
	double einf;				// ε∞
	pole_t pole[MAXPOLE];
	// 異方性テンソル (対称、成分順 xx, yy, zz, xy, yz, zx)。既定は等方性
	int    eps6_given;			// 1 : anisoeps が明示された (周波数掃引で上書きしない)
	double eps6[6];				// 比誘電率テンソル
	double mu6[6];				// 比透磁率テンソル
	// 導電率の温度依存 : σ(T) = σ0 / (1 + α (T - T0))  (金属の標準的な抵抗率モデル)
	double tempco;				// 抵抗率の温度係数 α [1/K] (0 : 温度依存なし)
	double temp0;				// 基準温度 T0 [degC]
	// εr / μr / B-H の温度依存 (いずれも 1 次係数、データシートの ppm/K 相当)
	//   εr(T) = εr0 (1 + αe (T - Te0))     epstempco
	//   μr(T) = μr0 (1 + αm (T - Tm0))     mutempco
	//   B(H; T) = B0(H) (1 + αb (T - Tb0)) bhtempco  (飽和磁束密度の温度低下)
	//
	// **εr / μr は読み出し時 (material_coef_pub) に掛ける。** σ のように入力解釈で
	// 一度だけ掛ける方式が使えないのは、material_freq() が分散材料の εr を
	// 毎回展開し直して上書きするため (周波数掃引で補正が消える)。
	// B-H 曲線は再計算されないので σ と同じく入力解釈で一度だけ掛ける。
	double epstempco, epstemp0;
	double mutempco,  mutemp0;
	double bhtempco,  bhtemp0;
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

// 格子の種類 (0 : 構造格子, 1 : 非構造格子 (四面体))
EXTERN int MeshMode;

// 非構造格子 (Gmsh ASCII 2.2 を読む)
EXTERN char MeshFile[BUFSIZ];
EXTERN int NNode;				// 節点数
EXTERN double *Xp, *Yp, *Zp;	// 節点座標
EXTERN int NTet;				// 四面体数
EXTERN int32_t *Tet;			// [4*NTet] 頂点の節点番号
EXTERN int *TetTag;				// [NTet] 物理タグ
EXTERN unsigned char *TetMat;	// [NTet] 材料番号
// 要素次数 (1 : 4 節点四面体 / 2 : 10 節点四面体)。2 次のとき Tet2 / Tri2 に
// 辺上の中間節点が入る。局所の辺の並びは Gmsh の tet10 / tri6 と同じ
//   Tet2 : (0,1) (1,2) (2,0) (3,0) (3,2) (3,1)
//   Tri2 : (0,1) (1,2) (2,0)
EXTERN int TetOrder;
EXTERN int32_t *Tet2;			// [6*NTet] 2 次のときのみ (1 次では NULL)
EXTERN int32_t *Tri2;			// [3*NTri] 2 次のときのみ (1 次では NULL)
// 辺 (3 次元渦電流の辺要素で使う)。辺は節点番号の小さい方から大きい方へ向ける
EXTERN int NEdge;
EXTERN int64_t *EdgePtr;		// [NNode+1] 節点 i を始点とする辺の範囲
EXTERN int32_t *EdgeTo;			// [NEdge] 終点 (昇順)
EXTERN int32_t *EdgeFrom;		// [NEdge] 始点 (節点番号の小さい方)
EXTERN int32_t *TetEdge;		// [6*NTet] 四面体の辺番号
EXTERN signed char *TetEdgeSgn;	// [6*NTet] 局所の並びと全体の向きの符号

EXTERN int NTri;				// 三角形数 (電極面の指定に使う)
EXTERN int32_t *Tri;			// [3*NTri]
EXTERN int *TriTag;				// [NTri] 物理タグ

// 格子の次元 (3 : 四面体、2 : 断面 2 次元の三角形)。2 のときは三角形が
// 体積要素になり、M / F (断面 2 次元の定式化) を非構造格子で解ける。
// 面は伝送線路軸 Tline に垂直で、その軸の座標は一定でなければならない
EXTERN int MeshDim;
EXTERN unsigned char *TriMat;	// [NTri] 材料番号 (2 次元格子)
EXTERN signed char *TriCond;	// [NTri] 導体番号 (2 次元格子、-1 = 導体でない)
EXTERN double *TriArea;			// [NTri] 面積 [m^2]

// 物理タグ -> 材料 / 導体の対応
#define MAXTAGMAP (64)
EXTERN int NRegion, RegionTag[MAXTAGMAP], RegionMat[MAXTAGMAP];
EXTERN int NElectrode, ElecTag[MAXTAGMAP], ElecCond[MAXTAGMAP];
// A の接線成分を 0 に固定する面 (B・n = 0 の対称面 / 磁束を通さない境界)
EXTERN int NAWall, AWallTag[MAXTAGMAP];
// 3 次元渦電流のゲージ固定 (tree-cotree)。0 = 固定しない (既定)
EXTERN int GaugeTree;

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
// 直接解法 (RCM + スカイライン Cholesky)。1 : 反復解法の代わりに使う (既定 0)。
// 実対称正定値の系 (C / L / R / M) 専用。渦電流 (F / A) の複素系は COCG のまま
EXTERN int Direct;

EXTERN int Analysis;			// ANALYSIS_* のビット OR
EXTERN char Tline;				// 'X'/'Y'/'Z' : 単位長あたりで出力、0 : 絶対値
EXTERN double LineLength;		// 等価回路を作るときの線路長 [m] (既定 : 解析領域長)
EXTERN double Volt;				// 励振電圧 [V]
EXTERN double Freq;				// 周波数 [Hz] (tanδ による誘電損の計算に使う、0 = 無効)
EXTERN double Curr;				// 静磁場解析の励振電流 [A]
EXTERN double Temperature;		// 動作温度 [degC] (tempco と併せて σ を補正する)

// 周波数掃引 (frequencysweep)。各点で分散材料を展開し直して解き直す
EXTERN int NFreqSweep;
EXTERN double FreqSweep[MAXSWEEP];
EXTERN int NSection;			// SPICE 等価回路の梯子段数

// 場の出力 (fieldout キー)。solve() の途中でしか手元に無い解を溜めて最後に書く
#define MAXFIELD (16)
EXTERN int FieldOut;			// 1 : ofe_field.vtk を書く (既定 0)
EXTERN int NFieldN, NFieldC;
EXTERN char FieldNName[MAXFIELD][64];
EXTERN double *FieldN[MAXFIELD];	// [num_node()] 節点スカラー
EXTERN char FieldCName[MAXFIELD][64];
EXTERN double *FieldC[MAXFIELD];	// [3*num_cell()] 要素ベクトル

// 入力解釈で出た警告 (ofe.log を開く前に走るので溜めておき、後でログにも出す)
#define MAXINWARN (16)
EXTERN int NInputWarn;
EXTERN char InputWarn[MAXINWARN][BUFSIZ];

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
EXTERN double CondTempco[MAXPORT];	// 抵抗率の温度係数 α [1/K] (0 = 温度依存なし)
EXTERN double CondTemp0[MAXPORT];	// 基準温度 T0 [degC]
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
// 鉄損 (bertotti)。ポート k を励振したときの損失を対角 [k][k] に入れる。
// 鉄損は B の非線形関数なので重ね合わせが効かず、非対角は定義できない
EXTERN double *Pfemat;			// 鉄損 [W/m] (単位長あたり) または [W]
EXTERN int HavePfe;
EXTERN double TlineLength;		// 伝送線路長 [m] (Tline 指定時)

#ifdef __cplusplus
}
#endif

#endif		// _FEM_H_
