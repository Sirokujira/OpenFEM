#!/bin/sh
# rlc_check.sh — 解析解との比較による検証 (CI 用)
#
# data/sample/ の検証ケースを実行し、ofe_post が出力する rlc.csv の値を
# 解析解と比較する。
#
#   parallel_plate : C = eps0 * epsr * A / d        (厳密、許容 1%)
#   resistor_bar   : R = d / (sigma * A)            (厳密、許容 1%)
#   coax           : C' = 2 pi eps / ln(b/a)        (階段近似、許容 5%)
#                    L' = mu0 / (2 pi) ln(b/a)      (階段近似、許容 5%)
#   coupled_microstrip : 唯一の多ポートケース。形状の対称性から C11=C22, L11=L22 が
#                    全桁一致すること + Maxwell 行列の符号 + 多ポート SPICE (K, CM)
#   microstrip     : Z0 が 40..60 ohm の範囲にあること (設計値 ~50 ohm)
#   plate_line_dc  : L'dc = mu0(d + 2t/3)/W         (1 次元厳密、許容 1%)
#                    Rs' = 2/(sigma W t)            (厳密、許容 1%)
#   coax_loss      : L'dc (内部インダクタンス込み)  (階段近似、許容 5%)
#                    Rs' (導体 DC 抵抗)             (階段近似、許容 3%)
#                    G' = omega C' tand             (均質誘電体では厳密、許容 0.1%)
#   plate_line_ac  : R(f), L(f) を 1 次元厳密解と比較 (1kHz / 10MHz、許容 2%)
#   plate_line_bh  : 非線形磁性体 L(I) を 1 次元厳密解と比較 (4 電流、許容 1%)
#   dispersive_plate : 多極分散 (Debye+Lorentz) の C, G  (厳密、許容 0.1%)
#   drude_plate    : Drude 媒質の C, G と、低周波で σ = ε0ωp^2/Γ の導体に収束すること
#   colecole_plate : Cole-Cole の C, G。α=0 が Debye に厳密一致することも見る
#   temp_resistor  : σ(T) = σ0/(1+α(T-T0)) で R が厳密に比例すること (4 温度)
#   sweep_plate    : 周波数掃引。各点が閉形式と一致し、かつ個別実行と完全一致すること
#   fieldout       : 書き出した場から集中定数を作り直して元の値と比べる
#                    (∫½ε|E|²dV = ½CV²、∫|J|²/(2σ)dV = ½Re(VI*))
#   input lint     : 選んだ解析が読まないキーを警告すること (5 つの罠) と、
#                    正しい 20 ケースで警告が 1 件も出ないこと
#   aniso_plate    : 異方性誘電体の C = eps0 εz A/d      (厳密、許容 0.1%)
#   aniso_rot      : 非対角テンソル (z 軸 30 度回転) でも同値 (厳密、許容 0.1%)
#   plate_line_bh_aniso : 軸毎 B-H。B は x のみなので X 曲線だけが効くこと
#   plate_line_ja  : Jiles-Atherton の履歴ループを、H 掃引で独立に積分した
#                    ODE 解と比較する (11 点、許容 0.5%)
#   box_tet        : 非構造格子 (四面体) の平行平板 C     (厳密、許容 0.1%)
#   coax_tet       : 円形境界に適合した四面体格子の同軸 C', L' (許容 1%)
#   edge_test      : Whitney 辺要素 (1 次 Nedelec) の自己検証 (機械精度)
#                    + ゲージ固定 (tree-cotree) と Hiptmair 前処理の検証
#   edge_test_aniso : 同じ自己検証を非対角成分まで詰まった異方性 ν で行う
#                    (等方性だけだと異方性項が一度も実行されない)
#   nodal_test_p1/p2 : 節点要素 (P1 / P2) の剛性行列を多項式再現の恒等式と
#                    比較する (φ^T K φ = ∫(∇φ)^T C (∇φ) dV、機械精度)
#   nodal_test_coax : 曲がった 2 次格子で積分した体積を円環の解析値と比較する
#                    (等パラメトリック写像のヤコビアンの検証、許容 0.1%)
#   box_p2         : 2 次要素を通した求解の一式 (平行平板 C、許容 0.1%)
#   coax_p2        : 粗い同軸 (nr=4, nt=12) を曲がった 2 次要素で (許容 0.3%)
#                    + 同じ格子の 1 次要素より 10 倍以上良いこと
#   mesh order     : 次数の混在・2 次格子への analysis=A・1 次三角形を弾くこと
#   direct         : 直接解法 (RCM + スカイライン Cholesky) が反復解法と
#                    完全一致すること (実対称の 9 ケース)
#   hdf5           : 系列の HDF5 出力 (任意依存)。既存の出口と同じ数字が入って
#                    いること (ofe_sweep.csv / ofe.log の履歴 / ofe_field.vtk)。
#                    HDF5 無しのビルドでは「弾かれること」だけ見て skip する
#   direct F/A     : 複素対称系の直接解法 (LDL^T)。反復解法と一致すること
#                    (F の構造格子・三角形格子、3 次元 A-φ) + 分解の残差が
#                    小さいこと + ゲージ固定なしの A を弾くこと
#   gmsh 4.1       : 同じ形状を Gmsh 2.2 と 4.1 で書いた結果が完全一致すること
#                    (+ $Entities 欠落と、ヘッダだけバイナリを名乗る
#                    ファイルを弾くこと)
#   plate2d_p2     : 断面 2 次元の 2 次要素 (6 節点三角形)。導体内の Az が厳密に
#                    2 次なので内部インダクタンスまで厳密に出る (1 次は -0.045%)
#   plate2d        : 断面 2 次元の三角形格子で M / F を解き、構造格子版と同じ
#                    1 次元厳密解と比較する (Ldc / Rs / R(f) / L(f)、許容 0.2%)
#   anisotropic mu : 面内の異方性 ν が **B に掛かる** こと (grad(Az) ではなく)。
#                    等方性では一致するので等方性ケースでは検出できない
#   post sweep     : rlc.csv が「掃引の最後の 1 点」であることを明記すること
#                    (掃引していないときは出ないことも見る)
#   hn_plate       : Havriliak-Negami 分散。β=1 で Cole-Cole、α=0 で Cole-Davidson、
#                    両方で Debye に厳密一致すること (極限が恒等式になる)
#   temp_material  : εr(T) で C が厳密に比例すること (4 温度)
#   temp_mur       : μr(T) で L が動くこと (内部インダクタンスは動かない)
#   bhtempco       : B 軸を k 倍すると (L - L_int) がちょうど k 倍になること
#   bertotti_core  : 鉄損 (Bertotti の損失分離)。1 次元厳密な B に対する閉形式と
#                    機械精度で一致すること + **3 項を指数で分離**して検査する
#   off-diagonal mu : 格子と材料テンソルを同じ角だけ面内で回すと答えが一致する
#                    こと (合同な離散問題なので厳密)。非対角成分を検査できるのは
#                    斜めに回した非対称断面だけ
#   bar_eddy       : 3 次元渦電流 (A-φ)。Z = γL/(2σW tanh(γt/2)) と比較
#                    (1e2 / 1e4 / 1e5 Hz、許容 0.5% / 0.5% / 2%)
#                    + gauge = 1 で Z が変わらないこと (ゲージ不変性)
#                    + 磁性導体 (mur = 50) の Z と、表皮深さ警告が出ること
#   bar_air        : 非導電層 (空気) を含む A-φ。空気層が Robin 条件に潰れる
#                    1 次元閉形式と比較 (1e2/1e4/1e5 Hz + mur=50 + gauge=1)
#
# 使い方 : rlc_check.sh <ofe 実行ファイル> <ofe_post 実行ファイル> [作業ディレクトリ]

set -e

OFE="$1"
OFE_POST="$2"
WORK="${3:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$OFE" ] || [ -z "$OFE_POST" ]; then
	echo "Usage: rlc_check.sh <ofe> <ofe_post> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"
status=0

# rlc.csv から行列 <name> の (1,1) 成分を取り出す
value_of() {
	awk -F, -v name="$1" '
		$1 == name && NF == 2 { found = 1; next }
		found && $1 == "1"    { print $2; exit }
	' "$WORK/rlc.csv"
}

# rlc.csv から行列 <name> の (i,j) 成分を取り出す (多ポートの非対角を見るため)
value_ij() {
	awk -F, -v name="$1" -v i="$2" -v j="$3" '
		$1 == name && NF == 2 { found = 1; next }
		found && ($1 == i)    { print $(j + 1); exit }
	' "$WORK/rlc.csv"
}

# スカラー項目 (例 "Z0 [ohm]") を取り出す
scalar_of() {
	awk -F, -v name="$1" '$1 == name { print $2; exit }' "$WORK/rlc.csv"
}

run_case() {
	cp "$SRC/$1.ofe" "$WORK/"
	# 非構造格子のケースはメッシュファイルも要る
	for m in "$SRC"/*.msh; do
		[ -f "$m" ] && cp "$m" "$WORK/"
	done
	(cd "$WORK" && "$OFE" -n 2 "$1.ofe" > /dev/null && "$OFE_POST" > /dev/null)
}

compare() {
	# compare <label> <value> <expected> <tolerance>
	res=$(awk -v v="$2" -v e="$3" -v tol="$4" \
		'BEGIN{ if (v == "") { print "NG (no value)"; exit }
		        d = (v - e) / e; a = (d < 0) ? -d : d;
		        printf "%s %+.2f%%", (a <= tol) ? "OK" : "NG", d * 100 }')
	echo "  $1 : value=$2 expected=$3 -> $res"
	case "$res" in NG*) status=1 ;; esac
}

# ofe_field.vtk から集計値を取り出す (vtkcheck.awk の出力を 1 項目だけ拾う)
vtk() { awk -v arr="$2" -v axis="${3:-0}" -v xcut="${4:-}" -v comp="${5:-}" \
	-f "$SRC/vtkcheck.awk" "$WORK/ofe_field.vtk" | awk -v k="$1" '$1 == k { print $2 }'; }

in_range() {
	# in_range <label> <value> <min> <max>
	res=$(awk -v v="$2" -v lo="$3" -v hi="$4" \
		'BEGIN{ if (v == "") { print "NG (no value)"; exit }
		        printf "%s", ((v >= lo) && (v <= hi)) ? "OK" : "NG" }')
	echo "  $1 : value=$2 range=[$3,$4] -> $res"
	case "$res" in NG*) status=1 ;; esac
}

echo "[parallel_plate] C = eps0 * 4.0 * 1e-6 / 0.2e-3"
run_case parallel_plate
compare "C [F]" "$(value_of C)" 1.7708376e-13 0.01

echo "[resistor_bar] R = 0.2e-3 / (1e3 * 1e-6)"
run_case resistor_bar
compare "R [ohm]" "$(value_of R)" 2.0e-1 0.01

echo "[coax] a = 0.5mm, b = 1.5mm, epsr = 2.1"
run_case coax
compare "C [F/m]" "$(value_of C)" 1.063395e-10 0.05
compare "L [H/m]" "$(value_of L)" 2.1972246e-07 0.05
compare "Z0 [ohm]" "$(scalar_of 'Z0 [ohm]')" 4.545953e+01 0.05

# 結合線路 (2 ポート)。**唯一の多ポートケース**なので、行列の非対角成分・
# 対称化の診断・多ポートの SPICE 出力 (相互インダクタンス K と結合容量 CM) を
# 通すのはここだけ。README に載っているのに検査が無く、壊れても気づけなかった。
#
# 使う恒等式は「形状の対称性」。2 本の線路は x = 0 について対称で xmesh も
# 対称なので、離散問題まで対称になり C11 = C22 と L11 = L22 が**全桁一致**する。
# ポート番号の取り違えはここで落ちる。
echo "[coupled_microstrip] 2-port coupled line (the only multi-port case)"
run_case coupled_microstrip
compare "C11 == C22 (structure is symmetric)" "$(value_ij C 1 1)" "$(value_ij C 2 2)" 1e-9
compare "L11 == L22 (structure is symmetric)" "$(value_ij L 1 1)" "$(value_ij L 2 2)" 1e-9
compare "C12 == C21 (matrix symmetry)" "$(value_ij C 1 2)" "$(value_ij C 2 1)" 1e-9
compare "L12 == L21 (matrix symmetry)" "$(value_ij L 1 2)" "$(value_ij L 2 1)" 1e-9
# Maxwell 容量行列の符号 : 対角 > 0、非対角 < 0、行和 > 0 (対地容量が正)。
# 相互インダクタンスは同方向電流なので正
res=$(awk -v c11="$(value_ij C 1 1)" -v c12="$(value_ij C 1 2)" -v l12="$(value_ij L 1 2)" \
	'BEGIN{ printf "%s", ((c11 > 0) && (c12 < 0) && ((c11 + c12) > 0) && (l12 > 0)) ? "OK" : "NG" }')
echo "  Maxwell signs (C diag>0, C offdiag<0, row sum>0, L offdiag>0) : $res"
case "$res" in NG*) status=1 ;; esac
# TEM の L は真空容量の逆行列 L = mu0 eps0 inv(C0) なので、**inv(L) も Maxwell 行列**
# (非対角が負) でなければならない。2x2 では inv([[a,b],[b,a]]) の非対角が -b/(a^2-b^2)
res=$(awk -v a="$(value_ij L 1 1)" -v b="$(value_ij L 1 2)" \
	'BEGIN{ d = (a * a) - (b * b)
	        printf "%s", ((d > 0) && ((-b / d) < 0)) ? "OK" : "NG" }')
echo "  inv(L) is a Maxwell matrix (L = mu0 eps0 inv(C0)) : $res"
case "$res" in NG*) status=1 ;; esac
# 多ポートの SPICE 出力 : 相互インダクタンス K と結合容量 CM が段数分だけ出ること
res=$(awk -v ns=4 '/^K12_/ { k++ } /^CM12_/ { c++ }
	END { printf "%s (K=%d, CM=%d, sections=%d)", ((k == ns) && (c == ns)) ? "OK" : "NG",
	             k, c, ns }' "$WORK/ofe_circuit.sp")
echo "  SPICE ladder has K and CM per section : $res"
case "$res" in NG*) status=1 ;; esac
# 結合係数 k = L12/L11 が妥当な範囲にあること (w = s = h なので 0.1..0.3)
in_range "coupling k = L12/L11" \
	"$(awk -v a="$(value_ij L 1 1)" -v b="$(value_ij L 1 2)" 'BEGIN{ printf "%.6f", b / a }')" \
	0.1 0.3
# 対称化の診断が出ていて、かつ十分小さいこと (多ポートでしか出ない行)
asym=$(awk '/electrostatic matrix asymmetry/ { print $(NF-1) }' "$WORK/ofe.log")
res=$(awk -v v="$asym" 'BEGIN{ if (v == "") { print "NG (no diagnostic)"; exit }
	                            printf "%s", (v < 1e-9) ? "OK" : "NG" }')
echo "  matrix asymmetry diagnostic = $asym -> $res"
case "$res" in NG*) status=1 ;; esac

echo "[microstrip] w/h = 0.75/0.4, epsr = 4.4"
run_case microstrip
in_range "Z0 [ohm]" "$(scalar_of 'Z0 [ohm]')" 40 60
in_range "eps_eff" "$(scalar_of 'eps_eff')" 2.5 4.0

echo "[plate_line_dc] W=1mm, d=0.2mm, t=0.05mm (magnetostatic, 1-D exact)"
run_case plate_line_dc
compare "Ldc [H/m]" "$(value_of Ldc)" 2.932153e-07 0.01
compare "Rs [ohm/m]" "$(value_of Rs)" 6.896552e-01 0.01

echo "[coax_loss] a=0.5mm, b=1.5mm, c=1.9mm, tand=0.001 @ 1GHz"
run_case coax_loss
compare "Ldc [H/m]" "$(value_of Ldc)" 2.873953e-07 0.05
compare "Rs [ohm/m]" "$(value_of Rs)" 2.598776e-02 0.03
# 均質誘電体では G' = omega * C' * tand が厳密に成り立つ (離散化誤差が相殺する)
gexp=$(awk -v c="$(value_of C)" -v f=1e9 -v td=0.001 \
	'BEGIN{ printf "%.8e", 2 * 3.14159265358979324 * f * c * td }')
compare "G [S/m]" "$(value_of G)" "$gexp" 0.001

# 渦電流 (表皮効果) : 平行平板線路の 1 次元厳密解と比較する
#   Z = 2*gamma/(sigma W) coth(gamma t) + j omega mu0 d/W,  gamma = sqrt(j omega mu0 sigma)
echo "[plate_line_ac] skin effect (1-D exact internal impedance)"
for pair in "1e3 6.896552e-01 2.932153e-07" "1e7 1.624296e+00 2.780550e-07"; do
	set -- $pair
	freq=$1
	rexp=$2
	lexp=$3
	sed "s/^frequency = .*/frequency = $freq/" "$SRC/plate_line_ac.ofe" > "$WORK/plate_line_ac_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 plate_line_ac_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R(f=$freq) [ohm/m]" "$(value_of Rf)" "$rexp" 0.02
	compare "L(f=$freq) [H/m]" "$(value_of Lf)" "$lexp" 0.02
done

echo "[dispersive_plate] Debye + Lorentz multi-pole at 1 GHz"
run_case dispersive_plate
compare "C [F]" "$(value_of C)" 2.123808563e-13 0.001
compare "G [S]" "$(value_of G)" 4.773870170e-04 0.001

# Drude 媒質 : 低周波では σ = ε0 ωp^2/Γ の導体に厳密に収束する
echo "[drude_plate] Drude medium (fp = 2 GHz, gamma = 5 GHz) at 1 GHz"
run_case drude_plate
compare "C [F]" "$(value_of C)" 1.702728426e-13 0.001
compare "G [S]" "$(value_of G)" 2.139711646e-04 0.001
# ω << Γ の極限で G -> σ A/d、σ = ε0 ωp^2/Γ = 4.4506002242e-02 S/m
sed "s/^frequency = .*/frequency = 1e6/" "$SRC/drude_plate.ofe" > "$WORK/drude_lf.ofe"
(cd "$WORK" && "$OFE" -n 2 drude_lf.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "G(f=1MHz) [S] -> sigma A/d" "$(value_of G)" 2.2253001121e-04 0.001
# ω < ωp では εr' < 0 になり準静的定式化が成り立たないので落ちること
sed -e "s/^drude = .*/drude = 2 4.0 2e9 1e8/" -e "s/^frequency = .*/frequency = 5e8/" \
    "$SRC/drude_plate.ofe" > "$WORK/drude_neg.ofe"
if (cd "$WORK" && "$OFE" -n 2 drude_neg.ofe > /dev/null 2>&1); then
	echo "  epsr' < 0 rejected : no -> NG" >&2
	status=1
else
	echo "  epsr' < 0 rejected : yes -> OK"
fi

# Cole-Cole 分散。α = 0 が Debye に厳密一致することと、ωτ = 1 で実部が α に
# 依らないことの 2 つを恒等式として使う
echo "[colecole_plate] Cole-Cole (alpha = 0.3) at 1 GHz"
run_case colecole_plate
compare "C [F]" "$(value_of C)" 1.549482868e-13 0.001
compare "G [S]" "$(value_of G)" 2.556873117e-04 0.001
sed "s/^colecole = .*/colecole = 2 2.0 3.0 1.591549431e-10 0.0/" \
    "$SRC/colecole_plate.ofe" > "$WORK/cc_a0.ofe"
(cd "$WORK" && "$OFE" -n 2 cc_a0.ofe > /dev/null && "$OFE_POST" > /dev/null)
cc0c=$(value_of C); cc0g=$(value_of G)
sed "s/^colecole = .*/debye = 2 2.0 3.0 1.591549431e-10/" \
    "$SRC/colecole_plate.ofe" > "$WORK/cc_db.ofe"
(cd "$WORK" && "$OFE" -n 2 cc_db.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "C (alpha=0) == Debye" "$cc0c" "$(value_of C)" 1e-12
compare "G (alpha=0) == Debye" "$cc0g" "$(value_of G)" 1e-12
# ωτ = 1 では実部が α に依らない (虚部だけ動く)
sed "s/^colecole = .*/colecole = 2 2.0 3.0 1.591549431e-10 0.6/" \
    "$SRC/colecole_plate.ofe" > "$WORK/cc_a6.ofe"
(cd "$WORK" && "$OFE" -n 2 cc_a6.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "C (alpha=0.6) == C (alpha=0)" "$(value_of C)" "$cc0c" 1e-12
compare "G (alpha=0.6) [S]" "$(value_of G)" 1.355707190e-04 0.001

# 導電率の温度依存 σ(T) = σ0/(1+α(T-T0))。R はこの因子にちょうど比例する
echo "[temp_resistor] sigma(T) = sigma0/(1+alpha(T-T0)), copper alpha"
for pair in "20 2.0000000e-01" "85 2.5109000e-01" "-40 1.5284000e-01" "125 2.8253000e-01"; do
	set -- $pair
	sed "s/^temperature = .*/temperature = $1/" "$SRC/temp_resistor.ofe" > "$WORK/temp_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 temp_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R(T=$1 degC) [ohm]" "$(value_of R)" "$2" 0.001
done
# CondSigma[] は Material[].sigma とは別系統なので、Rs 側も温度が効くこと
sed -e "s/^analysis = .*/conductortempco = 0 3.93e-3 20\nconductortempco = 1 3.93e-3 20\ntemperature = 85\n&/" \
    "$SRC/plate_line_dc.ofe" > "$WORK/temp_rs.ofe"
(cd "$WORK" && "$OFE" -n 2 temp_rs.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "Rs(T=85 degC) [ohm/m]" "$(value_of Rs)" 8.658275862e-01 0.001

echo "[aniso_plate] anisotropic eps (10, 5, 2), field along z"
run_case aniso_plate
compare "C [F]" "$(value_of C)" 8.854188e-14 0.001

echo "[aniso_rot] same tensor rotated 30 deg about z (off-diagonal exy)"
run_case aniso_rot
compare "C [F]" "$(value_of C)" 8.854188e-14 0.001

# 非線形磁性体 (B-H) : 平行平板線路の 1 次元厳密解と比較する
#   H = I/W が Ampere の法則から厳密に決まるので L = B(I/W)*d/I + 2*mu0*t/(3W)
echo "[plate_line_bh] nonlinear B-H (1-D exact, secant inductance)"
for pair in "0.5 2.889308e-04" "1.0 2.000419e-04" "2.0 1.022641e-04" "5.0 4.359744e-05"; do
	set -- $pair
	cur=$1
	lexp=$2
	sed "s/^current = .*/current = $cur/" "$SRC/plate_line_bh.ofe" > "$WORK/plate_line_bh_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 plate_line_bh_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "L(I=${cur}A) [H/m]" "$(value_of Ldc)" "$lexp" 0.01
	# 収束しなかった場合はログに警告が出る
	if grep -q "did not converge" "$WORK/ofe.log"; then
		echo "  *** B-H iteration did not converge (I=$cur)" >&2
		status=1
	fi
done

# 直交異方性の非線形磁性体 : B は x 方向のみなので X 軸の曲線だけが効く
# (Y/Z には 10 倍硬い曲線を与えてあるので、軸を取り違えると大きくずれる)
echo "[plate_line_bh_aniso] per-axis B-H (only the X curve must matter)"
for pair in "0.5 2.889308e-04" "1.0 2.000419e-04" "2.0 1.022641e-04" "5.0 4.359744e-05"; do
	set -- $pair
	cur=$1
	lexp=$2
	sed "s/^current = .*/current = $cur/" "$SRC/plate_line_bh_aniso.ofe" > "$WORK/bh_aniso_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 bh_aniso_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "L(I=${cur}A) [H/m]" "$(value_of Ldc)" "$lexp" 0.01
done

# 非構造格子 (四面体)。1 次四面体は 1 次元電界を厳密に表せるので構造格子と同値
echo "[box_tet] parallel plate on a tetrahedral mesh"
run_case box_tet
compare "C [F]" "$(value_of C)" 1.7708376e-13 0.001

# 円形境界に適合するので階段近似の誤差 (構造格子で -2.11%) が消える
echo "[coax_tet] coax on a conforming tetrahedral mesh"
run_case coax_tet
compare "C [F/m]" "$(value_of C)" 1.063417e-10 0.01
compare "L [H/m]" "$(value_of L)" 2.1972246e-07 0.01

# 辺要素 (Nedelec) の自己検証。3 次元渦電流 (A-φ) の基盤
# 勾配の零空間・一様場の質量・回転場の回転回転と質量・対称性を閉形式と比較する
edge_selftest() {
	run_case "$1"
	grad=$(awk  '/gradient null space/ { print $NF }' "$WORK/ofe.log")
	mass=$(awk  '/uniform field mass/  { print $NF }' "$WORK/ofe.log")
	curl0=$(awk '/uniform field curl/  { print $NF }' "$WORK/ofe.log")
	rot=$(awk   '/rotational field/    { print $NF }' "$WORK/ofe.log")
	rotm=$(awk  '/rotational mass/     { print $NF }' "$WORK/ofe.log")
	# (d) の行は S と T の 2 値を持つ ($6 が S、$NF が T)
	syms=$(awk '/symmetry  / { print $6  }' "$WORK/ofe.log")
	symt=$(awk '/symmetry  / { print $NF }' "$WORK/ofe.log")
	for pair in "gradient-null:$grad" "mass:$mass" "uniform-curl:$curl0" \
	            "rotational-curl:$rot" "rotational-mass:$rotm" \
	            "symmetry-S:$syms" "symmetry-T:$symt"; do
		nm=${pair%%:*}
		vv=${pair#*:}
		res=$(awk -v v="$vv" 'BEGIN{ if (v == "") { print "NG (no value)"; exit }
		                             printf "%s", (v < 1e-10) ? "OK" : "NG" }')
		echo "  $nm : $vv -> $res"
		case "$res" in NG*) status=1 ;; esac
	done
	if ! grep -q "edge element self test passed" "$WORK/ofe.log"; then
		echo "  *** edge element self test failed" >&2
		status=1
	fi
}

echo "[edge_test] Nedelec edge element self test"
edge_selftest edge_test
# ゲージ固定と前処理の要約を表示する (合否は自己検証の総合判定で見る)
awk '/\(e\) spanning tree/ { print "  " $0 }
     /gradient energy/       { print "  " $0 }
     /Jacobi   :/            { print "  " $0 }
     /Hiptmair :/            { print "  " $0 }
     /-> .*fewer iterations/ { print "  " $0 }' "$WORK/ofe.log"

# 等方性材料では ν の非対角成分が恒等的に 0 で異方性項が実行されないため、
# 6 成分すべてが非零な μ でもう一度回して成分順序と対称性を検証する
echo "[edge_test_aniso] the same self test with a fully anisotropic nu"
edge_selftest edge_test_aniso

# 節点要素 (P1 / P2) の自己検証。多項式の再現性の恒等式なので機械精度で合う。
# 右辺は要素毎の閉形式で、組み立てとは独立に計算している
nodal_selftest() {	# nodal_selftest <ケース> <2 次の場を検査するか yes|no>
	run_case "$1"
	null=$(awk '/constant null space/ { print $NF }' "$WORK/ofe.log")
	lin=$(awk  '/linear    field/     { print $NF }' "$WORK/ofe.log")
	quad=$(awk '/quadratic field .*err/ { print $NF }' "$WORK/ofe.log")
	sym=$(awk  '/\(d\) symmetry/      { print $NF }' "$WORK/ofe.log")
	# (e) は 1 次と 2 次の φ で 1 行ずつ出るので、悪い方 (最大) を採る
	grd=$(awk  '/centroid gradient/   { if ($NF > m) m = $NF } END { print m + 0 }' "$WORK/ofe.log")
	list="constant-null:$null linear:$lin symmetry:$sym centroid-grad:$grd"
	[ "$2" = yes ] && list="$list quadratic:$quad"
	for pair in $list; do
		nm=${pair%%:*}
		vv=${pair#*:}
		res=$(awk -v v="$vv" 'BEGIN{ if (v == "") { print "NG (no value)"; exit }
		                             printf "%s", (v < 1e-10) ? "OK" : "NG" }')
		echo "  $nm : $vv -> $res"
		case "$res" in NG*) status=1 ;; esac
	done
	if ! grep -q "result : passed" "$WORK/ofe.log"; then
		echo "  *** nodal element self test failed" >&2
		status=1
	fi
}

echo "[nodal_test_p1] P1 stiffness matrix vs the closed-form polynomial identity"
nodal_selftest nodal_test_p1 no
# 2 次の場は 1 次要素では補間が厳密でないので実行されないこと (検査の空回り防止)
if grep -q "quadratic field .*err" "$WORK/ofe.log"; then
	echo "  *** P1 must not run the quadratic identity" >&2
	status=1
fi

echo "[nodal_test_p2] P2 (10-node tetrahedra), including the quadratic field"
nodal_selftest nodal_test_p2 yes

# 曲がった 2 次格子。中間節点を円筒面に載せてあるので、積分した体積が
# 等パラメトリック写像のヤコビアンを直接検証する。頂点だけで積分すると
# 内接多角形の体積 (-4.51e-2) になるので 290 倍ずれて必ず落ちる
echo "[nodal_test_coax] isoparametric Jacobian on a curved P2 mesh"
nodal_selftest nodal_test_coax no
compare "curved volume [m^3]" \
	"$(awk '/integrated volume/ { print $4 }' "$WORK/ofe.log")" 6.283185307e-10 0.001
# 曲がっている格子では 2 次の場の検査を飛ばすこと (飛ばさないと誤検出になる)
if ! grep -q "quadratic field : skipped" "$WORK/ofe.log"; then
	echo "  *** the curved mesh must skip the quadratic identity" >&2
	status=1
fi

# 2 次要素を通した求解の一式。電極面は 6 節点三角形なので、辺上の中間節点まで
# 固定できていないと反作用から出る Q が合わない
echo "[box_p2] parallel plate on a second-order tetrahedral mesh"
run_case box_p2
compare "C [F]" "$(value_of C)" 1.7708376e-13 0.001

# 粗い同軸で 2 次要素の効き目を測る。1 次との対で見ないと「黙って 1 次に
# 退行しても閉形式との比較が通る」状態を検出できない
echo "[coax_p2] coarse coax (nr=4, nt=12), curved second-order elements"
run_case coax_p2
compare "C [F/m]" "$(value_of C)" 1.063417e-10 0.003
compare "L [H/m]" "$(value_of L)" 2.1972246e-07 0.003
cp2=$(value_of C)
echo "[coax_p1c] the same mesh with first-order elements (control)"
run_case coax_p1c
cp1=$(value_of C)
# 1 次の誤差が 2 次の 10 倍以上あること (2 次が実際に効いている証拠)
res=$(awk -v a="$cp1" -v b="$cp2" \
	'BEGIN{ ce = 1.063417e-10
	        e1 = (a - ce) / ce; if (e1 < 0) e1 = -e1
	        e2 = (b - ce) / ce; if (e2 < 0) e2 = -e2
	        printf "P1 %+.3f%% vs P2 %+.3f%% (ratio %.1f) -> %s",
	               100 * e1, 100 * e2, e1 / e2, ((e1 > 10 * e2) ? "OK" : "NG") }')
echo "  order 2 beats order 1 : $res"
case "$res" in *NG) status=1 ;; esac

# 格子の次数に関する入力の検査 (黙って誤答を出すより弾く)
mesh_reject() {	# mesh_reject <ラベル> <ofe ファイル>
	if (cd "$WORK" && "$OFE" -n 2 "$2" > /dev/null 2>&1); then
		echo "  $1 : accepted -> NG" >&2
		status=1
	else
		echo "  $1 : rejected -> OK"
	fi
}
echo "[mesh order] inputs the solver must reject"
# 辺要素 (A) は 1 次四面体の Whitney 形状関数に基づくので 2 次格子は使えない
# 導電率と周波数を入れないと「A には導体が要る」等の別の検査で弾かれてしまい、
# 次数の検査が空回りする (通っているように見えて何も見ていない状態になる)
sed -e 's/^analysis = .*/analysis = A/' \
    -e 's/^material = .*/material = 4.0 5.8e7/' \
    -e 's/^solver = .*/frequency = 1e4\nawall = 10/' \
    "$SRC/box_p2.ofe" > "$WORK/p2_edge.ofe"
mesh_reject "analysis A on an order-2 mesh" p2_edge.ofe
if ! (cd "$WORK" && "$OFE" -n 2 p2_edge.ofe 2>&1 | grep -q "need a first-order mesh"); then
	echo "  *** it was rejected for a different reason than the mesh order" >&2
	status=1
fi
# 次数の混在 (四面体を 1 個だけ 1 次に落とす)
awk '/^[0-9]+ 11 2 / && !done { print $1, 4, 2, $4, $5, $6, $7, $8, $9; done = 1; next } { print }' \
    "$SRC/box_p2.msh" > "$WORK/mixed.msh"
sed 's/^mesh = .*/mesh = mixed.msh/' "$SRC/box_p2.ofe" > "$WORK/p2_mixed.ofe"
mesh_reject "mixed element orders" p2_mixed.ofe
# 四面体だけ 2 次で三角形が 1 次 (電極面の中間節点が固定されない)
awk '/^[0-9]+ 9 2 / { print $1, 2, 2, $4, $5, $6, $7, $8; next } { print }' \
    "$SRC/box_p2.msh" > "$WORK/tri1.msh"
sed 's/^mesh = .*/mesh = tri1.msh/' "$SRC/box_p2.ofe" > "$WORK/p2_tri1.ofe"
mesh_reject "order-2 tetrahedra with order-1 triangles" p2_tri1.ofe

# 断面 2 次元の非構造格子 (三角形が体積要素) で M / F を解く。形状は構造格子版
# plate_line_dc / plate_line_ac と同じなので同じ 1 次元厳密解が使える。
# **伝送線路軸が違う** (構造格子は z、こちらは x) ので、面内 2 軸の取り方を
# 取り違えるとここで落ちる
echo "[plate2d] parallel plate line on a 2-D (cross-section) triangular mesh"
run_case plate2d_dc
compare "Ldc [H/m]" "$(value_of Ldc)" 2.932153e-07 0.002
compare "Rs [ohm/m]" "$(value_of Rs)" 6.896552e-01 0.001

echo "[plate2d_ac] the same cross-section with the eddy-current (F) analysis"
cp "$SRC/plate2d_ac.ofe" "$WORK/"
for m in "$SRC"/*.msh; do
	[ -f "$m" ] && cp "$m" "$WORK/"
done
# freq  R[ohm/m]       L[H/m]         許容 (格子が表皮深さを刻めているので厳しい)
for pair in "1e3 6.896552e-01 2.932153e-07" \
            "1e7 1.624296e+00 2.780550e-07"; do
	set -- $pair
	sed "s/^frequency = .*/frequency = $1/" "$SRC/plate2d_ac.ofe" > "$WORK/plate2d_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 plate2d_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R(f=$1) [ohm/m]" "$(value_of Rf)" "$2" 0.002
	compare "L(f=$1) [H/m]" "$(value_of Lf)" "$3" 0.002
done

# 直接解法 (RCM + スカイライン Cholesky)。**反復解法と同じ答えになること**を
# 恒等式に使う。既存のケースにそのまま適用でき、閉形式も許容誤差も要らない。
# 構造格子・非構造格子 (3 次元 / 断面 2 次元)・異方性・非線形をひと通り通す
echo "[direct] the direct solver must agree with the iterative one"
# direct_same <ケース> [nolog]
#   既定では「プロファイルの報告がログに出ること」も見て、黙って反復解法に
#   落ちていないことを確かめる。非線形 (B-H) の内側解法だけはログハンドルを
#   渡していないので報告が出ない。そのケースは nolog を付けて値の一致だけ見る
#   (経路が通っていることは他の 8 ケースで担保している)
direct_same() {
	run_case "$1"; grep -v '^title' "$WORK/rlc.csv" > "$WORK/it.csv"
	sed 's/^analysis = /direct = 1\nanalysis = /' "$SRC/$1.ofe" > "$WORK/dir.ofe"
	if ! (cd "$WORK" && "$OFE" -n 2 dir.ofe > /dev/null && "$OFE_POST" > /dev/null); then
		echo "  $1 : the direct run failed -> NG" >&2
		status=1
		return
	fi
	if [ "$2" != nolog ] && ! grep -q "direct : profile" "$WORK/ofe.log"; then
		echo "  $1 : the direct solver was not used -> NG" >&2
		status=1
		return
	fi
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/di.csv"
	if cmp -s "$WORK/it.csv" "$WORK/di.csv"; then
		echo "  $1 : identical -> OK"
	else
		echo "  $1 : differ -> NG" >&2
		diff "$WORK/it.csv" "$WORK/di.csv" | head -4 >&2
		status=1
	fi
}
direct_same parallel_plate
direct_same resistor_bar
direct_same aniso_plate
direct_same plate_line_dc
direct_same box_tet
direct_same coax_tet
direct_same plate2d_dc
direct_same plate2d_p2
# 非線形 (B-H) の内側でも使えること (内側解法はログを取らないので nolog)
direct_same plate_line_bh nolog
# **RCM が実際に効いていることも見る。** 並べ替えを外して密に詰めても答えは
# 正しいままなので、値の一致だけでは検出できない (実測: first[] を全部 0 に
# する変異が全緑で通った)。平均帯幅が節点数の 1/4 未満であることを assert する
# (密なら (n+1)/2 = 816 になる。RCM 有りの実測は 57.1)
sed 's/^analysis = /direct = 1\nanalysis = /' "$SRC/coax_tet.ofe" > "$WORK/dbw.ofe"
(cd "$WORK" && "$OFE" -n 2 dbw.ofe > /dev/null)
res=$(awk '/direct : profile/ { for (i = 1; i <= NF; i++) if ($i == "bandwidth") bw = $(i+2) }
	END { if (bw == "") { printf "NG (no report)" }
	      else printf "%s (mean bandwidth %.1f, dense would be 816)",
	                  ((bw < 1632 / 4) ? "OK" : "NG"), bw }' "$WORK/ofe.log")
echo "  RCM reduces the profile : $res"
case "$res" in NG*) status=1 ;; esac

# 複素対称系 (F / A) の直接解法 (スカイライン **LDL^T**)。
#
# 実対称のときと違い rlc.csv はバイト単位では一致しない。この系は
# ‖A‖‖x‖ >> ‖b‖ (低周波では ωM << K) なので、後退安定な分解でも相対残差は
# 1e-11 程度で頭打ちになり、COCG とは最後の 1 桁が違い得る。そこで
# 「同じ答えに落ちること」を相対誤差で見て、**分解が破綻していないことは
# ログの残差で別に assert する** (値の一致だけだと、両方が同じくらい
# 間違っている可能性を排除できない)
echo "[direct F/A] the complex-symmetric (LDL^T) direct solver"
direct_close() {	# direct_close <ラベル> <反復の .ofe> <直接の .ofe> <許容>
	if ! (cd "$WORK" && "$OFE" -n 2 "$2" > /dev/null && "$OFE_POST" > /dev/null); then
		echo "  $1 : the iterative run failed -> NG" >&2
		status=1
		return
	fi
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/it.csv"
	if ! (cd "$WORK" && "$OFE" -n 2 "$3" > /dev/null && "$OFE_POST" > /dev/null); then
		echo "  $1 : the direct run failed -> NG" >&2
		status=1
		return
	fi
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/di.csv"
	if ! grep -q "direct : profile" "$WORK/ofe.log"; then
		echo "  $1 : the direct solver was not used -> NG" >&2
		status=1
		return
	fi
	# 分解の残差 (ソルバーが自分で計算してログに出したもの)。黙って反復解法に
	# 落ちていれば行そのものが無いので、そこも NG にする
	res=$(awk '$NF == "direct" { r = $(NF - 1) + 0; if (r > m) m = r; n++ }
		END { if (n == 0) printf "NG (no residual reported)"
		      else printf "%s (max %.2e over %d solves)", ((m < 1e-8) ? "OK" : "NG"), m, n }' \
		"$WORK/ofe.log")
	echo "  $1 : factorization residual $res"
	case "$res" in NG*) status=1 ;; esac
	# 反復解法との一致 (rlc.csv の全成分の最大相対差)
	res=$(awk -F, -v tol="$4" '
		NR == FNR { a[FNR] = $0; next }
		{ n = split(a[FNR], p, ",")
		  if (n != NF) { bad = 1; exit }
		  for (i = 1; i <= NF; i++) {
			if (($i + 0 != $i) || (p[i] + 0 != p[i])) {
				if ($i != p[i]) bad = 1		# 名前の行が食い違った
				continue
			}
			d = $i - p[i]; if (d < 0) d = -d
			s = ($i < 0) ? -$i : $i
			t = (p[i] < 0) ? -p[i] : p[i]
			if (t > s) s = t
			if ((s > 0) && ((d / s) > m)) m = d / s
		  }
		}
		END { if (bad || (FNR == 0)) printf "NG (the two outputs differ in shape)"
		      else printf "%s (max %.1e)", ((m <= tol) ? "OK" : "NG"), m }' \
		"$WORK/it.csv" "$WORK/di.csv")
	echo "  $1 : agrees with the iterative solver $res"
	case "$res" in NG*) status=1 ;; esac
}
# (a) 構造格子の断面 2 次元 (F、表皮効果)
cp "$SRC/plate_line_ac.ofe" "$WORK/"
sed 's/^analysis = /direct = 1\nanalysis = /' "$SRC/plate_line_ac.ofe" > "$WORK/dir_f.ofe"
direct_close "plate_line_ac (structured, F)" plate_line_ac.ofe dir_f.ofe 1e-6
# (b) 非構造格子 (三角形) の F。断面 2 次元の CRS パターンを通す
cp "$SRC/plate2d_ac.ofe" "$WORK/"
sed 's/^analysis = /direct = 1\nanalysis = /' "$SRC/plate2d_ac.ofe" > "$WORK/dir_f2.ofe"
direct_close "plate2d_ac (triangles, F)" plate2d_ac.ofe dir_f2.ofe 1e-6
# (c) 3 次元 A-φ 連成 (辺 + 節点の混在した行列)。**ゲージ固定が要る** ので
#     反復側も gauge = 1 で回し、違いが解法だけになるようにする
sed 's/^analysis = /gauge = 1\nanalysis = /' "$SRC/bar_eddy.ofe" > "$WORK/g1_a.ofe"
sed 's/^analysis = /gauge = 1\ndirect = 1\nanalysis = /' "$SRC/bar_eddy.ofe" > "$WORK/dir_a.ofe"
direct_close "bar_eddy (3-D A-phi, gauge = 1)" g1_a.ofe dir_a.ofe 1e-6

# ゲージを固定しない A-φ は (A, φ) -> (A + Gψ, φ - jωψ) で不変な**特異系**。
# COCG は右辺が値域に入るので収束するが、分解はピボットが 0 になって成立しない。
# 黙って壊れた答えを返さずに入力段で弾くこと
sed 's/^analysis = /direct = 1\nanalysis = /' "$SRC/bar_eddy.ofe" > "$WORK/dir_a0.ofe"
mesh_reject "direct = 1 with analysis A and no gauge" dir_a0.ofe
if ! (cd "$WORK" && "$OFE" -n 2 dir_a0.ofe 2>&1 | grep -q "requires gauge = 1"); then
	echo "  *** it was rejected for a different reason than the missing gauge" >&2
	status=1
fi

# Gmsh ASCII 4.1 形式の読み込み。**同じ形状を 2.2 と 4.1 で書いた結果が
# 完全に一致すること**を恒等式に使う (閉形式は要らないうえ、これ以上厳密な
# 検査は無い)。4.1 は要素の行に物理タグが無く $Entities に付くので、
# そこを読み落とすと全要素のタグが 0 になって材料も電極も割り当たらない
echo "[gmsh 4.1] the same mesh in 2.2 and 4.1 must give identical results"
same_mesh() {	# same_mesh <ラベル> <ofe(2.2)> <ofe(4.1)>
	run_case "$2"; cp "$WORK/rlc.csv" "$WORK/a.csv"
	run_case "$3"; cp "$WORK/rlc.csv" "$WORK/b.csv"
	# title 行だけはファイル名由来で違うので外す
	grep -v '^title' "$WORK/a.csv" > "$WORK/a2.csv"
	grep -v '^title' "$WORK/b.csv" > "$WORK/b2.csv"
	if cmp -s "$WORK/a2.csv" "$WORK/b2.csv"; then
		echo "  $1 : identical -> OK"
	else
		echo "  $1 : differ -> NG" >&2
		diff "$WORK/a2.csv" "$WORK/b2.csv" | head -4 >&2
		status=1
	fi
}
same_mesh "3-D tetrahedral mesh (box_tet)" box_tet box_tet_41
same_mesh "2-D triangular mesh (plate2d)" plate2d_dc plate2d_41
# $Entities を落とすと物理タグが全部 0 になる。黙って解かずに落ちること
awk '/^\$Entities/ { skip = 1 } /^\$EndEntities/ { skip = 0; next } !skip { print }' \
    "$SRC/box_tet_41.msh" > "$WORK/noent.msh"
sed 's/^mesh = .*/mesh = noent.msh/' "$SRC/box_tet_41.ofe" > "$WORK/noent.ofe"
mesh_reject "Gmsh 4.1 without \$Entities" noent.ofe
# ヘッダだけバイナリを名乗って中身が ASCII のファイルは弾くこと。
# バイナリを読めるようになった今でも、これは黙って誤読してはいけない形
sed 's/^4.1 0 8$/4.1 1 8/' "$SRC/box_tet_41.msh" > "$WORK/bin.msh"
sed 's/^mesh = .*/mesh = bin.msh/' "$SRC/box_tet_41.ofe" > "$WORK/bin.ofe"
mesh_reject "ASCII content with a binary header" bin.ofe
if ! (cd "$WORK" && "$OFE" -n 2 bin.ofe 2>&1 | grep -q "content is ASCII"); then
	echo "  *** it was rejected for a different reason than the false header" >&2
	status=1
fi

# 六面体 (8 節点、等パラメトリック) の非構造格子。
#
# 検証は 3 段に分ける。**どの検査がどの誤りで落ちるか**が違うため:
#   (a) 直方体   : 閉形式と全桁一致 (要素・組み立て・電極が通っていること)
#   (b) 剛体回転 : (a) と全桁一致。回すと軸に平行でなくなるので、ヤコビアンの
#                  逆行列の転置はここでだけ落ちる (実測 -17.7%)
#   (c) ゆがみ   : analysis = P の線形場の恒等式。**平行六面体では
#                  ヤコビアンが要素内で一定**になるので、(a) も (b) も
#                  「J を要素中心で 1 回だけ評価する」誤りを検出できない
#                  (実測: どちらも 1 桁も動かない)。ゆがんだ格子でだけ 4.7% ずれる
#   (d) 同軸     : 曲がった (アフィンでない) 形状で閉形式と比較。ただし回転対称
#                  なので (c) の誤りはここでも検出できない (要素行列が 3% 違っても
#                  解が theta に依らず相殺する)
echo "[hexahedra] 8-node isoparametric hexahedral meshes"
run_case box_hex
compare "C (hex box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
grep -v '^title' "$WORK/rlc.csv" > "$WORK/hexa.csv"
# 剛体回転しても合同な離散問題なので全桁一致すること
sed 's/^mesh = .*/mesh = box_hex_rot.msh/' "$SRC/box_hex.ofe" > "$WORK/hexrot.ofe"
(cd "$WORK" && "$OFE" -n 2 hexrot.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/hexb.csv"
if cmp -s "$WORK/hexa.csv" "$WORK/hexb.csv"; then
	echo "  rigidly rotated mesh gives the identical answer : OK"
else
	echo "  rigidly rotated mesh differs -> NG" >&2
	diff "$WORK/hexa.csv" "$WORK/hexb.csv" | head -4 >&2
	status=1
fi
# ゆがんだ格子での自己検証 (analysis = P)。恒等式は機械精度で成り立つ
cp "$SRC/nodal_test_hex.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 nodal_test_hex.ofe > /dev/null)
res=$(awk '/linear field/ { for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
	/rel. diff/ { for (i = 1; i <= NF; i++) if ($i == "diff") v = $(i+2) }
	/^  warp/ { w = $3 }
	END { if ((e == "") || (v == "")) { printf "NG (no report)"; exit }
	      printf "%s (linear %s, volume %s, warp %s)",
	             (((e + 0) < 1e-12) && ((v + 0) < 1e-12) && ((w + 0) > 1e-3)) ? "OK" : "NG",
	             e, v, w }' "$WORK/ofe.log")
echo "  distorted mesh : the linear-field identity holds : $res"
case "$res" in NG*) status=1 ;; esac
# 曲がった (アフィンでない) 六面体を閉形式と比べる
run_case coax_hex
compare "C' (hex coax) [F/m]" "$(value_of C)" 1.063417e-10 0.01
compare "L' (hex coax) [H/m]" "$(value_of L)" 2.197225e-07 0.01
# 四面体と六面体の混在は弾くこと (要素種別で分岐しているので黙って通してはいけない)。
# 六面体の格子に四面体を 1 個だけ足す
awk '/^\$Elements/ { print; getline; print $1 + 1
                     print "9999 4 2 1 1 1 2 3 4"; next } { print }' \
    "$SRC/box_hex.msh" > "$WORK/mix.msh"
sed 's/^mesh = .*/mesh = mix.msh/' "$SRC/box_hex.ofe" > "$WORK/mix.ofe"
mesh_reject "a mesh mixing tetrahedra and hexahedra" mix.ofe
if ! (cd "$WORK" && "$OFE" -n 2 mix.ofe 2>&1 | grep -q "conformingly"); then
	echo "  *** it was rejected for a different reason than the conformity" >&2
	status=1
fi
# analysis は 1 トークン 1 文字。"CL" と続けて書くと先頭しか読まれないので弾く
sed 's/^analysis = .*/analysis = CL/' "$SRC/box_hex.ofe" > "$WORK/an.ofe"
mesh_reject "analysis = CL (letters not separated)" an.ofe
if ! (cd "$WORK" && "$OFE" -n 2 an.ofe 2>&1 | grep -q "separated by spaces"); then
	echo "  *** it was rejected for a different reason than the analysis spelling" >&2
	status=1
fi
# 辺要素 (E / A) は四面体の形状関数に基づくので六面体格子を弾くこと
sed -e 's/^analysis = .*/analysis = A/' \
    -e 's/^material = .*/material = 1.0 5.8e7/' \
    -e 's/^solver = .*/frequency = 1e4\nvoltage = 1.0/' "$SRC/box_hex.ofe" > "$WORK/hexedge.ofe"
mesh_reject "analysis A on a hexahedral mesh" hexedge.ofe
if ! (cd "$WORK" && "$OFE" -n 2 hexedge.ofe 2>&1 | grep -q "need a tetrahedral mesh"); then
	echo "  *** it was rejected for a different reason than the element type" >&2
	status=1
fi

# 角柱 (6 節点、三角形を押し出した等パラメトリック要素)。
# 検査の分け方は六面体と同じで、**ゆがんだ格子でだけ落ちる誤り**があるのが要点
echo "[prisms] 6-node isoparametric prism (wedge) meshes"
run_case box_prism
compare "C (prism box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
grep -v '^title' "$WORK/rlc.csv" > "$WORK/pra.csv"
sed 's/^mesh = .*/mesh = box_prism_rot.msh/' "$SRC/box_prism.ofe" > "$WORK/prrot.ofe"
(cd "$WORK" && "$OFE" -n 2 prrot.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/prb.csv"
if cmp -s "$WORK/pra.csv" "$WORK/prb.csv"; then
	echo "  rigidly rotated mesh gives the identical answer : OK"
else
	echo "  rigidly rotated mesh differs -> NG" >&2
	status=1
fi
cp "$SRC/nodal_test_prism.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 nodal_test_prism.ofe > /dev/null)
res=$(awk '/linear field/ { for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
	/rel. diff/ { for (i = 1; i <= NF; i++) if ($i == "diff") v = $(i+2) }
	/^  warp/ { w = $3 }
	END { if ((e == "") || (v == "")) { printf "NG (no report)"; exit }
	      printf "%s (linear %s, volume %s, warp %s)",
	             (((e + 0) < 1e-12) && ((v + 0) < 1e-12) && ((w + 0) > 1e-3)) ? "OK" : "NG",
	             e, v, w }' "$WORK/ofe.log")
echo "  distorted mesh : the linear-field identity holds : $res"
case "$res" in NG*) status=1 ;; esac

# 要素種別の混在 (角柱 + 四面体)。境界層を角柱、内部を四面体にした形で、
# 面は三角形どうしなので**適合する**。六面体と四面体の直接の隣接だけは
# 面上の解が食い違うので弾く (上の [hexahedra] 節で検査している)
echo "[mixed] mixing element types in one mesh"
run_case box_mixed
# 角柱側 er = 4、四面体側 er = 2 の直列。材料を種別ごとに変えてあるので、
# 要素番号 -> 材料の対応が種別の境目でずれると必ず落ちる
compare "C (prism er=4 + tet er=2) [F]" "$(value_of C)" 1.180558376e-13 1e-8
# 混在で一番危ないのは「種別をまたぐ面で節点が共有されない」形。
# 線形場の恒等式はそれが起きると破れるので、そこを見る
cp "$SRC/nodal_test_mixed.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 nodal_test_mixed.ofe > /dev/null)
res=$(awk '/linear field/ { for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
	/^Nodes =/ { for (i = 1; i <= NF; i++) {
	                 if ($i == "Tetrahedra") nt = $(i+2)+0
	                 if ($i == "Prisms") np = $(i+2)+0 } }
	END { if (e == "") { printf "NG (no report)"; exit }
	      printf "%s (%s, %d tets + %d prisms)",
	             (((e + 0) < 1e-12) && (nt > 0) && (np > 0)) ? "OK" : "NG", e, nt, np }' \
	"$WORK/ofe.log")
echo "  the linear-field identity holds across the interface : $res"
case "$res" in NG*) status=1 ;; esac

# Gmsh **バイナリ**形式の読み込み。
#
# 検証用のファイルは**本物の gmsh に作らせてある** (4.12.1):
#   gmsh -0 <格子> -o <出力> -format msh22|msh41 [-bin]
# 自作の書き手と読み手が同じ誤解を共有していると、テストは通るのに実際の
# ツールが出したファイルで落ちる。形式の検証では出所が本質なので、
# ここだけは自前で書き出したファイルを使わない。
#
# gmsh は変換のたびに節点を振り直すため、比較する ASCII 側も**同じ変換で
# 出したもの** (box_bin.msh 等) を置いてある。元の mkmesh.py 出力と直接
# 比べると節点の並びが違い、丸めの順序が変わって一致しない。
echo "[gmsh binary] the same mesh in ASCII and binary must give identical results"
for m in "$SRC"/*.msh; do
	[ -f "$m" ] && cp "$m" "$WORK/"
done
bin_same() {	# bin_same <ラベル> <.ofe> <ASCII の .msh> <バイナリの .msh>
	# **本当にバイナリかを確かめる。** ASCII に差し替えても一致検査は通るので、
	# それだけでは「バイナリを読めている」ことの証拠にならない
	if ! head -c 24 "$WORK/$4" | tr '\n' ' ' | grep -qa 'MeshFormat [0-9.]* 1 8'; then
		echo "  $1 : $4 is not a binary Gmsh file -> NG" >&2
		status=1
		return
	fi
	sed "s/^mesh = .*/mesh = $3/" "$SRC/$2" > "$WORK/bm_a.ofe"
	sed "s/^mesh = .*/mesh = $4/" "$SRC/$2" > "$WORK/bm_b.ofe"
	if ! (cd "$WORK" && "$OFE" -n 2 bm_a.ofe > /dev/null && "$OFE_POST" > /dev/null); then
		echo "  $1 : the ASCII run failed -> NG" >&2
		status=1
		return
	fi
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/bm_a.csv"
	if ! (cd "$WORK" && "$OFE" -n 2 bm_b.ofe > /dev/null && "$OFE_POST" > /dev/null); then
		echo "  $1 : the binary run failed -> NG" >&2
		status=1
		return
	fi
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/bm_b.csv"
	if cmp -s "$WORK/bm_a.csv" "$WORK/bm_b.csv"; then
		echo "  $1 : identical -> OK"
	else
		echo "  $1 : differ -> NG" >&2
		diff "$WORK/bm_a.csv" "$WORK/bm_b.csv" | head -4 >&2
		status=1
	fi
}
# 3 次元四面体 / 断面 2 次元三角形 / 2 次要素の 3 種を、2.2 と 4.1 の両方で。
# 2 次要素を入れてあるのは、バイナリでは節点数が型から決まる (10 / 6) ためで、
# 型ごとの節点数表を間違えるとここでずれる
bin_same "3-D tetrahedra, binary 2.2" box_tet.ofe box_bin.msh box_bin_22.msh
bin_same "3-D tetrahedra, binary 4.1" box_tet.ofe box_bin.msh box_bin_41.msh
bin_same "2-D triangles, binary 2.2" plate2d_dc.ofe plate2d_bin.msh plate2d_bin_22.msh
bin_same "2-D triangles, binary 4.1" plate2d_dc.ofe plate2d_bin.msh plate2d_bin_41.msh
bin_same "order-2 tetrahedra, binary 2.2" box_p2.ofe box_p2_bin.msh box_p2_bin_22.msh
bin_same "order-2 tetrahedra, binary 4.1" box_p2.ofe box_p2_bin.msh box_p2_bin_41.msh

# 逆エンディアンのファイルは**読み違えずに落ちること**。バイト入れ替えは
# 手元で検証できないので実装せず、はっきり断る方を選んでいる
sed 's/^mesh = .*/mesh = bin_swapped.msh/' "$SRC/box_tet.ofe" > "$WORK/binsw.ofe"
mesh_reject "binary Gmsh with the opposite byte order" binsw.ofe
if ! (cd "$WORK" && "$OFE" -n 2 binsw.ofe 2>&1 | grep -q "opposite byte order"); then
	echo "  *** it was rejected for a different reason than the byte order" >&2
	status=1
fi
# **知らない要素型はバイナリでは読み飛ばせない** (1 要素 1 行の ASCII と違い、
# 節点数が分からないと何バイト進めばよいか決まらない)。黙って誤読しないこと
sed 's/^mesh = .*/mesh = bin_unknown_type.msh/' "$SRC/box_tet.ofe" > "$WORK/binut.ofe"
mesh_reject "binary Gmsh with an unknown element type" binut.ofe
if ! (cd "$WORK" && "$OFE" -n 2 binut.ofe 2>&1 | grep -q "unknown element type"); then
	echo "  *** it was rejected for a different reason than the element type" >&2
	status=1
fi

# 断面 2 次元の 2 次要素 (6 節点三角形)。**導体内の Az は電流密度が一様なとき
# 厳密に 2 次**なので、内部インダクタンス 2t/3 を P2 は厳密に表せる。
# 1 次要素は -0.045% ずれるので、次数が効いていることがそのまま出る
echo "[plate2d_p2] second-order (6-node) triangles on the 2-D mesh"
run_case plate2d_p2
compare "Ldc [H/m]" "$(value_of Ldc)" 2.932153e-07 1e-5
compare "Rs [ohm/m]" "$(value_of Rs)" 6.896552e-01 1e-6
# 1 次との対比 (2 次が黙って 1 次に退行していないこと)
lp2=$(value_of Ldc)
run_case plate2d_dc
res=$(awk -v a="$(value_of Ldc)" -v b="$lp2" \
	'BEGIN{ e = 2.932153e-07
	        e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	        e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	        printf "P1 %+.4f%% vs P2 %+.4f%% -> %s", 100 * e1, 100 * e2,
	               ((e1 > 100 * e2) ? "OK" : "NG") }')
echo "  order 2 beats order 1 : $res"
case "$res" in *NG) status=1 ;; esac

echo "[plate2d_p2 F] the same mesh with the eddy-current analysis"
for pair in "1e3 6.896552e-01 2.932153e-07" "1e7 1.624296e+00 2.780550e-07"; do
	set -- $pair
	sed -e 's/plate2d.msh/plate2d_p2.msh/' -e "s/^frequency = .*/frequency = $1/" \
	    "$SRC/plate2d_ac.ofe" > "$WORK/p2ac.ofe"
	(cd "$WORK" && "$OFE" -n 2 p2ac.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R(f=$1) [ohm/m]" "$(value_of Rf)" "$2" 1e-4
	compare "L(f=$1) [H/m]" "$(value_of Lf)" "$3" 1e-4
done
# 場の出力 : 2 次三角形でも恒等式が成り立つこと + セル型が QUADRATIC_TRIANGLE
sed -e 's/plate2d.msh/plate2d_p2.msh/' -e 's/^analysis = F/fieldout = 1\nanalysis = F/' \
    "$SRC/plate2d_ac.ofe" > "$WORK/p2fld.ofe"
(cd "$WORK" && "$OFE" -n 2 p2fld.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "field area (P2 2-D) [m^2]" "$(vtk vol J_F_re_port1)" 3.0e-7 1e-9
pf=$(awk -v a="$(vtk int2 J_F_re_port1)" -v b="$(vtk int2 J_F_im_port1)" \
	'BEGIN{ printf "%.10e", (a + b) / (2 * 5.8e7) }')
pt=$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
	        printf "%.10e", 0.5 * R / ((R * R) + (X * X)) }')
compare "ohmic loss from the field (P2 2-D) [W/m]" "$pf" "$pt" 1e-4
res=$(awk '
	NF == 0       { next }
	/^CELLS/      { st = "c"; k = 0; next }
	/^CELL_TYPES/ { st = "t"; k = 0; next }
	/^[A-Z_]+ /   { st = ""; next }
	st == "t" { k++; if ($1 != 22) type++; next }
	st == "c" { if ($1 != 6) nn++; next }
	END { if (k == 0) { printf "NG (no cells)" }
	      else if (type || nn) { printf "NG (%d not type 22, %d not 6-node)", type, nn }
	      else { printf "OK" } }' "$WORK/ofe_field.vtk")
echo "  VTK cell type is QUADRATIC_TRIANGLE with 6 nodes : $res"
case "$res" in NG*) status=1 ;; esac
# 次数の混在は弾くこと
awk '/^[0-9]+ 9 2 / && !done { print $1, 2, 2, $4, $5, $6, $7, $8; done = 1; next } { print }' \
    "$SRC/plate2d_p2.msh" > "$WORK/mix2d.msh"
sed 's/^mesh = .*/mesh = mix2d.msh/' "$SRC/plate2d_p2.ofe" > "$WORK/p2mix.ofe"
mesh_reject "mixed triangle orders in a 2-D mesh" p2mix.ofe

# 異方性の μ。**B は ∇Az を 90 度回したものなので、ν を ∇Az にそのまま掛けると
# 面内 2 成分を取り違える。** 等方性では一致するので等方性ケースでは検出できない。
# 1 次元解では B が p 軸を向くので、効くのは μ_pp だけ (μ_qq は L を変えない)。
#   構造格子 : tline = Z -> p = x   (plate_line_dc の形状)
#   2 次元格子 : tline = X -> p = y (plate2d の形状)
echo "[plate2d_rot] the same cross-section rotated 90 degrees in the plane"
run_case plate2d_rot
# 面内で回しただけなので等方性の答えは変わらない (回転不変性)
compare "rotated mesh, isotropic Ldc [H/m]" "$(value_of Ldc)" 2.932153e-07 0.002

echo "[anisotropic mu] the in-plane tensor must act on B, not on grad(Az)"
aniso_mu() {	# aniso_mu <ラベル> <ofe> <sed で入れる anisomur 行> <期待値>
	sed "s/^region = 1 2/region = 1 2\n$3/" "$SRC/$2.ofe" > "$WORK/amu.ofe"
	(cd "$WORK" && "$OFE" -n 2 amu.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "$1" "$(value_of Ldc)" "$4" 0.002
}
# L'dc = mu0 (mur_pp d + 2t/3) / W
aniso_mu "2-D mesh, mu_yy = 5 (along B)"     plate2d_dc "anisomur = 2 1.0 5.0 1.0" 1.298525e-06
aniso_mu "2-D mesh, mu_zz = 5 (across B)"    plate2d_dc "anisomur = 2 1.0 1.0 5.0" 2.932153e-07
# **面内で 90 度回した格子も回すこと。** 1 次元解では面内勾配の片方しか立たない
# ので、片方の格子だけだと面内テンソルの一方の成分が一度も検査されない
# (実測: c_pp を壊す変異が回転前の格子だけでは素通りした)
aniso_mu "rotated mesh, mu_zz = 5 (along B)"  plate2d_rot "anisomur = 2 1.0 1.0 5.0" 1.298525e-06
aniso_mu "rotated mesh, mu_yy = 5 (across B)" plate2d_rot "anisomur = 2 1.0 5.0 1.0" 2.932153e-07

# 非対角成分 (ν_pq) はここでしか検査できない。
#
# 1 次元解では ∂Az/∂p = 0 なので非対角項が効かず、同軸は回転対称なので
# テンソルを回しても答えが変わらない (符号を反転してもビット単位で同じ結果に
# なることを確認済み)。**斜めに回した非対称断面**だけが効く。
#
# plate2d_r30.msh は plate2d.msh を面内に 30 度回しただけで位相は同一なので、
# 材料テンソルも同じ 30 度回せば**元と合同な離散問題**になり、答えは一致する
# はずである (実測で全桁一致)。これが恒等式になる。
#
#   mu = diag(5, 1) を 30 度回す:
#     mu_yy = 5c^2 + 1s^2 = 4,  mu_zz = 5s^2 + 1c^2 = 2,
#     mu_yz = (5-1) c s = 1.7320508075688772
#
# 非対角を落とすと L が -32.0%、符号を反転すると回転でなく鏡映になるので
# やはりずれる
echo "[off-diagonal mu] rotate the mesh and the tensor by the same angle"
run_case plate2d_r30
compare "rotated 30deg, isotropic Ldc [H/m]" "$(value_of Ldc)" 2.932153e-07 0.002
aniso_mu "rotated 30deg mesh + tensor" plate2d_r30 \
	"anisomur = 2 1.0 4.0 2.0 0.0 1.732050807568877 0.0" 1.298525e-06
# 対照 : 非対角を落とすと合同でなくなるので必ずずれること
sed 's/^region = 1 2/region = 1 2\nanisomur = 2 1.0 4.0 2.0 0.0 0.0 0.0/' \
    "$SRC/plate2d_r30.ofe" > "$WORK/amu0.ofe"
(cd "$WORK" && "$OFE" -n 2 amu0.ofe > /dev/null && "$OFE_POST" > /dev/null)
res=$(awk -v v="$(value_of Ldc)" -v e=1.298525e-06 \
	'BEGIN{ d = (v - e) / e; if (d < 0) d = -d
	        printf "%s (%.1f%% off)", ((d > 0.1) ? "OK" : "NG"), 100 * d }')
echo "  dropping mu_yz must change the answer : $res"
case "$res" in NG*) status=1 ;; esac
# 構造格子側は geometry で間隙に材料 2 を割り当てる
mkaniso() {
	sed -e "s/^conductorsigma = 0 5.8e7/material = 1.0 0\ngeometry = 2 1 0 1e-3 0.05e-3 0.25e-3 0 1e-4\n$1\nconductorsigma = 0 5.8e7/" \
	    "$SRC/plate_line_dc.ofe" > "$WORK/samu.ofe"
	(cd "$WORK" && "$OFE" -n 2 samu.ofe > /dev/null && "$OFE_POST" > /dev/null)
}
mkaniso "anisomur = 2 5.0 1.0 1.0"
compare "structured, mu_xx = 5 (along B)"  "$(value_of Ldc)" 1.298525e-06 0.002
mkaniso "anisomur = 2 1.0 5.0 1.0"
compare "structured, mu_yy = 5 (across B)" "$(value_of Ldc)" 2.932153e-07 0.002

# 2 次元格子で弾くべき入力
echo "[2-D mesh] inputs the solver must reject"
sed 's/^analysis = M/analysis = C/' "$SRC/plate2d_dc.ofe" > "$WORK/t2d_c.ofe"
mesh_reject "analysis C on a 2-D mesh" t2d_c.ofe
sed 's/^tline = X/tline = Y/' "$SRC/plate2d_dc.ofe" > "$WORK/t2d_plane.ofe"
mesh_reject "2-D mesh not normal to tline" t2d_plane.ofe
sed 's/^region = 1 2/region = 1 2\nbh = 2 100 0.5\nbh = 2 1000 1.5/' \
    "$SRC/plate2d_dc.ofe" > "$WORK/t2d_bh.ofe"
mesh_reject "nonlinear (bh) on a 2-D mesh" t2d_bh.ofe

# ofe_post が掃引を認識すること。rlc.csv は掃引の**最後の 1 点**しか持たないので、
# 黙って出すと「掃引の結果」と誤読される。点数・範囲・注記を出すこと。
# **誤検知も同じくらい悪い**ので、掃引していないケースで出ないことも見る
echo "[post sweep] rlc.csv must say it holds only the last sweep point"
run_case sweep_plate
if grep -q "^frequency sweep points,5$" "$WORK/rlc.csv" \
   && grep -q "^sweep range \[Hz\],1.00000000e+08,1.00000000e+10$" "$WORK/rlc.csv" \
   && grep -q "LAST sweep point" "$WORK/rlc.csv"; then
	echo "  sweep note in rlc.csv : OK"
else
	echo "  sweep note in rlc.csv : NG" >&2
	status=1
fi
run_case parallel_plate
if grep -q "sweep" "$WORK/rlc.csv"; then
	echo "  no sweep note without a sweep : NG" >&2
	status=1
else
	echo "  no sweep note without a sweep : OK"
fi

# Havriliak-Negami 分散。**3 つの極限がそのまま恒等式になる**ので、
# 閉形式を別に用意しなくても厳密に検査できる:
#   β = 1        -> Cole-Cole (同じ α)
#   α = 0        -> Cole-Davidson (同じ β)
#   α = 0, β = 1 -> Debye
# 「別のキーで書いた同じ物理」が全桁一致することを見る
echo "[hn_plate] Havriliak-Negami (alpha = 0.3, beta = 0.6) at 1 GHz"
run_case hn_plate
compare "C [F]" "$(value_of C)" 1.797624606e-13 1e-6
compare "G [S]" "$(value_of G)" 1.962350996e-04 1e-6
hn_same() {	# hn_same <ラベル> <極 1 の行> <極 2 の行>
	sed "s/^havriliak = .*/$2/" "$SRC/hn_plate.ofe" > "$WORK/hn1.ofe"
	(cd "$WORK" && "$OFE" -n 2 hn1.ofe > /dev/null && "$OFE_POST" > /dev/null)
	c1=$(value_of C); g1=$(value_of G)
	sed "s/^havriliak = .*/$3/" "$SRC/hn_plate.ofe" > "$WORK/hn2.ofe"
	(cd "$WORK" && "$OFE" -n 2 hn2.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "$1 : C" "$(value_of C)" "$c1" 1e-9
	compare "$1 : G" "$(value_of G)" "$g1" 1e-9
}
hn_same "beta = 1 equals Cole-Cole" \
	"havriliak = 2 2.0 3.0 1.591549431e-10 0.3 1.0" \
	"colecole = 2 2.0 3.0 1.591549431e-10 0.3"
hn_same "alpha = 0 equals Cole-Davidson" \
	"havriliak = 2 2.0 3.0 1.591549431e-10 0.0 0.6" \
	"coledavidson = 2 2.0 3.0 1.591549431e-10 0.6"
hn_same "alpha = 0, beta = 1 equals Debye" \
	"havriliak = 2 2.0 3.0 1.591549431e-10 0.0 1.0" \
	"debye = 2 2.0 3.0 1.591549431e-10"
# beta の範囲外は弾くこと
sed 's/^havriliak = .*/havriliak = 2 2.0 3.0 1.591549431e-10 0.3 1.5/' \
    "$SRC/hn_plate.ofe" > "$WORK/hn_bad.ofe"
mesh_reject "havriliak with beta > 1" hn_bad.ofe

# εr / μr / B-H の温度依存。σ(T) と同じ 1 次係数だが、**適用箇所が違う**:
#   εr / μr は読み出し時 (material_freq が εr を上書きするので入力時では消える)
#   B-H は入力時 (曲線は再計算されないので σ と同じでよい)
echo "[temp_material] epsr(T) = epsr0 (1 + alpha (T - T0)), C is exactly proportional"
for pair in "20 1.7708376e-13" "70 2.4849e-13" "120 2.2135470e-13" "-30 1.5495329e-13"; do
	set -- $pair
	sed "s/^temperature = .*/temperature = $1/" "$SRC/temp_material.ofe" > "$WORK/tm.ofe"
	(cd "$WORK" && "$OFE" -n 2 tm.ofe > /dev/null && "$OFE_POST" > /dev/null)
	# C = eps0 epsr0 (1 + alpha (T - T0)) A / d
	exp=$(awk -v t="$1" 'BEGIN{ printf "%.10e", 1.7708375600e-13 * (1 + 2.5e-3 * (t - 20)) }')
	compare "C(T=$1) [F]" "$(value_of C)" "$exp" 1e-6
done

# μr の温度依存。**L 全体は μr に比例しない** (内部インダクタンス 2t/3 は
# 導体側なので温度で動かない)。比例と誤って実装すると必ずずれる
echo "[temp_mur] mur(T) with the internal inductance held fixed"
for t in 20 120; do
	sed "s/^temperature = .*/temperature = $t/" "$SRC/temp_mur.ofe" > "$WORK/tu.ofe"
	(cd "$WORK" && "$OFE" -n 2 tu.ofe > /dev/null && "$OFE_POST" > /dev/null)
	# L' = mu0 (mur(T) d + 2t/3) / W
	exp=$(awk -v t="$t" 'BEGIN{ pi = 3.14159265358979324; mu0 = 4 * pi * 1e-7
		mur = 4.0 * (1 - 3.0e-3 * (t - 20))
		printf "%.10e", mu0 * ((mur * 0.2e-3) + (2 * 0.05e-3 / 3)) / 1e-3 }')
	compare "Ldc(T=$t) [H/m]" "$(value_of Ldc)" "$exp" 0.002
done

# B-H 曲線の温度依存。1 次元では H = I/W が Ampere の法則で固定されるので、
# B 軸を k 倍すると**間隙の鎖交磁束がちょうど k 倍**になる。内部インダクタンス
# L_int = mu0 2t/(3W) は動かないので (L - L_int) の比が厳密に k になる
echo "[bhtempco] scaling the B axis scales the flux exactly"
lint_bh=""
for t in 20 120; do
	sed 's/^analysis = M/bhtempco = 2 -2.0e-3 20\ntemperature = '"$t"'\nanalysis = M/' \
	    "$SRC/plate_line_bh.ofe" > "$WORK/tb.ofe"
	(cd "$WORK" && "$OFE" -n 2 tb.ofe > /dev/null && "$OFE_POST" > /dev/null)
	lint_bh="$lint_bh $(value_of Ldc)"
done
set -- $lint_bh
compare "(L - L_int) ratio at T = 120 vs 20" \
	"$(awk -v a="$2" -v b="$1" 'BEGIN{ pi = 3.14159265358979324; mu0 = 4 * pi * 1e-7
		li = mu0 * 2 * 0.05e-3 / (3 * 1e-3)
		printf "%.10e", (a - li) / (b - li) }')" 0.8 1e-5

# 鉄損 (Bertotti の損失分離)。板間を積層鉄心で埋めた平行平板線路で、
# 磁束密度が厳密に 1 次元 (B = mu0 I/W) なので閉形式と機械精度で一致する
echo "[bertotti_core] iron loss vs the closed form on a 1-D exact field"
pfe_of() {	# pfe_of : rlc.csv の Pfe (ポート 1)
	awk -F, '$1 == "Pfe" { f = 1; next } f && ($1 == "1") { print $2; exit }' "$WORK/rlc.csv"
}
run_case bertotti_core
# B = mu0 I/W、V = d W (単位長あたり)
exp_all=$(awk 'BEGIN{ pi = 3.14159265358979324; mu0 = 4 * pi * 1e-7
	B = mu0 * 1 / 1e-3; V = 0.2e-3 * 1e-3; f = 50
	ph = 200 * f * (B ^ 2) * V
	pc = (pi * pi * 2e6 * (0.35e-3) ^ 2 / 6) * f * f * B * B * V
	pe = 1.5 * ((f * B) ^ 1.5) * V
	printf "%.10e", ph + pc + pe }')
compare "Pfe (total) [W/m]" "$(pfe_of)" "$exp_all" 1e-9
compare "max|B| in the core [T]" \
	"$(awk '/max\|B\|/ { for (i = 1; i <= NF; i++) if ($i == "max|B|") print $(i+2) }' "$WORK/ofe.log")" \
	1.2566370614e-3 1e-6

# **3 項の分離は指数で検査する。** 係数だけ合わせても指数を取り違えると
# 別の周波数・電流で必ずずれるので、1 項ずつ残して f と I を 2 倍にする:
#   ヒステリシス f^1   B^alpha(=2) -> f x2、I x4
#   古典渦電流   f^2   B^2         -> f x4、I x4
#   異常 (過剰)  f^1.5 B^1.5       -> f x2.8284、I x2.8284
bert_ratio() {	# bert_ratio <ラベル> <bertotti 行> <変える行> <期待比>
	sed -e "s/^bertotti = .*/$2/" "$SRC/bertotti_core.ofe" > "$WORK/be0.ofe"
	(cd "$WORK" && "$OFE" -n 2 be0.ofe > /dev/null && "$OFE_POST" > /dev/null)
	p0=$(pfe_of)
	sed -e "s/^bertotti = .*/$2/" -e "$3" "$SRC/bertotti_core.ofe" > "$WORK/be1.ofe"
	(cd "$WORK" && "$OFE" -n 2 be1.ofe > /dev/null && "$OFE_POST" > /dev/null)
	p1=$(pfe_of)
	compare "$1" "$(awk -v a="$p1" -v b="$p0" 'BEGIN{ printf "%.10e", a / b }')" "$4" 1e-6
}
# ヒステリシス項だけ (ke = 0, d = 0)
bert_ratio "hysteresis: f x2"  "bertotti = 2 200.0 2.0 0.0 0.0" "s/^frequency = .*/frequency = 100/" 2.0
bert_ratio "hysteresis: I x2"  "bertotti = 2 200.0 2.0 0.0 0.0" "s/^current = .*/current = 2/"        4.0
# 古典渦電流項だけ (kh = 0, ke = 0)
bert_ratio "classical: f x2"   "bertotti = 2 0.0 2.0 0.0 0.35e-3" "s/^frequency = .*/frequency = 100/" 4.0
bert_ratio "classical: I x2"   "bertotti = 2 0.0 2.0 0.0 0.35e-3" "s/^current = .*/current = 2/"       4.0
# 異常 (過剰) 項だけ (kh = 0, d = 0)
bert_ratio "excess: f x2"      "bertotti = 2 0.0 2.0 1.5 0.0" "s/^frequency = .*/frequency = 100/" 2.8284271247
bert_ratio "excess: I x2"      "bertotti = 2 0.0 2.0 1.5 0.0" "s/^current = .*/current = 2/"       2.8284271247
# alpha が実際に効くこと (alpha = 2 -> 3 で I x2 の比が 4 -> 8)
bert_ratio "hysteresis: alpha = 3, I x2" "bertotti = 2 200.0 3.0 0.0 0.0" "s/^current = .*/current = 2/" 8.0

# 積層厚だけ与えて sigma が無いと古典項が黙って 0 になるので警告すること
sed -e 's/^material = 1.0 2e6/material = 1.0 0/' "$SRC/bertotti_core.ofe" > "$WORK/be_nosig.ofe"
if (cd "$WORK" && "$OFE" -n 2 be_nosig.ofe 2>&1 | grep -q "classical eddy-current term is 0"); then
	echo "  bertotti with sigma = 0 warns : OK"
else
	echo "  bertotti with sigma = 0 warns : NG" >&2
	status=1
fi
# frequency 無しは弾くこと
grep -v '^frequency = ' "$SRC/bertotti_core.ofe" > "$WORK/be_nofreq.ofe"
mesh_reject "bertotti without frequency" be_nofreq.ofe

# 3 次元渦電流 (A-φ、辺要素)。1 次元厳密解 Z = γL/(2σW tanh(γt/2)) と比較する。
# ω→0 では R が DC 抵抗 L/(σWt) に厳密一致しなければならない (連成系全体の検査)
echo "[bar_eddy] 3D eddy current (A-phi) vs the 1-D exact skin effect"
cp "$SRC/bar_eddy.ofe" "$WORK/"
for m in "$SRC"/*.msh; do
	[ -f "$m" ] && cp "$m" "$WORK/"
done
# freq  R[ohm]        L[H]           許容 (要素の最大寸法 / 表皮深さ で決まる)
for pair in "1e2 1.37931436e-04 8.37757344e-10 0.005" \
            "1e4 1.41899129e-04 8.30877185e-10 0.005" \
            "1e5 3.24859232e-04 5.34552067e-10 0.02"; do
	set -- $pair
	freq=$1
	sed "s/^frequency = .*/frequency = $freq/" "$SRC/bar_eddy.ofe" > "$WORK/bar_eddy_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 bar_eddy_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R(f=$freq) [ohm]" "$(value_of Rf)" "$2" "$4"
	compare "L(f=$freq) [H]" "$(value_of Lf)" "$3" "$4"
	if grep -q "NOT converged" "$WORK/ofe.log"; then
		echo "  *** A-phi solver did not converge at $freq Hz" >&2
		status=1
	fi
done

# ゲージ固定 (gauge = 1) は Z を変えてはいけない。awall の辺を優先して木に入れる
# 処理が壊れると端子電圧がずれるので、その不変性を直接検査する
# (既定 gauge = 0 のケースだけでは edge_tree_gauge() が一度も実行されない)
sed -e "s/^frequency = .*/frequency = 1e4/" -e "s/^awall = 20/awall = 20\ngauge = 1/" \
    "$SRC/bar_eddy.ofe" > "$WORK/bar_eddy_gauge.ofe"
# 木が awall を壊すと ofe 自身が電極の連結を検出して落とすので、
# その失敗を set -e の打ち切りではなく読める形で報告する
if (cd "$WORK" && "$OFE" -n 2 bar_eddy_gauge.ofe > /dev/null && "$OFE_POST" > /dev/null); then
	compare "R(gauge=1, f=1e4) [ohm]" "$(value_of Rf)" 1.41899129e-04 0.005
	compare "L(gauge=1, f=1e4) [H]" "$(value_of Lf)" 8.30877185e-10 0.005
else
	echo "  gauge=1 : the run failed -> NG" >&2
	sed -n 's/^\*\*\* /    /p' "$WORK/ofe.log" >&2
	status=1
fi

# 磁性導体 (mur != 1)。ν = (μ0 μr)^-1 の材料参照を通す唯一のケース。
# 厳密解は同じ閉形式で γ = sqrt(jω 50 μ0 σ) としたもの
sed -e "s/^frequency = .*/frequency = 1e2/" -e "s/^region = 1 2/region = 1 2\nmur = 2 50/" \
    "$SRC/bar_eddy.ofe" > "$WORK/bar_eddy_mur.ofe"
(cd "$WORK" && "$OFE" -n 2 bar_eddy_mur.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "R(mur=50, f=1e2) [ohm]" "$(value_of Rf)" 1.38932306e-04 0.005
compare "L(mur=50, f=1e2) [H]" "$(value_of Lf)" 4.18010403e-08 0.005
# 磁性導体では表皮深さに μr が効く。診断が μ0 だけで計算していると
# 「刻めている」と誤報するので、警告が出るべきケースで出ることを見る
sed -e "s/^frequency = .*/frequency = 1e5/" -e "s/^region = 1 2/region = 1 2\nmur = 2 50/" \
    "$SRC/bar_eddy.ofe" > "$WORK/bar_eddy_mur_coarse.ofe"
(cd "$WORK" && "$OFE" -n 2 bar_eddy_mur_coarse.ofe > /dev/null) || true
if grep -q "exceeds the skin depth" "$WORK/ofe.log"; then
	echo "  skin-depth warning (mur=50, f=1e5) : raised -> OK"
else
	echo "  skin-depth warning (mur=50, f=1e5) : missing -> NG" >&2
	status=1
fi

# 周波数掃引。**掃引の各点が個別実行とビット単位で一致すること**を見る。
# これは掃引間で状態が漏れる誤りを直接捕まえる (分散材料の再展開漏れで
# 実際に C が最初の周波数で凍結する不具合が出た)。
echo "[sweep_plate] frequency sweep : each point must equal a separate run"
cp "$SRC/sweep_plate.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 sweep_plate.ofe > /dev/null)
# 閉形式との比較 (Lorentz 共鳴をまたぐので C が 3.6 倍変わる)
n=0
for pair in "1.00000000e+08 2.64414642e-13 8.61164622e-06" \
            "3.16227766e+08 2.54612385e-13 7.95141286e-05" \
            "1.00000000e+09 2.12380856e-13 4.77387017e-04" \
            "3.16227766e+09 7.30183021e-14 9.03125045e-04" \
            "1.00000000e+10 8.80172244e-14 8.32245471e-04"; do
	set -- $pair
	n=$((n + 1))
	row=$(awk -F, -v r="$n" 'NR == r + 1 { print }' "$WORK/ofe_sweep.csv")
	compare "sweep f=$1 : f [Hz]" "$(echo "$row" | cut -d, -f1)" "$1" 1e-6
	compare "sweep f=$1 : C [F]"  "$(echo "$row" | cut -d, -f2)" "$2" 1e-6
	compare "sweep f=$1 : G [S]"  "$(echo "$row" | cut -d, -f3)" "$3" 1e-6
done
# 掃引の 1 点を個別に回して**完全一致**すること (状態の漏れの検出)
for f in 1e8 1e10; do
	sed -e "s/^frequencysweep = .*/frequency = $f/" "$SRC/sweep_plate.ofe" > "$WORK/sweep_one.ofe"
	(cd "$WORK" && "$OFE" -n 2 sweep_one.ofe > /dev/null && "$OFE_POST" > /dev/null)
	r=$(awk -F, -v ff="$f" 'NR > 1 { if (($1 + 0) == (ff + 0)) print }' "$WORK/ofe_sweep.csv")
	compare "sweep vs single run (f=$f) C" "$(echo "$r" | cut -d, -f2)" "$(value_of C)" 1e-12
	compare "sweep vs single run (f=$f) G" "$(echo "$r" | cut -d, -f3)" "$(value_of G)" 1e-12
done
# log / lin の生成が正しいこと
sed "s/^frequencysweep = .*/frequencysweep = lin 1e9 5e9 3/" "$SRC/sweep_plate.ofe" > "$WORK/sweep_lin.ofe"
(cd "$WORK" && "$OFE" -n 2 sweep_lin.ofe > /dev/null)
compare "lin sweep point 2 [Hz]" "$(awk -F, 'NR == 3 { print $1 }' "$WORK/ofe_sweep.csv")" 3.0e9 1e-9
# currentsweep との併用は禁止。ja がある正当なケース (plate_line_ja) に
# frequencysweep を足す。ja 抜きで currentsweep を書くと「ja が要る」という
# 別の検査に引っかかり、併用禁止が外れても気づけない
sed "s/^analysis = .*/frequencysweep = 1e3 1e4\nanalysis = M/" \
    "$SRC/plate_line_ja.ofe" > "$WORK/sweep_bad.ofe"
if (cd "$WORK" && "$OFE" -n 2 sweep_bad.ofe > /dev/null 2>&1); then
	echo "  frequencysweep + currentsweep rejected : no -> NG" >&2
	status=1
else
	echo "  frequencysweep + currentsweep rejected : yes -> OK"
fi

# 場の出力 (fieldout = 1)。**書き出した場から集中定数を作り直して**元の抽出値と
# 比べる。厳密な恒等式なので、場が壊れていれば必ず落ちる
#   ∫ ½ε|E|² dV = ½CV²          (静電界)
#   ∫ |J|²/(2σ) dV = ½Re(V I*)  (渦電流。J が一様な低周波では厳密)
echo "[fieldout] the written field must reproduce the extracted lumped values"

# (0) 向き・分母・並べ替えの検査。直列 2 材料で E が区分一様になり、
#     3 軸とも長さも分割数も違い、通電方向は不等間隔にしてある。
#     体積加重平均は**符号つき**なので E = -∇φ の向きの誤りが落ちる。
#     xcut で界面の両側に分けた平均は、セルと座標の対応が崩れると入れ替わる。
for pair in "x 0" "y 1"; do
	set -- $pair
	cp "$SRC/field_probe_$1.ofe" "$WORK/"
	(cd "$WORK" && "$OFE" -n 2 "field_probe_$1.ofe" > /dev/null)
	compare "mean E_$1 (signed) [V/m]" "$(vtk mean$2 E_R_port1 $2 0.4e-3)" -1.0e3 1e-6
	compare "E_$1 below the interface [V/m]" "$(vtk lo E_R_port1 $2 0.4e-3)" -1.8181818182e3 1e-6
	compare "E_$1 above the interface [V/m]" "$(vtk hi E_R_port1 $2 0.4e-3)" -4.5454545455e2 1e-6
	o=$(vtk amax0 E_R_port1); p=$(vtk amax1 E_R_port1); q=$(vtk amax2 E_R_port1)
	res=$(awk -v a="$o" -v b="$p" -v c="$q" -v ax="$2" \
		'BEGIN{ m = (ax == 0) ? (b > c ? b : c) : (a > c ? a : c)
		        printf "%s", (m < 1e-6) ? "OK" : "NG" }')
	echo "  E transverse to $1 ~ 0 : $res"
	case "$res" in NG*) status=1 ;; esac
done

# (1) 構造格子 + 静電界 : E は誘電体内で厳密に一様、C は恒等式で戻る

awk '/^analysis = /{print "fieldout = 1"} {print}' "$SRC/parallel_plate.ofe" > "$WORK/fld_pp.ofe"
(cd "$WORK" && "$OFE" -n 2 fld_pp.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "field vol [m^3]" "$(vtk vol E_C_port1)" 2.0e-10 1e-9
compare "|E| min [V/m]" "$(vtk vmin E_C_port1)" 5.0e3 1e-9
compare "|E| max [V/m]" "$(vtk vmax E_C_port1)" 5.0e3 1e-9
compare "mean E_z (signed) [V/m]" "$(vtk mean2 E_C_port1)" -5.0e3 1e-9
res=$(awk -v a="$(vtk amax0 E_C_port1)" -v b="$(vtk amax1 E_C_port1)" \
	'BEGIN{ printf "%s", ((a < 1e-3) && (b < 1e-3)) ? "OK" : "NG" }')
echo "  E transverse ~ 0 : max|Ex|=$(vtk amax0 E_C_port1) -> $res"
case "$res" in NG*) status=1 ;; esac
cfld=$(awk -v i2="$(vtk int2 E_C_port1)" 'BEGIN{ printf "%.10e", 8.8541878128e-12 * 4 * i2 }')
compare "C from the field [F]" "$cfld" "$(value_of C)" 1e-6

# (2) 非構造格子 + 静電界 : 同じ恒等式が四面体でも成り立つこと
awk '/^analysis = /{print "fieldout = 1"} {print}' "$SRC/box_tet.ofe" > "$WORK/fld_box.ofe"
(cd "$WORK" && "$OFE" -n 2 fld_box.ofe > /dev/null && "$OFE_POST" > /dev/null)
cfld=$(awk -v i2="$(vtk int2 E_C_port1)" 'BEGIN{ printf "%.10e", 8.8541878128e-12 * 4 * i2 }')
compare "C from the field (tet) [F]" "$cfld" "$(value_of C)" 1e-6
# |E|^2 の恒等式は符号に無感なので、四面体側も符号つきで見る (E = -∇φ の向き)
compare "mean E_z (tet, signed) [V/m]" "$(vtk mean2 E_C_port1)" -5.0e3 1e-9

# (2b) 2 次要素 : 同じ恒等式に加えて、セルの型と節点数が VTK の 2 次四面体に
#      なっていること。∇N は重心で評価するので、E が一様な平行平板では厳密
awk '/^analysis = /{print "fieldout = 1"} {print}' "$SRC/box_p2.ofe" > "$WORK/fld_p2.ofe"
(cd "$WORK" && "$OFE" -n 2 fld_p2.ofe > /dev/null && "$OFE_POST" > /dev/null)
cfld=$(awk -v i2="$(vtk int2 E_C_port1)" 'BEGIN{ printf "%.10e", 8.8541878128e-12 * 4 * i2 }')
compare "C from the field (P2) [F]" "$cfld" "$(value_of C)" 1e-6
compare "mean E_z (P2, signed) [V/m]" "$(vtk mean2 E_C_port1)" -5.0e3 1e-9
# VTK_QUADRATIC_TETRA = 24、1 セルあたり 10 節点。あわせて中間節点の並べ替え
# (Gmsh -> VTK) も見る。VTK の規約では中間節点 m はその辺の中点に載るはずで、
# 辺 (1,3) と (2,3) を入れ替えると 2 個ずつ外れる。型だけ見ても並びの誤りは
# 検出できない (ParaView は黙って歪んだ要素を描く)
res=$(awk '
	NF == 0       { next }
	/^POINTS/     { st = "p"; k = 0; next }
	/^CELLS/      { st = "c"; k = 0; next }
	/^CELL_TYPES/ { st = "t"; k = 0; next }
	/^[A-Z_]+ /   { st = ""; next }
	st == "p" { px[k] = $1; py[k] = $2; pz[k] = $3; k++; next }
	st == "t" { k++; if ($1 != 24) type++; next }
	st == "c" {
		if ($1 != 10) { nn++; next }
		# VTK の並び : 頂点 0..3、中間 (0,1)(1,2)(0,2)(0,3)(1,3)(2,3)
		split("0 1 1 2 0 2 0 3 1 3 2 3", ev, " ")
		for (m = 0; m < 6; m++) {
			a = $(2 + ev[(2 * m) + 1]); b = $(2 + ev[(2 * m) + 2]); c = $(6 + m)
			d = ((px[c] - ((px[a] + px[b]) / 2)) ^ 2) \
			  + ((py[c] - ((py[a] + py[b]) / 2)) ^ 2) \
			  + ((pz[c] - ((pz[a] + pz[b]) / 2)) ^ 2)
			if (d > 1e-24) off++
		}
	}
	END { if (k == 0) { printf "NG (no cells)" }
	      else if (type || nn) { printf "NG (%d cells not type 24, %d not 10-node)", type, nn }
	      else if (off) { printf "NG (%d mid-nodes off their edge)", off }
	      else { printf "OK" } }' "$WORK/ofe_field.vtk")
echo "  VTK quadratic tetra (type, node count, mid-node order) : $res"
case "$res" in NG*) status=1 ;; esac

# (3) 3 次元渦電流 : 場の損失と端子電力が一致すること。
#     J は重心で評価するので、J が一様な低周波でのみ厳密になる (100 Hz で 8 桁)
sed -e "s/^frequency = .*/frequency = 1e2/" \
    -e "s/^analysis = A/fieldout = 1\nanalysis = A/" "$SRC/bar_eddy.ofe" > "$WORK/fld_bar.ofe"
(cd "$WORK" && "$OFE" -n 2 fld_bar.ofe > /dev/null)
compare "field vol (tet) [m^3]" "$(vtk vol J_A_re)" 5.0e-10 1e-9
# DC の電流密度は σ V / L = 5.8e7 / 2e-3 = 2.9e10 [A/m^2] で一様
compare "|J| max [A/m^2]" "$(vtk vmax J_A_re)" 2.9e10 1e-4
compare "|J| min [A/m^2]" "$(vtk vmin J_A_re)" 2.9e10 1e-4
compare "mean J_x (signed) [A/m^2]" "$(vtk mean0 J_A_re)" -2.9e10 1e-4
pf=$(awk -v r="$(vtk int2 J_A_re)" -v m="$(vtk int2 J_A_im)" \
	'BEGIN{ printf "%.10e", (r + m) / (2 * 5.8e7) }')
pt=$(awk '/port 1 : I =/{ printf "%.10e", 0.5 * $6 }' "$WORK/ofe.log")
compare "ohmic loss from the field [W]" "$pf" "$pt" 1e-6

# (4) 2 次元渦電流 (F) の場。オーム損が端子から見た ½Re(Y) と一致すること。
#     F は V' = 1 [V/m] 励振なので Y = 1/Z より ½Re(Y) = ½R/(R²+X²)
awk '/^analysis = /{print "fieldout = 1"} {print}' "$SRC/plate_line_ac.ofe" > "$WORK/fld_ac.ofe"
(cd "$WORK" && "$OFE" -n 2 fld_ac.ofe > /dev/null && "$OFE_POST" > /dev/null)
tl=$(awk '/Transmission line axis/ { for (i = 1; i <= NF; i++) if ($i == "=") l = $(i+1); print l }' \
	"$WORK/ofe.log" | tail -1)
pf=$(awk -v a="$(vtk int2 J_F_re_port1)" -v b="$(vtk int2 J_F_im_port1)" -v t="$tl" \
	'BEGIN{ printf "%.10e", (a + b) / (2 * 5.8e7) / t }')
pt=$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
	        printf "%.10e", 0.5 * R / ((R * R) + (X * X)) }')
compare "ohmic loss from the field (F) [W/m]" "$pf" "$pt" 1e-6
# 磁気エネルギーの恒等式 : ω∫|B|²/μ0 dA = X/(R²+X²)。
# **B の実部と虚部の両方**が効くので、片方を取り違えると落ちる
qf=$(awk -v a="$(vtk int2 B_F_re_port1)" -v b="$(vtk int2 B_F_im_port1)" -v t="$tl" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3
	        printf "%.10e", om * (a + b) / (4e-7 * 3.14159265358979324) / t }')
qt=$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
	        printf "%.10e", X / ((R * R) + (X * X)) }')
compare "magnetic energy from the field (F)" "$qf" "$qt" 1e-6
# 導体毎の電流。|J|² は符号に無感なので、**符号つきで**端子アドミタンスと比べる。
# 2 枚の板は y で分かれ、電流は tline 軸 (z) 向きなので分割軸と成分を分ける
compare "Re(I) of the driven conductor" \
	"$(awk -v v="$(vtk ihi J_F_re_port1 1 0.15e-3 2)" -v t="$tl" 'BEGIN{ printf "%.10e", v / t }')" \
	"$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
		'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
		        printf "%.10e", R / ((R * R) + (X * X)) }')" 1e-6
compare "Im(I) of the driven conductor" \
	"$(awk -v v="$(vtk ihi J_F_im_port1 1 0.15e-3 2)" -v t="$tl" 'BEGIN{ printf "%.10e", v / t }')" \
	"$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
		'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
		        printf "%.10e", -X / ((R * R) + (X * X)) }')" 1e-6
# B の向き。板は x 方向に一様なので ∂Az/∂x = 0、B = ∇×(Az ẑ) は **x 向き**になる。
# 回転を取らず -∇Az にすると同じ大きさで y 向きになるため、|B|² の恒等式では
# 見えず、この成分比較でだけ落ちる
res=$(awk -v bx="$(vtk amax0 B_F_re_port1)" -v by="$(vtk amax1 B_F_re_port1)" \
	'BEGIN{ printf "%s", ((bx > 0) && (by < 1e-6 * bx)) ? "OK" : "NG" }')
echo "  B_F is x-directed (curl, not gradient) : $res"
case "$res" in NG*) status=1 ;; esac
# 帰路の導体はちょうど符号が逆 (伝送線路の前提)
compare "Re(I) of the return conductor" \
	"$(awk -v v="$(vtk ilo J_F_re_port1 1 0.15e-3 2)" 'BEGIN{ printf "%.10e", -v }')" \
	"$(vtk ihi J_F_re_port1 1 0.15e-3 2)" 1e-9

# (4b) 同じ恒等式を**断面 2 次元の非構造格子**でも通す。三角形セルなので
#      vtkcheck.awk はセル面積を使う (単位長あたりなので面積 = 「体積」)。
#      この経路は場の出力を足したときに検証が無いまま出荷していた
awk '/^analysis = F/{print "fieldout = 1"} {print}' "$SRC/plate2d_ac.ofe" > "$WORK/fld2d.ofe"
(cd "$WORK" && "$OFE" -n 2 fld2d.ofe > /dev/null && "$OFE_POST" > /dev/null)
# 断面積 = W (2t + d) = 1e-3 * 0.3e-3
compare "field area (2-D mesh) [m^2]" "$(vtk vol J_F_re_port1)" 3.0e-7 1e-9
# VTK_TRIANGLE = 5、1 セルあたり 3 節点。**vtkcheck.awk は先頭の節点数で
# 判定するのでセル型の誤りは見えない**ので、型そのものを別に assert する
# (型が違うと ParaView が開けない)
res=$(awk '
	NF == 0       { next }
	/^CELLS/      { st = "c"; k = 0; next }
	/^CELL_TYPES/ { st = "t"; k = 0; next }
	/^[A-Z_]+ /   { st = ""; next }
	st == "t" { k++; if ($1 != 5) type++; next }
	st == "c" { if ($1 != 3) nn++; next }
	END { if (k == 0) { printf "NG (no cells)" }
	      else if (type || nn) { printf "NG (%d not type 5, %d not 3-node)", type, nn }
	      else { printf "OK" } }' "$WORK/ofe_field.vtk")
echo "  VTK cell type is TRIANGLE with 3 nodes : $res"
case "$res" in NG*) status=1 ;; esac
pf=$(awk -v a="$(vtk int2 J_F_re_port1)" -v b="$(vtk int2 J_F_im_port1)" \
	'BEGIN{ printf "%.10e", (a + b) / (2 * 5.8e7) }')
pt=$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
	        printf "%.10e", 0.5 * R / ((R * R) + (X * X)) }')
compare "ohmic loss from the field (2-D) [W/m]" "$pf" "$pt" 1e-4
qf=$(awk -v a="$(vtk int2 B_F_re_port1)" -v b="$(vtk int2 B_F_im_port1)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3
	        printf "%.10e", om * (a + b) / (4e-7 * 3.14159265358979324) }')
qt=$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
	        printf "%.10e", X / ((R * R) + (X * X)) }')
compare "magnetic energy from the field (2-D)" "$qf" "$qt" 1e-3
# B の向き : 板は y (幅方向) に一様なので ∂Az/∂y = 0、B = ∇×(Az x̂) は **y 向き**。
# 回転を取らず -∇Az にすると同じ大きさで z 向きになるので、大きさの恒等式では
# 見えずこの比較でだけ落ちる。
#
# **最大値ではなく符号つきの体積加重平均で見ること。** 四角形を対角線で 2 つの
# 三角形に割っているので、z にだけ変化する場でも重心での ∂Az/∂y が対角線の
# 向きに応じて交番し、max|B_z| は max|B_y| の 0.9% ほど残る (実測 1.6e-5 対
# 1.8e-3)。符号つき平均ではこれが打ち消し合って 1e-11 まで落ちる
res=$(awk -v by="$(vtk mean1 B_F_re_port1)" -v bz="$(vtk mean2 B_F_re_port1)" \
	'BEGIN{ if (by < 0) by = -by; if (bz < 0) bz = -bz
	        printf "%s", ((by > 0) && (bz < 1e-6 * by)) ? "OK" : "NG" }')
echo "  B_F is y-directed on the 2-D mesh (curl, not gradient) : $res"
case "$res" in NG*) status=1 ;; esac
# 符号つきの導体電流。2 枚の板は z で分かれ、電流は tline 軸 (x) 向きなので
# 分割軸 (z=2) と積分する成分 (x=0) を分ける。界面は z = 0.1e-3 (間隙の中央)
compare "Re(I) of the return conductor (2-D)" \
	"$(awk -v v="$(vtk ilo J_F_re_port1 2 0.1e-3 0)" 'BEGIN{ printf "%.10e", -v }')" \
	"$(vtk ihi J_F_re_port1 2 0.1e-3 0)" 1e-6


# 非導電領域 (空気) を含む 3 次元渦電流。空気層は界面で Robin 条件に潰れるので
# 1 次元の閉形式が残る。**空気が効いていることを見るには L を見る必要がある**
# (空気の辺を A=0 に固定する変異は R をほとんど動かさず L だけ -43% ずらす)
echo "[bar_air] A-phi with a non-conducting layer : 1-D exact with an air gap"
cp "$SRC/bar_air.ofe" "$WORK/"
for m in "$SRC"/*.msh; do
	[ -f "$m" ] && cp "$m" "$WORK/"
done
for pair in "1e2 1.3793210584e-04 1.6755132157e-09 0.005" \
            "1e4 1.4831436443e-04 1.6479186641e-09 0.005" \
            "1e5 4.7930357804e-04 9.7419014024e-10 0.02"; do
	set -- $pair
	sed "s/^frequency = .*/frequency = $1/" "$SRC/bar_air.ofe" > "$WORK/bar_air_run.ofe"
	(cd "$WORK" && "$OFE" -n 2 bar_air_run.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R(f=$1) [ohm]" "$(value_of Rf)" "$2" "$4"
	compare "L(f=$1) [H]" "$(value_of Lf)" "$3" "$4"
	if grep -q "NOT converged" "$WORK/ofe.log"; then
		echo "  *** A-phi solver did not converge at $1 Hz (air)" >&2
		status=1
	fi
done
# 磁性導体 + 空気。**1 つの格子に 2 つの異なる ν が同居する唯一のケース**なので、
# 要素毎の材料参照を取り違える誤りを捕まえられるのはここだけ
sed -e "s/^frequency = .*/frequency = 1e2/" -e "s/^region = 1 2/region = 1 2\nmur = 2 50/" \
    "$SRC/bar_air.ofe" > "$WORK/bar_air_mur.ofe"
(cd "$WORK" && "$OFE" -n 2 bar_air_mur.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "R(mur=50, f=1e2) [ohm]" "$(value_of Rf)" 1.3893378002e-04 0.005
compare "L(mur=50, f=1e2) [H]" "$(value_of Lf)" 4.3045141586e-08 0.005
# 空気があってもゲージ不変性は保たれること
sed -e "s/^frequency = .*/frequency = 1e4/" -e "s/^awall = 20/awall = 20\ngauge = 1/" \
    "$SRC/bar_air.ofe" > "$WORK/bar_air_gauge.ofe"
if (cd "$WORK" && "$OFE" -n 2 bar_air_gauge.ofe > /dev/null && "$OFE_POST" > /dev/null); then
	compare "R(gauge=1, air) [ohm]" "$(value_of Rf)" 1.4831436443e-04 0.005
	compare "L(gauge=1, air) [H]" "$(value_of Lf)" 1.6479186641e-09 0.005
else
	echo "  gauge=1 (air) : the run failed -> NG" >&2
	status=1
fi

# 未使用キーの警告。σ の読み出しは Material[].sigma (R/A/E) と CondSigma[] (Rs/F) の
# 2 系統に分かれていて、取り違えても値が 0 になるだけで黙って通る。
# 警告が (a) 罠のケースで必ず出て (b) 正しいケースでは 1 件も出ないことを両方見る
echo "[input lint] keys the selected analysis cannot read"
lint_expect() {	# lint_expect <label> <ofe> <yes|no>
	nw=$(cd "$WORK" && "$OFE" -n 2 "$2" 2>&1 | grep -c "warning" || true)
	if [ "$3" = yes ]; then
		if [ "$nw" -gt 0 ]; then echo "  $1 : warned -> OK"
		else echo "  $1 : no warning -> NG" >&2; status=1; fi
	else
		if [ "$nw" -eq 0 ]; then echo "  $1 : silent -> OK"
		else echo "  $1 : $nw spurious warning(s) -> NG" >&2; status=1; fi
	fi
}
cp "$SRC"/*.ofe "$WORK/" 2>/dev/null || true
for m in "$SRC"/*.msh; do [ -f "$m" ] && cp "$m" "$WORK/"; done
awk '/^analysis = /{print "tempco = 0 3.93e-3 20"} {print}' \
    "$SRC/plate_line_ac.ofe" > "$WORK/lint_tempco_f.ofe"
awk '/^analysis = /{print "conductortempco = 1 3.93e-3 20"} {print}' \
    "$SRC/bar_eddy.ofe" > "$WORK/lint_ctempco_a.ofe"
awk '/^analysis = /{print "mur = 2 100"} {print}' \
    "$SRC/parallel_plate.ofe" > "$WORK/lint_mur_c.ofe"
sed "s/^material = .*/material = 4.0 1e3/" \
    "$SRC/parallel_plate.ofe" > "$WORK/lint_sigma_c.ofe"
sed "s/^analysis = .*/analysis = L/" "$SRC/coax.ofe" > "$WORK/lint_l_only.ofe"
lint_expect "tempco + analysis=F"          lint_tempco_f.ofe  yes
lint_expect "conductortempco + analysis=A" lint_ctempco_a.ofe yes
lint_expect "mur + analysis=C"             lint_mur_c.ofe     yes
lint_expect "material sigma + analysis=C"  lint_sigma_c.ofe   yes
lint_expect "analysis=L alone"             lint_l_only.ofe    yes
# 正しいケースでは 1 件も出ないこと (誤検知は警告を無視させるので同じくらい悪い)
for c in parallel_plate resistor_bar coax microstrip plate_line_dc coax_loss \
         plate_line_ac plate_line_bh dispersive_plate drude_plate colecole_plate \
         temp_resistor aniso_plate box_tet coax_tet edge_test bar_eddy \
         box_p2 coax_p2 nodal_test_p2 plate2d_dc plate2d_ac plate2d_r30 \
         bertotti_core temp_material temp_mur hn_plate plate2d_p2 \
         box_tet_41 plate2d_41 coupled_microstrip; do
	lint_expect "$c (clean)" "$c.ofe" no
done

# ヒステリシス (Jiles-Atherton) : H = I/W が Ampere の法則で厳密に決まるので、
# FEM の結果はスカラー J-A モデルを H 掃引で積分した ODE 解と一致しなければならない
# (期待値は data/sample/plate_line_ja.ofe のコメント参照)
echo "[plate_line_ja] hysteresis loop vs independent ODE integration"
run_case plate_line_ja
n=0
for bexp in 0.355138 0.912045 1.342832 1.260954 0.713845 \
            -0.908562 -1.342754 -1.260921 -0.713844 0.908562 1.342754; do
	n=$((n + 1))
	bval=$(awk -v s="$n" '/Jiles-Atherton hysteresis/ { inb = 1; next }
	                       inb && ($1 == s) && (NF == 6) { print $4; exit }' "$WORK/ofe.log")
	compare "B(step $n) [T]" "$bval" "$bexp" 0.005
done
if grep -q "did not converge" "$WORK/ofe.log"; then
	echo "  *** hysteresis iteration did not converge" >&2
	status=1
fi

# 系列の HDF5 出力 (hdf5 = 1)。**HDF5 は任意依存**なので、
#   ・入っていないビルド : 「hdf5 = 1 が弾かれること」だけを見て残りは skip する
#     (黙って通すと、書けないのに気づかないまま緑になる)
#   ・入っているビルド   : 既存の出口 (ofe_sweep.csv / ofe.log / ofe_field.vtk) と
#     **同じ数字が入っていること**を恒等式にする。系列専用の閉形式は要らないし、
#     「2 つの出口が食い違う」という一番起きやすい壊れ方をそのまま検出できる
echo "[hdf5] the time-series output (ofe_series.h5)"
sed 's/^analysis = /hdf5 = 1\nanalysis = /' "$SRC/parallel_plate.ofe" > "$WORK/h5probe.ofe"
h5msg=$(cd "$WORK" && "$OFE" -n 2 h5probe.ofe 2>&1 || true)
case "$h5msg" in
*"needs a build with HDF5"*)
	echo "  built without HDF5 : hdf5 = 1 is rejected -> OK (the rest is skipped)"
	;;
*)
	if ! command -v h5dump > /dev/null 2>&1; then
		echo "  *** the binary supports HDF5 but h5dump is missing; cannot verify" >&2
		status=1
	else

	# データセットを 1 行 1 値で取り出す。h5v <名前> [start] [count]
	#
	# **桁数は比べる相手に合わせること** ($H5FMT)。csv は %.8e、VTK は %.9e、
	# ログは %13.6e で書かれているので、同じ倍精度値を同じ桁で丸めれば
	# 文字列として一致する。取り出しと比較で桁が違うと二重丸めで
	# 最後の桁がずれ、正しい実装が落ちる (実際に踏んだ)
	H5FMT="%.8e"
	h5v() {
		if [ -n "$2" ]; then
			h5dump -d "$1" -s "$2" -c "$3" -m "$H5FMT" -y -w 1 -A 0 "$WORK/ofe_series.h5"
		else
			h5dump -d "$1" -m "$H5FMT" -y -w 1 -A 0 "$WORK/ofe_series.h5"
		fi 2>/dev/null | sed -n '/DATA {/,/}/p' | tr -d ' ,' | grep -E '^-?[0-9]' || true
		# データセットが無いときは空を返す。ここで落とすと set -e で
		# スクリプトごと死に、どの検査が NG なのか分からなくなる
	}
	# 2 つの数値列が (指定した桁で) 一致すること
	same_nums() {	# same_nums <ラベル> <file a> <file b> <printf 書式>
		# 空のファイルを awk の NR == FNR に食わせると 2 つ目が 1 つ目として
		# 読まれ、件数の報告が入れ替わる。先に弾いて診断を正直にする
		if [ ! -s "$2" ]; then
			echo "  $1 : NG (no values in the HDF5 file)" >&2
			status=1
			return
		fi
		if [ ! -s "$3" ]; then
			echo "  $1 : NG (no reference values)" >&2
			status=1
			return
		fi
		res=$(awk -v f="${4:-%.8e}" '
			NR == FNR { a[FNR] = sprintf(f, $1 + 0); na = FNR; next }
			{ b[FNR] = sprintf(f, $1 + 0); nb = FNR }
			END { if ((na == 0) || (na != nb)) { printf "NG (%d vs %d values)", na, nb; exit }
			      for (i = 1; i <= na; i++) if (a[i] != b[i]) {
			          printf "NG (row %d : %s vs %s)", i, a[i], b[i]; exit }
			      printf "OK (%d values)", na }' "$2" "$3")
		echo "  $1 : $res"
		case "$res" in NG*) status=1 ;; esac
	}

	# (a) 周波数掃引 : /sweep が ofe_sweep.csv と完全に一致すること
	sed 's/^analysis = /hdf5 = 1\nanalysis = /' "$SRC/sweep_plate.ofe" > "$WORK/h5sw.ofe"
	(cd "$WORK" && "$OFE" -n 2 h5sw.ofe > /dev/null)
	h5v /sweep/frequency > "$WORK/h5f.txt"
	awk -F, 'NR > 1 { print $1 }' "$WORK/ofe_sweep.csv" > "$WORK/csvf.txt"
	same_nums "sweep frequency == ofe_sweep.csv" "$WORK/h5f.txt" "$WORK/csvf.txt"
	# 行列は 1 ポートなので (点数 x 1 x 1)。csv の列と 1 対 1 に並ぶ
	q=2
	for m in C G R; do
		h5v "/sweep/$m" > "$WORK/h5m.txt"
		awk -F, -v q="$q" 'NR > 1 { print $q }' "$WORK/ofe_sweep.csv" > "$WORK/csvm.txt"
		same_nums "sweep $m == ofe_sweep.csv" "$WORK/h5m.txt" "$WORK/csvm.txt"
		q=$((q + 1))
	done

	# (b) ヒステリシス履歴 : /hysteresis が ofe.log の表と一致すること。
	# **順序そのものが物理** (B は履歴依存) なので、並べ替えたら別の答えになる
	sed 's/^analysis = /hdf5 = 1\nanalysis = /' "$SRC/plate_line_ja.ofe" > "$WORK/h5ja.ofe"
	(cd "$WORK" && "$OFE" -n 2 h5ja.ofe > /dev/null)
	H5FMT="%.6e"			# ログの表は %13.6e
	awk '$1 ~ /^[0-9]+$/ && NF == 6 { print $3 }' "$WORK/ofe.log" > "$WORK/logh.txt"
	awk '$1 ~ /^[0-9]+$/ && NF == 6 { print $4 }' "$WORK/ofe.log" > "$WORK/logb.txt"
	awk '$1 ~ /^[0-9]+$/ && NF == 6 { print $2 }' "$WORK/ofe.log" > "$WORK/logi.txt"
	h5v /hysteresis/H > "$WORK/h5h.txt"
	h5v /hysteresis/B > "$WORK/h5b.txt"
	h5v /hysteresis/current > "$WORK/h5i.txt"
	# ログは %13.6e なので 6 桁で比べる (HDF5 側は倍精度そのもの)
	same_nums "hysteresis H == ofe.log" "$WORK/h5h.txt" "$WORK/logh.txt" "%.6e"
	same_nums "hysteresis B == ofe.log" "$WORK/h5b.txt" "$WORK/logb.txt" "%.6e"
	same_nums "hysteresis current == ofe.log" "$WORK/h5i.txt" "$WORK/logi.txt" "%.6e"
	# 履歴の向きが保たれていること (B は行って戻る。並べ替えると必ず壊れる)
	res=$(awk 'NR <= 3 { up = ($1 > prev) } { prev = $1 }
		END { printf "%s", (up ? "OK" : "NG") }' "$WORK/h5b.txt")
	echo "  hysteresis order is the magnetisation history : $res"
	case "$res" in NG*) status=1 ;; esac

	# (c) 場のスナップショット。**表皮効果のある F の掃引**を使う:
	# 一様誘電体の C 掃引では φ が周波数に依らない (ε が一様なら ∇・(ε∇φ)=0 の
	# 解は ε に依らない) ので、「点ごとに場が違う」ことを検査できない
	sed -e 's/^analysis = /hdf5 = 1\nfieldout = 1\nanalysis = /' \
	    -e 's/^frequency = .*/frequencysweep = 1e3 1e5 1e7/' \
	    "$SRC/plate_line_ac.ofe" > "$WORK/h5ac.ofe"
	(cd "$WORK" && "$OFE" -n 2 h5ac.ofe > /dev/null)
	H5FMT="%.8e"
	h5v /field/axis > "$WORK/h5ax.txt"
	printf '1e3\n1e5\n1e7\n' > "$WORK/wantax.txt"
	same_nums "field axis == the swept frequencies" "$WORK/h5ax.txt" "$WORK/wantax.txt"
	# 最後のスナップショットが ofe_field.vtk (= 最後の点) と一致すること。
	# 並べ替え (構造格子は x が最内) と点の取り違えはここで落ちる
	awk '/^SCALARS Az_F_re_port1/ { f = 1; getline; next }
	     f && /^[A-Z]/ { f = 0 } f && (NF == 1) { print }' "$WORK/ofe_field.vtk" > "$WORK/vtkaz.txt"
	nn=$(wc -l < "$WORK/vtkaz.txt")
	H5FMT="%.9e"			# ofe_field.vtk は %.9e
	h5v /field/Az_F_re_port1 "2,0" "1,$nn" > "$WORK/h5az2.txt"
	same_nums "last field snapshot == ofe_field.vtk" "$WORK/h5az2.txt" "$WORK/vtkaz.txt" "%.9e"
	# **点ごとに中身が違うこと。** 全点に同じ配列を書く壊れ方は、点数も
	# 最後の点との一致も通ってしまうのでこれだけが捕まえる
	h5v /field/Az_F_re_port1 "0,0" "1,$nn" > "$WORK/h5az0.txt"
	res=$(paste "$WORK/h5az0.txt" "$WORK/h5az2.txt" | awk '
		{ d = $1 - $2; if (d < 0) d = -d; if (d > m) m = d
		  a = ($1 < 0) ? -$1 : $1; if (a > x) x = a }
		END { if (x <= 0) { printf "NG (the field is identically zero)"; exit }
		      printf "%s (max|d| / max|v| = %.2f)", ((m > 0.1 * x) ? "OK" : "NG"), m / x }')
	echo "  the snapshots differ from point to point (skin effect) : $res"
	case "$res" in NG*) status=1 ;; esac

	# (d) 要素種別が混在する格子の /mesh。均一な格子では cells が
	# [ncell][nen] の矩形だが、混在では要素ごとに節点数が違うので平坦化して
	# cell_offset[ncell+1] を添える。**この分岐は混在格子のファイルが無いと
	# 一度も実行されない** (形式の分岐は通るファイルが無いと死んだコードになる)。
	#
	# 閉形式は要らない。同じ連結を ofe_field.vtk にも書いているので、
	# **2 つの出口が食い違わないこと**を恒等式にすれば足りる (一番起きやすい
	# 壊れ方がまさにそれ)。VTK の CELLS は行頭に節点数が付くので、それを
	# 落とした並びが /mesh/cells に、累積和が /mesh/cell_offset に一致する。
	sed 's/^analysis = /hdf5 = 1\nfieldout = 1\nanalysis = /' \
	    "$SRC/box_mixed.ofe" > "$WORK/h5mix.ofe"
	(cd "$WORK" && "$OFE" -n 2 h5mix.ofe > /dev/null)
	awk '/^CELLS /{f = 1; next} f && /^CELL_TYPES/{exit}
	     f && NF { for (i = 2; i <= NF; i++) print $i }' \
	    "$WORK/ofe_field.vtk" > "$WORK/vtkcells.txt"
	awk '/^CELLS /{f = 1; s = 0; print s; next} f && /^CELL_TYPES/{exit}
	     f && NF { s += $1; print s }' "$WORK/ofe_field.vtk" > "$WORK/vtkoff.txt"
	H5FMT="%d"
	h5v /mesh/cells > "$WORK/h5cells.txt"
	h5v /mesh/cell_offset > "$WORK/h5off.txt"
	same_nums "mixed mesh: /mesh/cells == ofe_field.vtk" \
		"$WORK/h5cells.txt" "$WORK/vtkcells.txt" "%d"
	same_nums "mixed mesh: /mesh/cell_offset == ofe_field.vtk" \
		"$WORK/h5off.txt" "$WORK/vtkoff.txt" "%d"
	# **節点数が要素ごとに違うことそのもの**を見る。ここを見ないと、
	# 混在なのに全要素を同じ節点数として書く誤りが上の 2 つを通ってしまう
	# (offset も cells も自分自身と突き合わせているだけになるため)
	res=$(awk 'NR > 1 { d = $1 - p; if (d != f && NR > 2) mixed = 1; f = d }
	           { p = $1 } END { printf "%s (%d cells)", (mixed ? "OK" : "NG"), NR - 1 }' \
		"$WORK/h5off.txt")
	echo "  mixed mesh: the cells really have different node counts : $res"
	case "$res" in NG*) status=1 ;; esac

	fi
	;;
esac

if [ "$status" -ne 0 ]; then
	echo "*** RLC validation FAILED" >&2
else
	echo "RLC validation passed"
fi
exit $status
