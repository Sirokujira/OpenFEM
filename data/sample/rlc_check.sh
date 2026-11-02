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
# conductortempcotable (CondSigma 系統の折れ線)。Rs は抵抗率に厳密に比例するので、
# 直線の表は conductortempco と一致し、折れ線の表では Rs(85) = f(85) Rs(20)
rs_lin=$(value_of Rs)
sed -e "s/^analysis = .*/conductortempcotable = 0 20 1 220 1.786\nconductortempcotable = 1 20 1 220 1.786\ntemperature = 85\n&/" \
    "$SRC/plate_line_dc.ofe" > "$WORK/temp_rs_tab.ofe"
(cd "$WORK" && "$OFE" -n 2 temp_rs_tab.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "conductortempcotable (a straight line) == conductortempco, Rs(85) [ohm/m]" "$(value_of Rs)" "$rs_lin" 1e-12
sed -e "s/^analysis = .*/conductortempcotable = 0 20 1 50 1.2 100 2.0\nconductortempcotable = 1 20 1 50 1.2 100 2.0\ntemperature = 85\n&/" \
    "$SRC/plate_line_dc.ofe" > "$WORK/temp_rs_tab2.ofe"
(cd "$WORK" && "$OFE" -n 2 temp_rs_tab2.ofe > /dev/null && "$OFE_POST" > /dev/null)
rs_tab=$(value_of Rs)
(cd "$WORK" && cp "$SRC/plate_line_dc.ofe" . && "$OFE" -n 2 plate_line_dc.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "conductortempcotable : Rs(85) = 1.76 Rs(20) [ohm/m]" "$rs_tab" \
	"$(awk -v r="$(value_of Rs)" 'BEGIN { printf "%.12e", 1.76 * r }')" 1e-7
sed -e "s/^analysis = .*/conductortempco = 0 3.93e-3 20\nconductortempcotable = 0 20 1 100 1.3\n&/" \
    "$SRC/plate_line_dc.ofe" > "$WORK/temp_rs_both.ofe"
if (cd "$WORK" && "$OFE" -n 2 temp_rs_both.ofe > "$WORK/temp_rs_both.txt" 2>&1); then
	echo "  conductortempco and conductortempcotable together : accepted -> NG" >&2; status=1
elif ! grep -q "both conductortempco and conductortempcotable" "$WORK/temp_rs_both.txt"; then
	echo "  conductortempco and conductortempcotable together : rejected for a different reason -> NG" >&2; status=1
else
	echo "  conductortempco and conductortempcotable together : rejected -> OK"
fi

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
# 辺要素 (E / A) は四面体ならどの次数でも解ける ([A order 2] / [A order 3+])。高次の六面体・角柱は弾く (h3x_e)
# 次数の混在 (四面体を 1 個だけ 1 次に落とす)
awk '/^[0-9]+ 11 2 / && !done { print $1, 4, 2, $4, $5, $6, $7, $8, $9; done = 1; next } { print }' \
    "$SRC/box_p2.msh" > "$WORK/mixed.msh"
sed 's/^mesh = .*/mesh = mixed.msh/' "$SRC/box_p2.ofe" > "$WORK/p2_mixed.ofe"
#   次数の違う要素は最高次数に上げて解く ([mixed orders])。直線の要素なので上げると一様な 2 次と同じ離散問題になる
mx_same() {	# mx_same <ラベル> <元の入力> <次数を混ぜた入力> : rlc.csv が一致し、ログに Mixed orders の行がある
	rm -f "$WORK/rlc.csv"
	(cd "$WORK" && "$OFE" -n 2 "$2" > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1)
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/mxs_a.csv" 2> /dev/null
	rm -f "$WORK/rlc.csv"
	(cd "$WORK" && "$OFE" -n 2 "$3" > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1)
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/mxs_b.csv" 2> /dev/null
	if [ -s "$WORK/mxs_a.csv" ] && cmp -s "$WORK/mxs_a.csv" "$WORK/mxs_b.csv" && grep -q "^Mixed orders : " "$WORK/ofe.log"; then
		echo "  $1 : raised to the highest order, rlc.csv identical to the uniform mesh -> OK"
	else echo "  $1 : not raised, or differs from the uniform mesh -> NG" >&2; status=1; fi
}
cp "$SRC/box_p2.ofe" "$WORK/"
mx_same "mixed element orders (one 4-node tetrahedron among 10-node ones)" box_p2.ofe p2_mixed.ofe
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
direct_same plate2d_p3
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
# **Sloan と RCM の小さい方を使っていること。** 並べ替えは性能のためだけなので
# 値の一致では何も検査できない。使ったプロファイル (分解の報告) が両者の最小値に
# 一致することを、Sloan が勝つ格子 (構造格子 3 次元) と RCM が勝つ格子 (2 次四面体)
# の両方で見る。性能そのものも assert する: 構造格子では RCM の 0.85 倍未満
# (実測 0.80。Sloan の優先度の更新を 1 つ落とすと 0.89、距離の符号を逆にすると
# 0.90 で、どちらも答えは正しいまま)、四面体 (box_tet) ではプロファイルが 19611 以下
# (並べ替えは浮動小数点を使わない整数の計算なので OS に依らず決まる。擬似周辺節点の
# 探索を外すと 20011。**構造格子では次数最小の節点 (角) が最初から周辺なので、
# この誤りは非構造格子でしか出ない**。並べ替えを改良して小さくなるのは構わない)
ord_check() {	# ord_check <ケース> <期待する勝者> [<上限の比>] [<プロファイルの上限>]
	sed 's/^analysis = /direct = 1\nanalysis = /' "$SRC/$1.ofe" > "$WORK/dord.ofe"
	if ! (cd "$WORK" && "$OFE" -n 2 dord.ofe > /dev/null 2>&1); then
		echo "  ordering on $1 : the run failed -> NG" >&2; status=1; return
	fi
	res=$(awk -v want="$2" -v lim="${3:-0}" -v cap="${4:-0}" '
		/ordering : RCM profile/ && (w == "") { for (i = 1; i <= NF; i++) { if (($i == "RCM") && (r == "")) r = $(i+3); if (($i == "Sloan") && (s == "")) s = $(i+3) }
		                                        w = $NF }
		/direct : profile =/ && (u == "") { for (i = 1; i <= NF; i++) if ($i == "profile") u = $(i+2) }
		END { if ((r == "") || (u == "")) { printf "NG (no report)"; exit }
		      m = ((s + 0) < (r + 0)) ? s : r
		      ok1 = ((u + 0) == (m + 0)) && (w == want)
		      ok2 = ((lim == 0) || ((u + 0) < lim * (r + 0))) && ((cap == 0) || ((u + 0) <= (cap + 0)))
		      ok = ok1 && ok2
		      printf "%s (RCM %d, Sloan %d, used %d by %s)", ok ? "OK" : "NG", r, s, u, w }' "$WORK/ofe.log")
	echo "  ordering on $1 : the smaller of RCM / Sloan : $res"
	case "$res" in NG*) status=1 ;; esac
}
ord_check parallel_plate Sloan 0.85
for m in "$SRC"/*.msh; do [ -f "$m" ] && cp "$m" "$WORK/"; done
ord_check box_tet Sloan 0 19611
ord_check box_p2 RCM

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

# ピラミッド (5 節点、四角形の底面 + 頂点)。六面体と四面体をつなぐ遷移要素で、
# 四角形面のトレースが双 1 次 (六面体と適合)、三角形面が 1 次 (四面体と適合)。
#
#   (a) 純ピラミッド格子 : 底面がすべて長方形 (平行四辺形) なので剛性が厳密。
#       閉形式と全桁一致
#   (b) 六面体 -> ピラミッド -> 四面体の遷移格子 : 直列容量の閉形式と全桁一致。
#       種別ごとに物理タグを分けてあるので、要素番号 -> 材料の対応が
#       種別の境目でずれると必ず落ちる
#   (c) ゆがんだ純ピラミッド格子 : analysis = P の線形場の恒等式 (機械精度) と
#       8 点 / 27 点則の体積一致。底面が平行四辺形でなくなると剛性は
#       有理式の積分近似になる (warp > 0 でないと検査にならない)
#   (d) ゆがんだ遷移格子 : 混在の恒等式 (種別をまたぐ面の節点共有)
#   (e) 非適合な接続 : ピラミッドがあっても、四角形面を三角形 2 枚で覆った
#       面が残っていれば弾く。**節点は全部共有されるので線形場の恒等式では
#       検出できず**、読み込み時の位相の検査だけが防波堤になる
echo "[pyramids] 5-node pyramid (hex-to-tet transition) meshes"
run_case box_pyr
compare "C (pyramid box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
run_case box_hexpyrtet
compare "C (hex er=4 + pyr/tet er=2) [F]" "$(value_of C)" 1.180558376e-13 1e-8
cp "$SRC/nodal_test_pyr.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 nodal_test_pyr.ofe > /dev/null)
res=$(awk '/linear field/ { for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
	/rel. diff/ { for (i = 1; i <= NF; i++) if ($i == "diff") v = $(i+2) }
	/^  warp/ { w = $3 }
	END { if ((e == "") || (v == "")) { printf "NG (no report)"; exit }
	      printf "%s (linear %s, volume %s, warp %s)",
	             (((e + 0) < 1e-12) && ((v + 0) < 1e-12) && ((w + 0) > 1e-3)) ? "OK" : "NG",
	             e, v, w }' "$WORK/ofe.log")
echo "  distorted mesh : the linear-field identity holds : $res"
case "$res" in NG*) status=1 ;; esac
cp "$SRC/nodal_test_hexpyrtet.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 nodal_test_hexpyrtet.ofe > /dev/null)
res=$(awk '/linear field/ { for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
	/^  nodes =/ { for (i = 1; i <= NF; i++) {
	                 if ($i == "tetrahedra") nt = $(i+2)+0
	                 if ($i == "hexahedra") nh = $(i+2)+0
	                 if ($i == "pyramids") np = $(i+2)+0 } }
	END { if (e == "") { printf "NG (no report)"; exit }
	      printf "%s (%s, %d tets + %d hexes + %d pyramids)",
	             (((e + 0) < 1e-12) && (nt > 0) && (nh > 0) && (np > 0)) ? "OK" : "NG",
	             e, nt, nh, np }' "$WORK/ofe.log")
echo "  the identity holds across the hex-pyramid-tet interfaces : $res"
case "$res" in NG*) status=1 ;; esac
# (e) 四角形面を三角形 2 枚で覆った非適合の接続。六面体の上面 (節点 5 6 7 8) を
# 四面体 2 個の三角形が対角線で割って覆う。ピラミッドが 1 個あるので
# 「hex + tet はピラミッド無しでは混ぜられない」の検査は通過し、
# 位相の検査 (elem_face_check) だけがこれを捕まえる
cat > "$WORK/badjunction.msh" <<'EOF'
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
14
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0 0 1
6 1 0 1
7 1 1 1
8 0 1 1
9 0.5 0.5 2
10 3 0 0
11 4 0 0
12 4 1 0
13 3 1 0
14 3.5 0.5 1
$EndNodes
$Elements
5
1 3 2 10 10 1 2 3 4
2 5 2 1 1 1 2 3 4 5 6 7 8
3 4 2 3 3 5 6 7 9
4 4 2 3 3 5 7 8 9
5 7 2 2 2 10 11 12 13 14
$EndElements
EOF
sed 's/^mesh = .*/mesh = badjunction.msh/' "$SRC/box_pyr.ofe" > "$WORK/badj.ofe"
mesh_reject "a quad face covered by two triangles (pyramids present)" badj.ofe
if ! (cd "$WORK" && "$OFE" -n 2 badj.ofe 2>&1 | grep -q "covered by triangular"); then
	echo "  *** it was rejected for a different reason than the face conformity" >&2
	status=1
fi

# 2 次の六面体・角柱 (Gmsh の型 12 = hex27 / 17 = hex20 / 13 = prism18 /
# 18 = prism15)。局所節点の並びは gmsh 4.12.1 の出力を分類して実測した表を使う。
#
#   (a) 長方形格子 (gmsh 製) : 平行平板 C が閉形式と全桁一致 (4 つの型すべて)。
#       電極の中間節点 (9/8 節点四角形・6 節点三角形) の固定もここで通る
#   (b) 2 次場の恒等式 (analysis = P、長方形格子) : **2 次特有の基底 (辺・面・
#       体心) と積分則を検査できるのはこれだけ**。線形場では grad phi が定数に
#       なり、2 次の係数の誤りの多くが素通りする (tet10 と同じ構図)
#   (c) ゆがんだ格子 : 線形場の恒等式 + 3/4 点則の体積一致 (機械精度)
#   (d) 曲がった同軸 (全節点が円筒面上) : 閉形式と 0.1% + **同じ格子の 1 次より
#       30 倍以上良いこと** (2 次が黙って 1 次に退行すると比の検査で落ちる)
#   (e) 次数の混在・1 次の電極面・2 次どうしの種別混在は理由の文字列まで見て弾く
echo "[order-2 hex/prism] 20/27-node hexahedra and 15/18-node prisms"
run_case box_hex2
compare "C (hex27 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
run_case box_hex2s
compare "C (hex20 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
run_case box_prism2
compare "C (prism18 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
run_case box_prism2s
compare "C (prism15 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
# (b) 2 次場の恒等式 (4 つの型すべて。skipped になっていないことも見る)
for m in box_hex2 box_hex2s box_prism2 box_prism2s; do
	sed "s/^mesh = .*/mesh = $m.msh/" "$SRC/nodal_test_hex2r.ofe" > "$WORK/p2q.ofe"
	(cd "$WORK" && "$OFE" -n 2 p2q.ofe > /dev/null)
	res=$(awk '/quadratic field :/ {
			if ($0 ~ /skipped/) { printf "NG (skipped)"; exit }
			for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
		END { if (e == "") { printf "NG (no report)"; exit }
		      printf "%s (%s)", ((e + 0) < 1e-10) ? "OK" : "NG", e }' "$WORK/ofe.log")
	echo "  quadratic-field identity ($m) : $res"
	case "$res" in NG*) status=1 ;; esac
done
# (c) ゆがんだ格子
for f in nodal_test_hex2 nodal_test_prism2; do
	cp "$SRC/$f.ofe" "$WORK/"
	(cd "$WORK" && "$OFE" -n 2 "$f.ofe" > /dev/null)
	res=$(awk '/linear field/ { for (i = 1; i <= NF; i++) if ($i == "error") e = $(i+2) }
		/rel. diff/ { for (i = 1; i <= NF; i++) if ($i == "diff") v = $(i+2) }
		END { if ((e == "") || (v == "")) { printf "NG (no report)"; exit }
		      printf "%s (linear %s, volume %s)",
		             (((e + 0) < 1e-12) && ((v + 0) < 1e-12)) ? "OK" : "NG", e, v }' \
		"$WORK/ofe.log")
	echo "  distorted mesh ($f) : $res"
	case "$res" in NG*) status=1 ;; esac
done
# (d) 曲がった同軸と 1 次の対照 (同じ nr=4, nt=12)
run_case coax_hex2
compare "C' (hex27 coax) [F/m]" "$(value_of C)" 1.063417e-10 0.001
compare "L' (hex27 coax) [H/m]" "$(value_of L)" 2.197225e-07 0.001
c2=$(value_of C)
run_case coax_h1c
c1=$(value_of C)
res=$(awk -v a="$c1" -v b="$c2" 'BEGIN {
	e = 1.0634166122760362e-10
	e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	r = ((e2 > 0) ? (e1 / e2) : 1e9)
	printf "%s (order 1 %+.2f%%, order 2 %+.4f%%, ratio %.0f)",
	       ((r >= 30) ? "OK" : "NG"), e1 * 100, e2 * 100, r }')
echo "  order 2 beats order 1 on the same coarse mesh : $res"
case "$res" in NG*) status=1 ;; esac
# (e) 弾かれるべき入力 (理由の文字列まで見る)
# (次数の違う六面体の混在は弾かず最高次数に上げる。検査は [mixed orders] の (e)。以前ここにあった「1 次の格子に
#  でたらめな節点の 27 節点六面体を足すと弾く」検査は、次数の混在で弾かなくなったので外した)
# 2 次の六面体格子に 1 次の四角形 (電極面) を足す -> 中間節点が固定されない
awk '/^\$Elements/ { print; getline; print $1 + 1
                     print "9998 3 2 10 10 1 2 3 4"; next } { print }' \
    "$SRC/box_hex2.msh" > "$WORK/q1.msh"
sed 's/^mesh = .*/mesh = q1.msh/' "$SRC/box_hex2.ofe" > "$WORK/q1.ofe"
mesh_reject "an order-2 hex mesh with an order-1 quadrilateral" q1.ofe
if ! (cd "$WORK" && "$OFE" -n 2 q1.ofe 2>&1 | grep -q "regenerate the mesh with -order 2"); then
	echo "  *** it was rejected for a different reason than the face order" >&2
	status=1
fi
# 27 節点の六面体に 8 節点四角形の電極面を混ぜる。辺の中間節点は埋まるので
# 「先頭だけ見る」検査は素通りし、面心 32 個が未固定のまま C が -41% になる
# (レビューの実測)。面心の有無まで見る検査だけがこれを弾く
awk '{ if (($2 == 10) && (NF == 14)) { $2 = 16; NF = 13 } print }' \
    "$SRC/box_hex2.msh" > "$WORK/q8.msh"
sed 's/^mesh = .*/mesh = q8.msh/' "$SRC/box_hex2.ofe" > "$WORK/q8.ofe"
mesh_reject "a 27-node hex mesh with 8-node quadrilaterals" q8.ofe
if ! (cd "$WORK" && "$OFE" -n 2 q8.ofe 2>&1 | grep -q "no centre node"); then
	echo "  *** it was rejected for a different reason than the missing centre node" >&2
	status=1
fi
# (2 次の六面体と 1 次の角柱の混在も最高次数に上げる。検査は [mixed high order] の (d) の同軸)

# 六面体・角柱の辺要素 (最低次 Nedelec)。参照要素の辺基底を共変 Piola 変換
# (W = J^-T W^, ∇×W = J c^ / det J) で物理空間に写して数値積分する。
#
#   (a)-(f) 自己検証 (analysis = E) : 四面体と同じ恒等式。**アフィンな要素
#       (平行六面体 / 平行移動の角柱) が前提**で、剛体回転した格子も回す
#       (回転で J の非対角成分が効くため、Piola の取り違えはそこでだけ落ちる)
#   異方性 ν でも回す (等方性だけだとテンソルの非対角項が死んだコードになる)
#   A 解析 : bar_eddy と同じ導体棒・同じ 1 次元厳密解を六面体・角柱で
#       (awall の四角形面 -> 4 辺の固定もここで通る)
#   gauge = 1 が Z を変えないこと (六面体でも tree-cotree が成立すること)
#   混在格子・ピラミッド格子・2 次格子は理由の文字列まで見て弾く
echo "[edge hex/prism] Nedelec edge elements on hexahedra and prisms"
edge_pass() {	# edge_pass <ラベル> <ofe ファイル>
	(cd "$WORK" && "$OFE" -n 2 "$2" > /dev/null 2>&1)
	res=$(awk '/rotational field/ { for (i = 1; i <= NF; i++) if ($i == "err") e = $(i+2) }
		/rotational mass/ { for (i = 1; i <= NF; i++) if ($i == "err") e2 = $(i+2) }
		/self test passed/ { ok = 1 }
		END { if (!ok) { printf "NG (failed)"; exit }
		      printf "OK (curl %s, mass %s)", e, e2 }' "$WORK/ofe.log")
	echo "  $1 : $res"
	case "$res" in NG*) status=1 ;; esac
}
cp "$SRC/edge_test_hex.ofe" "$SRC/edge_test_prism.ofe" "$WORK/"
edge_pass "hexahedra (box)" edge_test_hex.ofe
sed 's/^mesh = .*/mesh = box_hex_rot.msh/' "$SRC/edge_test_hex.ofe" > "$WORK/ehrot.ofe"
edge_pass "hexahedra (rigidly rotated)" ehrot.ofe
edge_pass "prisms (box)" edge_test_prism.ofe
sed 's/^mesh = .*/mesh = box_prism_rot.msh/' "$SRC/edge_test_prism.ofe" > "$WORK/eprot.ofe"
edge_pass "prisms (rigidly rotated)" eprot.ofe
# 異方性 ν (edge_test_aniso と同じテンソル)。等方性では非対角項が実行されない
sed 's|^material = .*|material = 1.0 1e6\nanisomur = 2 4.0 3.0 2.0 0.5 0.7 1.1|' \
    "$SRC/edge_test_hex.ofe" > "$WORK/ehan.ofe"
edge_pass "hexahedra (anisotropic nu)" ehan.ofe
sed 's|^material = .*|material = 1.0 1e6\nanisomur = 2 4.0 3.0 2.0 0.5 0.7 1.1|' \
    "$SRC/edge_test_prism.ofe" > "$WORK/epan.ofe"
edge_pass "prisms (anisotropic nu)" epan.ofe
# A 解析 (bar_eddy と同じ 1 次元厳密解。許容も同じ物理根拠 : 要素寸法/表皮深さ)
# **六面体と角柱は R/L の印字全桁が一致する** (1 次元解が両空間で同じに
# 離散化される) ので、格子が入れ替わっても値の比較では検出できない。
# ofe.log の要素種別の行まで grep して「本当にその種別で解いたこと」を見る
for m in bar_hex bar_prism; do
	cp "$SRC/$m.ofe" "$WORK/"
	case "$m" in
		bar_hex)   kind="trilinear hexahedra" ;;
		bar_prism) kind="6-node prisms" ;;
	esac
	for pair in "1e2 1.37931436e-04 8.37757344e-10 0.005" \
	            "1e4 1.41899129e-04 8.30877185e-10 0.005" \
	            "1e5 3.24859232e-04 5.34552067e-10 0.02"; do
		set -- $pair
		sed "s/^frequency = .*/frequency = $1/" "$SRC/$m.ofe" > "$WORK/bar_kind.ofe"
		(cd "$WORK" && "$OFE" -n 2 bar_kind.ofe > /dev/null && "$OFE_POST" > /dev/null)
		compare "R($m, f=$1) [ohm]" "$(value_of Rf)" "$2" "$4"
		compare "L($m, f=$1) [H]" "$(value_of Lf)" "$3" "$4"
		if grep -q "NOT converged" "$WORK/ofe.log"; then
			echo "  *** A-phi solver did not converge ($m, $1 Hz)" >&2
			status=1
		fi
		if ! grep -q "$kind" "$WORK/ofe.log"; then
			echo "  *** $m did not actually solve on $kind (mesh mixed up?)" >&2
			status=1
		fi
	done
done
# gauge = 1 は Z を変えてはいけない (awall の四角形の辺を優先して木に入れる
# 処理が壊れると端子電圧がずれる)
sed -e "s/^frequency = .*/frequency = 1e4/" "$SRC/bar_hex.ofe" > "$WORK/bh0.ofe"
sed -e "s/^frequency = .*/frequency = 1e4/" -e "s/^awall = 20/awall = 20\ngauge = 1/" \
    "$SRC/bar_hex.ofe" > "$WORK/bh1.ofe"
(cd "$WORK" && "$OFE" -n 2 bh0.ofe > /dev/null && "$OFE_POST" > /dev/null)
r0=$(value_of Rf); l0=$(value_of Lf)
(cd "$WORK" && "$OFE" -n 2 bh1.ofe > /dev/null && "$OFE_POST" > /dev/null)
r1=$(value_of Rf); l1=$(value_of Lf)
res=$(awk -v a="$r0" -v b="$r1" -v c="$l0" -v d="$l1" 'BEGIN {
	er = (a - b) / a; if (er < 0) er = -er
	el = (c - d) / c; if (el < 0) el = -el
	printf "%s (dR %.1e, dL %.1e)", ((er < 1e-6) && (el < 1e-6)) ? "OK" : "NG", er, el }')
echo "  gauge = 1 leaves Z unchanged (hex) : $res"
case "$res" in NG*) status=1 ;; esac
# 弾かれるべき入力
sed -e 's/^analysis = .*/analysis = A/' -e 's/^material = .*/material = 1.0 5.8e7/' \
    -e 's/^solver = .*/frequency = 1e4\nvoltage = 1.0/' "$SRC/box_pyr.ofe" > "$WORK/pyredge.ofe"
mesh_reject "analysis A on a pyramid mesh" pyredge.ofe
if ! (cd "$WORK" && "$OFE" -n 2 pyredge.ofe 2>&1 | grep -q "do not support pyramids"); then
	echo "  *** it was rejected for a different reason than the element kind" >&2
	status=1
fi
sed -e 's/^analysis = .*/analysis = E/' -e 's/^material = .*/material = 1.0 1e6/' \
    "$SRC/box_hex2.ofe" > "$WORK/h2edge.ofe"
mesh_reject "analysis E on an order-2 hex mesh" h2edge.ofe

# ---- 辺要素 : 種別の混在した格子 ----
#
# 最低次 Nedelec の接線トレースは三角形面では三角形の Whitney 空間、四角形面では
# 四角形の最低次 Nedelec 空間になり、どちらもその面の辺の自由度だけで決まるので、
# 面を共有する要素の種別が違っても接線連続性が保たれる。接する面の形が揃う向きに
# しか割れないので、組み合わせで格子の割り方が変わる:
#   六面体 + 角柱 : x で分割 (界面は四角形面)
#   角柱 + 四面体 : z で分割 (界面は三角形面)
#
# **エネルギーの恒等式 (b) (c) (c2) は接線連続性を検査できない。** {a + b×r} の
# 場は要素毎に厳密に補間されるので、隣の要素とトレースが食い違っていても
# 要素積分の和は合ってしまう。面をまたぐ跳びを直接見る (g) が要る。
echo "[edge mixed] Nedelec edge elements on meshes that mix element kinds"
cp "$SRC/edge_mix_hp.ofe" "$SRC/edge_mix_pt.ofe" "$WORK/"
edge_pass "hexahedra + prisms (quadrilateral interface)" edge_mix_hp.ofe
# 種別をまたぐ面が本当にあったか (0 面なら検査は素通りしているのと同じ)
res=$(awk '/tangential continuity/ {
		for (i = 1; i <= NF; i++) {
			if ($i == "faces") { n = $(i-2); m = $(i+1); sub(/^\(/, "", m) }
			if ($i == "=") j = $(i+1)
		}
	}
	END { if (n == "") { printf "NG (no (g) line in the log)"; exit }
	      if (m + 0 < 1) { printf "NG (no face between different kinds)"; exit }
	      if (j + 0 > 1e-10) { printf "NG (jump %s)", j; exit }
	      printf "OK (%s faces, %s between kinds, max jump %s)", n, m, j }' "$WORK/ofe.log")
echo "  (g) tangential trace across kind-mixed faces : $res"
case "$res" in NG*) status=1 ;; esac
edge_pass "prisms + tetrahedra (triangular interface)" edge_mix_pt.ofe
res=$(awk '/tangential continuity/ {
		for (i = 1; i <= NF; i++) if ($i == "faces") { m = $(i+1); sub(/^\(/, "", m) }
	}
	END { printf "%s (%s faces between different kinds)",
	             ((m + 0 > 0) ? "OK" : "NG"), (m == "" ? "no" : m) }' "$WORK/ofe.log")
echo "  (g) the prism-tetrahedron interface is really there : $res"
case "$res" in NG*) status=1 ;; esac
# 渦電流 (A) : 混在させても 1 次元厳密解は同じ。純格子と同じ許容で見る
for pair in "bar_hexprism 0 hexahedra 576 prisms 1152" \
            "bar_prismtet 3456 tetrahedra 0 prisms 1152"; do
	set -- $pair
	m=$1
	cp "$SRC/$m.ofe" "$WORK/"
	sed "s/^frequency = .*/frequency = 1e5/" "$SRC/$m.ofe" > "$WORK/bar_mix.ofe"
	(cd "$WORK" && "$OFE" -n 2 bar_mix.ofe > /dev/null && "$OFE_POST" > /dev/null)
	compare "R($m, f=1e5) [ohm]" "$(value_of Rf)" "3.24859232e-04" "0.02"
	compare "L($m, f=1e5) [H]" "$(value_of Lf)" "5.34552067e-10" "0.02"
	if grep -q "NOT converged" "$WORK/ofe.log"; then
		echo "  *** A-phi solver did not converge ($m)" >&2
		status=1
	fi
	# **本当に混在した格子で解いたか。** 1 次元解はどの空間でも同じように
	# 離散化されるので、値だけでは純格子と区別できない (実測: bar_hexprism の
	# R は純六面体と印字全桁が一致する)。要素数の内訳をログから確かめる
	if ! grep -q "($2 tet / $4 hex / $6 prism)" "$WORK/ofe.log"; then
		echo "  *** $m did not solve on a mixed mesh ($2 tet / $4 hex / $6 prism)" >&2
		status=1
	fi
done
# 弾かれるべき格子
# (1) z に積んだ六面体 + 角柱 : 四角形面を三角形 2 枚が覆う非適合な界面。
#     **どのエネルギー恒等式も (g) も素通りする** (裂けた面はもはや共有面では
#     ないので跳びを比べる相手がいない) ので、位相の検査だけが頼りになる
sed 's/^mesh = .*/mesh = edge_nonconf.msh/' "$SRC/edge_mix_hp.ofe" > "$WORK/enonconf.ofe"
mesh_reject "z-stacked hexahedra + prisms (quadrilateral face covered by triangles)" enonconf.ofe
if ! (cd "$WORK" && "$OFE" -n 2 enonconf.ofe 2>&1 | grep -q "covered by triangular faces"); then
	echo "  *** it was rejected for a different reason than the non-conforming face" >&2
	status=1
fi
# (2) 未マージの格子 (界面の節点を複製したもの、mkmesh.py -mix 1 -tear 1)。
#     これも**すべての恒等式を素通りする** (実測: (a)〜(d) と (g) がすべて
#     機械精度で通り、連結成分が 2 になったことだけが痕跡だった) ので、
#     格子そのものを見るしかない
sed 's/^mesh = .*/mesh = edge_torn.msh/' "$SRC/edge_mix_hp.ofe" > "$WORK/etorn.ofe"
mesh_reject "a torn (unmerged) mesh with duplicated interface nodes" etorn.ofe
if ! (cd "$WORK" && "$OFE" -n 2 etorn.ofe 2>&1 | grep -q "were not merged"); then
	echo "  *** it was rejected for a different reason than the duplicated nodes" >&2
	status=1
fi
# 場の出力 (B / J) を混在格子でも出す。**セル毎のベクトルを作るループは
# 種別ごとに式が違う** (四面体は Whitney の閉形式、六面体・角柱は Piola 変換した
# 要素中心の基底) ので、混在格子ではその分岐を両方通る。
# セル型と、セル数が要素数と種別ごとに合っていることを見る
# (vtkcheck.awk は四面体と三角形しか測れないので体積の恒等式は使えない —
#  六面体・角柱のセルは節点数で三角形と区別が付かず、型を読む必要がある)
awk '/^analysis = A/{print "fieldout = 1"} {print}' "$SRC/bar_hexprism.ofe" > "$WORK/fldmix.ofe"
(cd "$WORK" && "$OFE" -n 2 fldmix.ofe > /dev/null && "$OFE_POST" > /dev/null)
res=$(awk '
	NF == 0            { next }
	/^CELL_TYPES/      { st = "t"; next }
	st == "t" && NF == 1 { n[$1]++ }
	END { printf "%s (%d hexahedra, %d wedges)",
	             ((n[12] == 576) && (n[13] == 1152)) ? "OK" : "NG", n[12], n[13] }' "$WORK/ofe_field.vtk")
echo "  VTK cells are 576 HEXAHEDRON + 1152 WEDGE : $res"
case "$res" in NG*) status=1 ;; esac
# **四面体を含む混在格子でも場を出す。** bar_hexprism は四面体が 0 個なので、
# 角柱の節点勾配を種別ごとの添字 (e - NTet - NHex) で引く経路が NTet > 0 で
# 一度も実行されない。bar_prismtet はそこを通る (四面体 3456 + 角柱 1152)
awk '/^analysis = A/{print "fieldout = 1"} {print}' "$SRC/bar_prismtet.ofe" > "$WORK/fldmix2.ofe"
(cd "$WORK" && "$OFE" -n 2 fldmix2.ofe > /dev/null && "$OFE_POST" > /dev/null)
res=$(awk '
	NF == 0            { next }
	/^CELL_TYPES/      { st = "t"; next }
	st == "t" && NF == 1 { n[$1]++ }
	END { printf "%s (%d tetrahedra, %d wedges)",
	             ((n[10] == 3456) && (n[13] == 1152)) ? "OK" : "NG", n[10], n[13] }' "$WORK/ofe_field.vtk")
echo "  VTK cells are 3456 TETRA + 1152 WEDGE : $res"
case "$res" in NG*) status=1 ;; esac
# **向きの恒等式はこの形状では書けない。** 棒の電流は x 方向なので B は y-z 面内を
# 回り、y 端の近くでは B_z が物理的に 0 でない (実測: 符号つき平均で
# |Bx|/|By| = 9.1e-3、|Bz|/|By| = 4.0e-3、最大値では 6.4e-2 / 2.1e-2)。
# 根拠のあるしきい値が引けないので数値の恒等式は置かず、種別ごとのセル型と
# セル数だけを見る。場の値そのものの恒等式は vtkcheck.awk が六面体・角柱の
# 体積を測れるようになったら (∫J dV = 端子電流の形で) 書ける
# ---- 場の出力の体積 (vtkcheck.awk のセル型対応) ----
#
# vtkcheck.awk は**セル型 (CELL_TYPES) でセルを見分け**、六面体・角柱・
# ピラミッドの体積を**ソルバーと同じ求積**で出す (等パラメトリック写像の体積は
# ∫det J であって四面体に割った和ではない)。節点数だけで見分けていた頃は
# 角柱の 6 節点を 2 次三角形と解釈して面積を足していた (実測: 5.0e-10 のところ
# 6.0e-06)。
#
# 曲げた (warp) 格子でも**外形は変わらない** (内部節点だけ動かすため) ので、
# 六面体・角柱・ピラミッドの ∫det J の和は箱の体積に厳密に一致する。
# 四面体が混ざると直線の要素が角を削るので一致しない — そこは**ソルバー自身が
# 出した体積と突き合わせる** (同じ積分の独立な 2 実装の一致)。
echo "[vtk volume] cell volumes in the field output (hexahedra / prisms / pyramids)"
for pair in "bar_hexprism B_A_re 5.0e-10 mixed hex+prism (A)" \
            "bar_prismtet B_A_re 5.0e-10 mixed prism+tet (A)" \
            "bar_hex B_A_re 5.0e-10 hexahedra (A)" \
            "bar_prism B_A_re 5.0e-10 prisms (A)"; do
	set -- $pair
	awk '/^analysis/{print "fieldout = 1"} {print}' "$SRC/$1.ofe" > "$WORK/vv.ofe"
	(cd "$WORK" && "$OFE" -n 2 vv.ofe > /dev/null 2>&1)
	compare "field volume ($1) [m^3]" "$(vtk vol "$2")" "$3" 1e-12
done
# 曲げた格子 (内部節点だけ動かすので外形 = 箱のまま)
for m in box_hex_warp box_prism_warp box_pyr_warp; do
	# **sed の "i" は GNU 拡張** (BSD sed は i\ + 改行を要求するので macOS で
	# 落ちる)。行の挿入は awk で書く
	awk -v m="$m.msh" '/^mesh = /{print "mesh = " m; next}
		/^analysis/{print "fieldout = 1"} {print}' "$SRC/box_hex.ofe" > "$WORK/vv.ofe"
	(cd "$WORK" && "$OFE" -n 2 vv.ofe > /dev/null 2>&1)
	compare "field volume ($m, warped) [m^3]" "$(vtk vol E_C_port1)" 2.0e-10 1e-12
done
# 四面体が混ざる曲げた格子 : ソルバー自身の体積と突き合わせる
(cd "$WORK" && "$OFE" -n 2 nodal_test_mixed.ofe > /dev/null 2>&1)
vref=$(awk '/volume/ { for (i = 1; i <= NF; i++) if ($i == "=") { print $(i+1); exit } }' "$WORK/ofe.log")
awk -v m="box_mixed_warp.msh" '/^mesh = /{print "mesh = " m; next}
	/^analysis/{print "fieldout = 1"} {print}' "$SRC/box_hex.ofe" > "$WORK/vv.ofe"
(cd "$WORK" && "$OFE" -n 2 vv.ofe > /dev/null 2>&1)
compare "field volume (box_mixed_warp) == the solver's own [m^3]" \
	"$(vtk vol E_C_port1)" "$vref" 1e-12
# ピラミッドは**底面の向きで det J の第 1 項が 0 になる**ので、行列式を
# 1 行で書けていないと 1/3 のセルが 0 になる (実測: 総体積が 2/3 になった)。
# 六面体・四面体だけでは検出できないので、ピラミッド格子を必ず 1 つ通すこと
awk '/^analysis/{print "fieldout = 1"} {print}' "$SRC/box_pyr.ofe" > "$WORK/vv.ofe"
(cd "$WORK" && "$OFE" -n 2 vv.ofe > /dev/null 2>&1)
compare "field volume (box_pyr) [m^3]" "$(vtk vol E_C_port1)" 2.0e-10 1e-12
res=$(awk '/^CELL_TYPES/ { st = "t"; next }
	st == "t" && NF == 1 { n[$1]++ }
	END { printf "%s (%d pyramids)", ((n[14] == 162) ? "OK" : "NG"), n[14] }' "$WORK/ofe_field.vtk")
echo "  the pyramid case really has VTK_PYRAMID cells : $res"
case "$res" in NG*) status=1 ;; esac

# 誤検知の対向検査 : 曲面に載せた 2 次格子は**別々の辺**の中点が同じ点に落ちる
# (実測 coax_p2 で 48 組)。頂点だけを見ているのでこれは弾いてはいけない
res=$( (cd "$WORK" && "$OFE" -n 2 nodal_test_coax.ofe 2>&1 | grep -c "were not merged") || true)
echo "  a curved order-2 mesh with coincident mid-side nodes : $([ "$res" = 0 ] && echo "accepted -> OK" || echo "NG (false positive)")"
[ "$res" = 0 ] || status=1
if ! (cd "$WORK" && "$OFE" -n 2 h2edge.ofe 2>&1 | grep -q "first-order mesh"); then
	echo "  *** it was rejected for a different reason than the order" >&2
	status=1
fi

# MPI (任意依存、WITH_MPI=ON のビルドだけ)。並列化の粒度は**周波数掃引の点の
# 分散**で、各点の演算は直列と完全に同一。したがって検証は閉形式ではなく
# **「mpirun -np N の出力が直列と バイト単位で一致すること」**が恒等式になる
# (同じ物理を別の書き方で 2 通り与えて一致を見る、の並列版)。
#
#   mpirun が無い環境 (Windows / macOS の CI) では節ごと skip する。
#   無効ビルド (WITH_MPI=OFF) では「mpirun 配下に置かれたら弾かれること」だけ
#   見て skip する (黙って N 個の直列プロセスが走ると出力を壊し合う)。
echo "[mpi] frequency-sweep points distributed over processes"
if ! command -v mpirun > /dev/null 2>&1; then
	echo "  mpirun not found : skipped"
else
	# root で走る CI コンテナと少コア環境のための OpenMPI の環境変数
	# (他の MPI 実装は無視するだけなので無害)。**過剰割り当ての変数は OpenMPI 4 と
	# 5 で違う**: 5 (PRRTE) は OMPI_MCA_rmaps_base_oversubscribe を読まない。
	# 実測: CI の runner が OpenMPI 5.0.10 に上がった日から、2 コアの runner で
	# np = 3 が「スロット不足」で起動に失敗し、mpirun がこのスクリプトの
	# シェルごと落として NG の行すら出なかった
	export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
	export OMPI_MCA_rmaps_base_oversubscribe=1
	export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
	cp "$SRC/sweep_plate.ofe" "$WORK/"
	mpimsg=$( (cd "$WORK" && mpirun -np 2 "$OFE" -n 2 sweep_plate.ofe 2>&1) || true)
	case "$mpimsg" in
	*"built without MPI"*)
		echo "  built without MPI : running under mpirun is rejected -> OK (the rest is skipped)"
		;;
	*)
		# 直列の基準 (同じバイナリを mpirun なしで)
		(cd "$WORK" && "$OFE" -n 2 sweep_plate.ofe > /dev/null && "$OFE_POST" > /dev/null)
		cp "$WORK/ofe_sweep.csv" "$WORK/mpi_ref_sweep.csv"
		cp "$WORK/rlc.csv" "$WORK/mpi_ref_rlc.csv"
		cp "$WORK/ofe.out" "$WORK/mpi_ref.out"
		for np in 2 3; do
			rm -f "$WORK/ofe_sweep.csv" "$WORK/rlc.csv" "$WORK/ofe.out"
			(cd "$WORK" && mpirun -np $np "$OFE" -n 2 sweep_plate.ofe > /dev/null 2>&1 \
				&& "$OFE_POST" > /dev/null)
			if cmp -s "$WORK/mpi_ref_sweep.csv" "$WORK/ofe_sweep.csv" \
			&& cmp -s "$WORK/mpi_ref_rlc.csv" "$WORK/rlc.csv" \
			&& cmp -s "$WORK/mpi_ref.out" "$WORK/ofe.out"; then
				echo "  np = $np gives byte-identical ofe_sweep.csv / rlc.csv / ofe.out : OK"
			else
				echo "  np = $np differs from the serial run -> NG" >&2
				status=1
			fi
		done
		# 掃引の無い入力と fieldout = 1 は黙って直列で走らず弾くこと
		cp "$SRC/parallel_plate.ofe" "$WORK/"
		(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 parallel_plate.ofe > /dev/null 2>&1) \
			&& { echo "  no-sweep input under mpirun : accepted -> NG" >&2; status=1; } \
			|| echo "  no-sweep input under mpirun : rejected -> OK"
		if ! grep -q "fewer than two" "$WORK/ofe.log"; then
			echo "  *** it was rejected for a different reason than the missing sweep" >&2
			status=1
		fi
		sed 's/^analysis = /fieldout = 1\nanalysis = /' "$SRC/sweep_plate.ofe" > "$WORK/mfs.ofe"
		(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 mfs.ofe > /dev/null 2>&1) \
			&& { echo "  fieldout = 1 under mpirun : accepted -> NG" >&2; status=1; } \
			|| echo "  fieldout = 1 under mpirun : rejected -> OK"
		if ! grep -q "cannot be combined with MPI" "$WORK/ofe.log"; then
			echo "  *** it was rejected for a different reason than fieldout" >&2
			status=1
		fi
		# F の掃引 (Rf / Lf と Have* の別の組み合わせが pack の復元を通る。
		# sweep_plate は C だけなので、これが無いと Rfmat / Lfmat の復元が
		# 死んだコードになる)
		sed 's/^frequency = .*/frequencysweep = 1e3 1e5 1e7/' \
		    "$SRC/plate_line_ac.ofe" > "$WORK/mfsw.ofe"
		(cd "$WORK" && "$OFE" -n 2 mfsw.ofe > /dev/null 2>&1)
		cp "$WORK/ofe_sweep.csv" "$WORK/mpi_reff.csv"
		rm -f "$WORK/ofe_sweep.csv"
		(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 mfsw.ofe > /dev/null 2>&1)
		if cmp -s "$WORK/mpi_reff.csv" "$WORK/ofe_sweep.csv"; then
			echo "  np = 2 gives a byte-identical eddy-current (F) sweep : OK"
		else
			echo "  np = 2 differs from the serial F sweep -> NG" >&2
			status=1
		fi
		# HDF5 も有効なビルドなら、系列 (ofe_series.h5) の内容一致も見る。
		# ファイルの cmp はオブジェクトヘッダのタイムスタンプで揺れるので
		# h5dump の内容で比べる
		sed 's/^analysis = /hdf5 = 1\nanalysis = /' "$SRC/sweep_plate.ofe" > "$WORK/mh5.ofe"
		h5m=$( (cd "$WORK" && "$OFE" -n 2 mh5.ofe 2>&1) || true)
		if command -v h5dump > /dev/null 2>&1 && [ "${h5m#*needs a build with HDF5}" = "$h5m" ]; then
			(cd "$WORK" && h5dump ofe_series.h5 > mpi_ref.h5dump 2> /dev/null)
			rm -f "$WORK/ofe_series.h5"
			(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 mh5.ofe > /dev/null 2>&1 \
				&& h5dump ofe_series.h5 > mpi_np2.h5dump 2> /dev/null)
			if cmp -s "$WORK/mpi_ref.h5dump" "$WORK/mpi_np2.h5dump"; then
				echo "  np = 2 gives an identical HDF5 series (h5dump) : OK"
			else
				echo "  np = 2 HDF5 series differs from the serial run -> NG" >&2
				status=1
			fi
		else
			echo "  HDF5 series comparison : skipped (no HDF5 build or no h5dump)"
		fi
		# 行ブロックの行列ベクトル積 (mpisolve = 1)。各 rank が crs_spmv の自分の行だけ計算して
		# Allgatherv で揃え、他の演算は全 rank が直列の順で重複して行うので、出力は直列とバイト一致。
		# 構造格子の CG (C)、A-φ の COCG (A)、熱の過渡 (Picard + θ 法)、掃引、場の出力 (rank 0 だけが書く)
		cp "$SRC/parallel_plate.ofe" "$SRC/bar_eddy.ofe" "$SRC/bar_tet.msh" "$SRC/thermal_transient.ofe" "$SRC/sweep_plate.ofe" "$WORK/"
		for c in parallel_plate bar_eddy thermal_transient sweep_plate; do
			awk '/^analysis/{print "mpisolve = 1"} {print}' "$WORK/$c.ofe" > "$WORK/ms_$c.ofe"
		done
		sed 's/^analysis = /fieldout = 1\nanalysis = /' "$WORK/ms_parallel_plate.ofe" > "$WORK/ms_field.ofe"
		for c in parallel_plate bar_eddy thermal_transient sweep_plate field; do
			rm -rf "$WORK/msref"; mkdir -p "$WORK/msref"
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv" "$WORK/ofe_sweep.csv" "$WORK/ofe_field.vtk"
			(cd "$WORK" && "$OFE" -n 2 "ms_$c.ofe" > /dev/null 2>&1) || { echo "  mpisolve $c : the serial run failed -> NG" >&2; status=1; }
			for o in ofe.out ofe_transient.csv ofe_sweep.csv ofe_field.vtk; do
				[ -f "$WORK/$o" ] && cp "$WORK/$o" "$WORK/msref/"
			done
			nps="2 3"; [ "$c" = thermal_transient ] && nps=2
			for np in $nps; do
				rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv" "$WORK/ofe_sweep.csv" "$WORK/ofe_field.vtk"
				(cd "$WORK" && mpirun -np $np "$OFE" -n 2 "ms_$c.ofe" > /dev/null 2>&1) || true
				r=OK; nf=0
				for o in ofe.out ofe_transient.csv ofe_sweep.csv ofe_field.vtk; do
					if [ -f "$WORK/msref/$o" ]; then
						nf=$((nf + 1))
						cmp -s "$WORK/msref/$o" "$WORK/$o" || r=NG
					fi
				done
				grep -q "matrix-vector products split by rows over $np ranks" "$WORK/ofe.log" || r=NG
				if [ "$r" = OK ] && [ "$nf" -ge 1 ]; then
					echo "  mpisolve = 1, np = $np : $c gives byte-identical output ($nf files) : OK"
				else
					echo "  mpisolve = 1, np = $np : $c differs from the serial run (or not split) -> NG" >&2
					status=1
				fi
			done
		done
		# 行列の分散 (mpisolve = 2、C / L / R / M / F / A と熱連成。自己検証の E / P だけは入力エラー)。各 rank は行ブロック (構造格子は x 方向に切った板) の行だけを
		# 組み立てて持ち、対角・反作用は持ち主から集める。行ごとの和の順は直列と同じなので出力はバイト一致。
		# 「分散したこと」は rank 0 が持つ行の数 (= n / np) で見る (分散し忘れても答えは同じになるので)
		cp "$SRC/resistor_bar.ofe" "$SRC/coupled_microstrip.ofe" "$WORK/"
		for c in parallel_plate resistor_bar coupled_microstrip sweep_plate; do
			awk '/^analysis/{print "mpisolve = 2"} {print}' "$WORK/$c.ofe" > "$WORK/md_$c.ofe"
			rm -f "$WORK/ofe.out" "$WORK/ofe_sweep.csv"; rm -rf "$WORK/mdref"; mkdir -p "$WORK/mdref"
			(cd "$WORK" && "$OFE" -n 2 "md_$c.ofe" > /dev/null 2>&1) || { echo "  mpisolve = 2 $c : the serial run failed -> NG" >&2; status=1; }
			for o in ofe.out ofe_sweep.csv; do [ -f "$WORK/$o" ] && cp "$WORK/$o" "$WORK/mdref/"; done
			for np in 2 3; do
				rm -f "$WORK/ofe.out" "$WORK/ofe_sweep.csv"
				(cd "$WORK" && mpirun -np $np "$OFE" -n 2 "md_$c.ofe" > /dev/null 2>&1) || true
				r=OK; nf=0
				for o in ofe.out ofe_sweep.csv; do
					if [ -f "$WORK/mdref/$o" ]; then nf=$((nf + 1)); cmp -s "$WORK/mdref/$o" "$WORK/$o" || r=NG; fi
				done
				d=$(awk -v np=$np '/the matrix is distributed by rows/ { for (i = 1; i <= NF; i++) if ($i == "rows") { split($(i + 1), a, "\\.\\."); n = $(i + 3) }
					ok = ((a[1] == 0) && (a[2] + 1 == int(n / np))) ? "split" : "bad"; print ok, a[2] + 1, "of", n; exit }' "$WORK/ofe.log")
				case "$d" in split*) ;; *) r=NG ;; esac
				if [ "$r" = OK ] && [ "$nf" -ge 1 ]; then
					echo "  mpisolve = 2, np = $np : $c gives byte-identical output ($nf files), rank 0 holds ${d#split } rows : OK"
				else
					echo "  mpisolve = 2, np = $np : $c differs from the serial run (or not distributed: $d) -> NG" >&2
					status=1
				fi
			done
		done
		# 非構造格子 (四面体・3 次の曲がった同軸・六面体 + ピラミッド + 四面体の混在) の C、四面体の R、1 次の A-φ
		# (四面体・空気層・六面体。実部 Kr と虚部 Ki を同じ分け方で詰める)。A は節点の行列と A-φ の行列の 2 つが分散する
		cp "$SRC/box_tet.ofe" "$SRC/box_tet.msh" "$SRC/coax_p3.ofe" "$SRC/coax_p3.msh" "$SRC/box_hexpyrtet.ofe" "$SRC/box_hexpyrtet.msh" \
		   "$SRC/bar_air.ofe" "$SRC/bar_air.msh" "$SRC/bar_hex.ofe" "$SRC/bar_hex.msh" "$WORK/"
		grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^frequency/{next} /^awall/{next} /^analysis/{print "analysis = R"; next} {print}' > "$WORK/bar_r.ofe"
		for c in box_tet coax_p3 box_hexpyrtet bar_r bar_eddy bar_air bar_hex; do
			awk '/^analysis/{print "mpisolve = 2"} {print}' "$WORK/$c.ofe" > "$WORK/mu_$c.ofe"
			rm -f "$WORK/ofe.out"
			(cd "$WORK" && "$OFE" -n 2 "mu_$c.ofe" > /dev/null 2>&1) || { echo "  mpisolve = 2 $c : the serial run failed -> NG" >&2; status=1; }
			cp "$WORK/ofe.out" "$WORK/mu_ref.out" 2> /dev/null || true
			nd=1; case "$c" in bar_eddy|bar_air|bar_hex) nd=2 ;; esac
			rm -f "$WORK/ofe.out"
			(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 "mu_$c.ofe" > /dev/null 2>&1) || true
			k=$(grep -c "the matrix is distributed by rows" "$WORK/ofe.log" || true)
			if cmp -s "$WORK/mu_ref.out" "$WORK/ofe.out" && [ "$k" = "$nd" ]; then
				echo "  mpisolve = 2, np = 2 : $c (unstructured) gives a byte-identical ofe.out, $k distributed matrices : OK"
			else
				echo "  mpisolve = 2, np = 2 : $c (unstructured) differs from the serial run (or $k of $nd matrices distributed) -> NG" >&2
				status=1
			fi
		done
		# 2 次 (eddy3d_p2.c) と 3 次 (eddy3d_pr.c) の A、1 次と 3 次の A の線形の過渡解析 (K と M を分散、ofe_transient.csv も一致)
		cp "$SRC/bar_tet2.msh" "$SRC/bar_ho3.msh" "$WORK/"
		for m in bar_tet2 bar_ho3; do
			grep -v '^#' "$WORK/bar_eddy.ofe" | awk -v m="$m.msh" '/^mesh/{print "mesh = " m; next} /^analysis/{print "mpisolve = 2"} {print}' > "$WORK/mu_$m.ofe"
		done
		grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^frequency/{next} /^analysis/{print "mpisolve = 2"; print "transient = 2e-4 8"} {print}' > "$WORK/mu_tr.ofe"
		awk '/^mesh/{print "mesh = bar_ho3.msh"; next} {print}' "$WORK/mu_tr.ofe" > "$WORK/mu_tr3.ofe"
		for c in bar_tet2 bar_ho3 tr tr3; do
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv" "$WORK/mu_ref.csv"
			(cd "$WORK" && "$OFE" -n 2 "mu_$c.ofe" > /dev/null 2>&1) || { echo "  mpisolve = 2 $c : the serial run failed -> NG" >&2; status=1; }
			cp "$WORK/ofe.out" "$WORK/mu_ref.out" 2> /dev/null || true
			[ -f "$WORK/ofe_transient.csv" ] && cp "$WORK/ofe_transient.csv" "$WORK/mu_ref.csv"
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv"
			(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 "mu_$c.ofe" > /dev/null 2>&1) || true
			k=$(grep -c "the matrix is distributed by rows" "$WORK/ofe.log" || true)
			r=OK
			cmp -s "$WORK/mu_ref.out" "$WORK/ofe.out" || r=NG
			if [ -f "$WORK/mu_ref.csv" ]; then cmp -s "$WORK/mu_ref.csv" "$WORK/ofe_transient.csv" || r=NG; fi
			if [ "$r" = OK ] && [ "$k" = 2 ]; then
				echo "  mpisolve = 2, np = 2 : A ($c) gives byte-identical output, $k distributed matrices : OK"
			else
				echo "  mpisolve = 2, np = 2 : A ($c) differs from the serial run (or $k of 2 matrices distributed) -> NG" >&2
				status=1
			fi
		done
		# 断面 2 次元の M (線形、構造格子・三角形) と F (構造格子・三角形・B-H の実効磁気抵抗率法)、熱連成 (構造格子・
		# 四面体・過渡解析)。M / F は K と M (σ) を、熱連成は電気・熱・熱容量の行列を同じ分け方で詰める
		cp "$SRC/plate_line_dc.ofe" "$SRC/plate2d_dc.ofe" "$SRC/plate2d.msh" "$SRC/plate_line_ac.ofe" "$SRC/plate2d_ac.ofe" \
		   "$SRC/plate_line_bh_ac.ofe" "$SRC/thermal_bar.ofe" "$SRC/thermal_tet.ofe" "$SRC/thermal_transient.ofe" "$WORK/"
		for c in plate_line_dc plate2d_dc plate_line_ac plate2d_ac plate_line_bh_ac thermal_bar thermal_tet thermal_transient; do
			[ -f "$WORK/$c.ofe" ] || { echo "  mpisolve = 2 $c : no input -> NG" >&2; status=1; continue; }
			m=$(awk '/^mesh/{print $3}' "$WORK/$c.ofe"); [ -n "$m" ] && cp "$SRC/$m" "$WORK/" 2> /dev/null
			awk '/^analysis/{print "mpisolve = 2"} {print}' "$WORK/$c.ofe" > "$WORK/mv_$c.ofe"
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv" "$WORK/mv_ref.csv"
			(cd "$WORK" && "$OFE" -n 2 "mv_$c.ofe" > /dev/null 2>&1) || { echo "  mpisolve = 2 $c : the serial run failed -> NG" >&2; status=1; }
			cp "$WORK/ofe.out" "$WORK/mv_ref.out" 2> /dev/null || true
			[ -f "$WORK/ofe_transient.csv" ] && cp "$WORK/ofe_transient.csv" "$WORK/mv_ref.csv"
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv"
			(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 "mv_$c.ofe" > /dev/null 2>&1) || true
			k=$(grep -c "the matrix is distributed by rows" "$WORK/ofe.log" || true)
			r=OK
			cmp -s "$WORK/mv_ref.out" "$WORK/ofe.out" || r=NG
			if [ -f "$WORK/mv_ref.csv" ]; then cmp -s "$WORK/mv_ref.csv" "$WORK/ofe_transient.csv" || r=NG; fi
			if [ "$r" = OK ] && [ "$k" -ge 2 ]; then
				echo "  mpisolve = 2, np = 2 : $c gives byte-identical output, $k distributed matrices : OK"
			else
				echo "  mpisolve = 2, np = 2 : $c differs from the serial run (or $k matrices distributed) -> NG" >&2
				status=1
			fi
		done
		# 非線形の M (Newton の B-H・J-A・渦電流 + J-A の過渡) と、非線形の A の過渡解析 (1 次の J-A、2 次の B-H。接線の
		# 行列 Kt も同じ分け方で詰める)
		cp "$SRC/plate_line_bh.ofe" "$SRC/plate_line_ja.ofe" "$SRC/lamination_step.ofe" "$SRC/bar_tet2c.msh" "$WORK/"
		grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^frequency/{next} /^voltage/{print "voltage = 1e-4"; next}
			/^analysis/{print "transient = 2e-4 10"; print "ja = 2 1.6e6 1100 1.6e-3 400 0.2"; print "nlsolver = 200 1e-9 1"} {print}' > "$WORK/a1ja.ofe"
		grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_tet2c.msh"; next} /^frequency/{next} /^voltage/{print "voltage = 1e-4"; next}
			/^analysis/{print "transient = 2e-4 10"; print "bh = 2 100 0.5"; print "bh = 2 1000 1.0"; print "nlsolver = 50 1e-9 1"} {print}' > "$WORK/a2bh.ofe"
		for c in plate_line_bh plate_line_ja lamination_step a1ja a2bh; do
			awk '/^analysis/{print "mpisolve = 2"} {print}' "$WORK/$c.ofe" > "$WORK/mn_$c.ofe"
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv" "$WORK/mn_ref.csv"
			(cd "$WORK" && "$OFE" -n 2 "mn_$c.ofe" > /dev/null 2>&1) || { echo "  mpisolve = 2 $c : the serial run failed -> NG" >&2; status=1; }
			cp "$WORK/ofe.out" "$WORK/mn_ref.out" 2> /dev/null || true
			[ -f "$WORK/ofe_transient.csv" ] && cp "$WORK/ofe_transient.csv" "$WORK/mn_ref.csv"
			rm -f "$WORK/ofe.out" "$WORK/ofe_transient.csv"
			(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 "mn_$c.ofe" > /dev/null 2>&1) || true
			k=$(grep -c "the matrix is distributed by rows" "$WORK/ofe.log" || true)
			r=OK
			cmp -s "$WORK/mn_ref.out" "$WORK/ofe.out" || r=NG
			if [ -f "$WORK/mn_ref.csv" ]; then cmp -s "$WORK/mn_ref.csv" "$WORK/ofe_transient.csv" || r=NG; fi
			if [ "$r" = OK ] && [ "$k" -ge 1 ]; then
				echo "  mpisolve = 2, np = 2 : nonlinear $c gives byte-identical output, $k distributed matrices logged : OK"
			else
				echo "  mpisolve = 2, np = 2 : nonlinear $c differs from the serial run (or $k matrices distributed) -> NG" >&2
				status=1
			fi
		done
		# 自己検証 (E / P) は分散しない (入力エラー)
		awk '/^analysis/{print "mpisolve = 2"} {print}' "$SRC/edge_test.ofe" > "$WORK/mdx_e.ofe"
		cp "$SRC/box_tet.msh" "$WORK/"
		(cd "$WORK" && mpirun -np 2 "$OFE" -n 2 mdx_e.ofe > /dev/null 2>&1) \
			&& { echo "  mpisolve = 2 with analysis E : accepted -> NG" >&2; status=1; } \
			|| echo "  mpisolve = 2 with analysis E : rejected -> OK"
		grep -q "mpisolve = 2 (the matrix distributed over the ranks) supports analysis C / L / R / M / F / A" "$WORK/ofe.log" \
			|| { echo "  *** it was rejected for a different reason" >&2; status=1; }
		;;
	esac
fi
# mpisolve = 1 は MPI の無いビルドでは黙って直列に落ちず入力エラー (MPI のビルドでは直列でも走る)
awk '/^analysis/{print "mpisolve = 1"} {print}' "$SRC/parallel_plate.ofe" > "$WORK/msx.ofe"
msmsg=$( (cd "$WORK" && "$OFE" -n 2 msx.ofe 2>&1) || true)
case "$msmsg" in
*"needs a build with -DWITH_MPI=ON"*)
	echo "  [mpi] mpisolve = 1 in a build without MPI : rejected -> OK" ;;
*)
	if grep -q "normal end" "$WORK/ofe.log"; then
		echo "  [mpi] mpisolve = 1 without mpirun (MPI build) : runs serially -> OK"
	else
		echo "  [mpi] mpisolve = 1 : neither rejected nor run -> NG" >&2; status=1
	fi ;;
esac
awk '/^analysis/{print "mpisolve = 2"} {print}' "$SRC/parallel_plate.ofe" > "$WORK/msx2.ofe"
msmsg=$( (cd "$WORK" && "$OFE" -n 2 msx2.ofe 2>&1) || true)
case "$msmsg" in
*"needs a build with -DWITH_MPI=ON"*)
	echo "  [mpi] mpisolve = 2 in a build without MPI : rejected -> OK" ;;
*)
	if grep -q "normal end" "$WORK/ofe.log"; then
		echo "  [mpi] mpisolve = 2 without mpirun (MPI build) : runs serially -> OK"
	else
		echo "  [mpi] mpisolve = 2 : neither rejected nor run -> NG" >&2; status=1
	fi ;;
esac

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
# ピラミッドを含む混在格子 (型 2/3/4/5/7 が 1 つのファイルに同居する)。
# バイナリでは節点数が型の表から決まるので、型 7 = 5 節点の表の誤りはここで
# ずれる。ファイルは本物の gmsh (4.12.1) に box_hexpyrtet.msh を変換させたもの
bin_same "hex-pyramid-tet mix, binary 2.2" box_hexpyrtet.ofe box_hexpyrtet_bin.msh box_hexpyrtet_bin_22.msh
bin_same "hex-pyramid-tet mix, binary 4.1" box_hexpyrtet.ofe box_hexpyrtet_bin.msh box_hexpyrtet_bin_41.msh
# 2 次の六面体 (型 12 と 9 節点四角形 = 型 10)。バイナリでは節点数が型の表から
# 決まるので、27 / 9 の表の誤りはここでずれる。ファイルは本物の gmsh 製
bin_same "order-2 hexahedra, binary 2.2" box_hex2.ofe box_hex2_bin.msh box_hex2_bin_22.msh
bin_same "order-2 hexahedra, binary 4.1" box_hex2.ofe box_hex2_bin.msh box_hex2_bin_41.msh
# ASCII 4.1 のピラミッド (read_elements_v41 の節点数表と $Entities 経由の
# 物理タグ)。比較の基準は**同じ gmsh 変換で出した 2.2 ASCII** (gmsh は変換の
# たびに節点を振り直すので、元の box_hexpyrtet.msh と直接比べると丸めの順序が
# 変わって一致しない)
sed 's/^mesh = .*/mesh = box_hexpyrtet_bin.msh/' "$SRC/box_hexpyrtet.ofe" > "$WORK/pm_a.ofe"
sed 's/^mesh = .*/mesh = box_hexpyrtet_41.msh/' "$SRC/box_hexpyrtet.ofe" > "$WORK/pm_b.ofe"
(cd "$WORK" && "$OFE" -n 2 pm_a.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/pm_a.csv"
(cd "$WORK" && "$OFE" -n 2 pm_b.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/pm_b.csv"
if cmp -s "$WORK/pm_a.csv" "$WORK/pm_b.csv"; then
	echo "  hex-pyramid-tet mix, ASCII 4.1 : identical -> OK"
else
	echo "  hex-pyramid-tet mix, ASCII 4.1 : differ -> NG" >&2
	diff "$WORK/pm_a.csv" "$WORK/pm_b.csv" | head -4 >&2
	status=1
fi
# 2 次の六面体の ASCII 4.1 (read_elements_v41 の節点数表)
sed 's/^mesh = .*/mesh = box_hex2_bin.msh/' "$SRC/box_hex2.ofe" > "$WORK/h2_a.ofe"
sed 's/^mesh = .*/mesh = box_hex2_41.msh/' "$SRC/box_hex2.ofe" > "$WORK/h2_b.ofe"
(cd "$WORK" && "$OFE" -n 2 h2_a.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/h2_a.csv"
(cd "$WORK" && "$OFE" -n 2 h2_b.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/h2_b.csv"
if cmp -s "$WORK/h2_a.csv" "$WORK/h2_b.csv"; then
	echo "  order-2 hexahedra, ASCII 4.1 : identical -> OK"
else
	echo "  order-2 hexahedra, ASCII 4.1 : differ -> NG" >&2
	status=1
fi

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
cp "$SRC/plate2d_p2.ofe" "$WORK/"
mx_same "mixed triangle orders in a 2-D mesh (one 3-node triangle among 6-node ones)" plate2d_p2.ofe p2mix.ofe

# ---- 断面 2 次元の 3 次三角形 (10 節点、gmsh の tri10) ----
#
# 局所節点の並びは gmsh 4.15.2 の tri10 を重心座標で分類した実測の表 (unstruct.c の
# TRI3_A)。形状関数は 3 次四面体と同じ多重指数の Lagrange の式、積分は Duffy + 4 点。
#   (a) DC : 導体内の Az は 2 次なので Ldc・Rs は P2 と同じく閉形式に厳密
#   (b) F (10 MHz、t/δ = 2.4) : 1 次元厳密解に R +6.7e-8 / L +1.8e-7 (実測)。許容 1e-6。
#       **3 次特有の基底を検査できるのは表皮効果の解だけ** (DC の 2 次の Az は P2 の
#       部分空間で厳密に表せる)。同じ頂点の 2 次格子 (plate2d_p3o2.msh、gmsh の -order 2)
#       より R の誤差が 100 倍以上小さいこと (実測 7200 倍) を assert する。L は 43 倍 (> 10)
#   (c) 2.2 と 4.1 (同じ gmsh の変換で出したもの) で rlc.csv が完全一致
#   (d) 場の出力 : 面積と損失の恒等式、セルは頂点だけの三角形 (VTK_TRIANGLE、3 節点)
#   (e) 次数の混在は弾く
echo "[plate2d_p3] third-order (10-node) triangles on the 2-D mesh (gmsh tri10)"
run_case plate2d_p3
compare "Ldc (P3) [H/m]" "$(value_of Ldc)" 2.932153e-07 1e-5
compare "Rs (P3) [ohm/m]" "$(value_of Rs)" 6.896552e-01 1e-6
p3z() {	# p3z <格子> <周波数> : plate2d_ac.ofe の格子と周波数を差し替えて解く
	sed -e "s/plate2d.msh/$1/" -e "s/^frequency = .*/frequency = $2/" "$SRC/plate2d_ac.ofe" > "$WORK/p3ac.ofe"
	(cd "$WORK" && "$OFE" -n 2 p3ac.ofe > /dev/null && "$OFE_POST" > /dev/null)
}
#     閉形式 (awk で全桁) : R = 2k (S + s)、X = 2k (S - s) + ωμ0 d/W、
#     k = 1/(δσW (cosh 2a - cos 2a))、S = sinh 2a、s = sin 2a、a = t/δ
set -- $(awk 'BEGIN { pi = atan2(0, -1); mu = 4e-7 * pi; s = 5.8e7; W = 1e-3; t = 0.05e-3; d = 0.2e-3
	om = 2 * pi * 1e7; de = sqrt(2 / (om * mu * s)); a = t / de
	S = (exp(2 * a) - exp(-2 * a)) / 2; C = (exp(2 * a) + exp(-2 * a)) / 2
	k = 1 / (de * s * W * (C - cos(2 * a)))
	printf "%.12e %.12e", 2 * k * (S + sin(2 * a)), (2 * k * (S - sin(2 * a)) + om * mu * d / W) / om }')
rex=$1; lex=$2
p3z plate2d_p3.msh 1e7
r3=$(value_of Rf); l3=$(value_of Lf)
compare "R(f = 10 MHz, P3) vs coth closed form [ohm/m]" "$r3" "$rex" 1e-6
compare "L(f = 10 MHz, P3) vs coth closed form [H/m]" "$l3" "$lex" 1e-6
p3z plate2d_p3o2.msh 1e7
res=$(awk -v a="$(value_of Rf)" -v b="$r3" -v c="$(value_of Lf)" -v d="$l3" -v r="$rex" -v l="$lex" 'BEGIN {
	e2 = (a - r) / r; if (e2 < 0) e2 = -e2
	e3 = (b - r) / r; if (e3 < 0) e3 = -e3
	f2 = (c - l) / l; if (f2 < 0) f2 = -f2
	f3 = (d - l) / l; if (f3 < 0) f3 = -f3
	ok = (e2 > 100 * e3) && (f2 > 10 * f3)
	printf "%s (R : P2 %.2e vs P3 %.2e, L : P2 %.2e vs P3 %.2e)", ok ? "OK" : "NG", e2, e3, f2, f3 }')
echo "  order 3 beats order 2 on the same vertices (10 MHz) : $res"
case "$res" in NG*) status=1 ;; esac
#     1 kHz (表皮効果なし) では DC と同じく P2 / P3 とも厳密
p3z plate2d_p3.msh 1e3
compare "R(f = 1 kHz, P3) [ohm/m]" "$(value_of Rf)" 6.896552e-01 1e-6
compare "L(f = 1 kHz, P3) [H/m]" "$(value_of Lf)" 2.932153e-07 1e-5
# (c) 2.2 と 4.1
sed 's/^mesh = .*/mesh = plate2d_p3_41.msh/' "$SRC/plate2d_p3.ofe" > "$WORK/p3_41.ofe"
(cd "$WORK" && "$OFE" -n 2 p3_41.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/p3_41.csv"
run_case plate2d_p3
grep -v '^title' "$WORK/rlc.csv" > "$WORK/p3_22.csv"
if cmp -s "$WORK/p3_22.csv" "$WORK/p3_41.csv"; then echo "  10-node triangles, ASCII 2.2 vs 4.1 : identical -> OK"
else echo "  10-node triangles, ASCII 2.2 vs 4.1 : differ -> NG" >&2; status=1; fi
# (d) 場の出力
sed -e 's/plate2d.msh/plate2d_p3.msh/' "$SRC/plate2d_ac.ofe" | awk '/^analysis/{print "fieldout = 1"} {print}' > "$WORK/p3fld.ofe"
(cd "$WORK" && "$OFE" -n 2 p3fld.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "field area (P3 2-D) [m^2]" "$(vtk vol J_F_re_port1)" 3.0e-7 1e-9
pf=$(awk -v a="$(vtk int2 J_F_re_port1)" -v b="$(vtk int2 J_F_im_port1)" \
	'BEGIN{ printf "%.10e", (a + b) / (2 * 5.8e7) }')
pt=$(awk -v R="$(value_of Rf)" -v L="$(value_of Lf)" \
	'BEGIN{ om = 2 * 3.14159265358979324 * 1e3; X = om * L
	        printf "%.10e", 0.5 * R / ((R * R) + (X * X)) }')
compare "ohmic loss from the field (P3 2-D) [W/m]" "$pf" "$pt" 1e-4
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
echo "  VTK cells are vertex triangles (type 5, 3 nodes) : $res"
case "$res" in NG*) status=1 ;; esac
# (e) 次数の混在 (tri10 を 1 つだけ頂点の tri3 にする) は 3 次に上げて解く。直線の三角形なので一様な 3 次と同じ離散問題
#     (tri10 の先頭 6 節点を tri6 にすると辺の節点が 1/3 の位置の歪んだ要素になり、上げてもその形のままなので一致しない)
awk '/^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 }
	e && ($2 == 21) && !done { s = $1 " 2"; for (i = 3; i <= 3 + $3 + 3; i++) s = s " " $i; print s; done = 1; next } { print }' \
	"$SRC/plate2d_p3.msh" > "$WORK/mix3d2.msh"
sed 's/^mesh = .*/mesh = mix3d2.msh/' "$SRC/plate2d_p3.ofe" > "$WORK/p3mix.ofe"
cp "$SRC/plate2d_p3.ofe" "$WORK/"
mx_same "mixed triangle orders in a 2-D mesh (one 3-node triangle among 10-node ones)" plate2d_p3.ofe p3mix.ofe

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

# 温度係数の表 (epstemptable / mutemptable / bhtemptable = <mid> <T1> <f1> ...)。
# 1 次係数の代わりに折れ線の倍率を使い、キュリー点近傍の非線形な落ち込みを書く。
# 適用箇所は 1 次係数と同じ (εr / μr は読み出し時、B-H は入力時) で、temp_factor() に
# 集約してある
echo "[temptable] piecewise-linear temperature factors (Curie-point drop)"
tt_run() {	# tt_run <dir> : $WORK/<dir>/x.ofe を解いて ofe_post まで
	(cd "$WORK/$1" && "$OFE" -n 2 x.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1)
}
tt_val() { awk -F, -v name="$2" '$1 == name && NF == 2 { f = 1; next } f && $1 == "1" { print $2; exit }' "$WORK/$1/rlc.csv"; }
#   (a) 1 次と同じ直線の表は 1 次係数と一致する (μr、T = 70 / 120)
for t in 70 120; do
	mkdir -p "$WORK/ttl$t" "$WORK/ttt$t"
	sed "s/^temperature = .*/temperature = $t/" "$SRC/temp_mur.ofe" > "$WORK/ttl$t/x.ofe"
	sed "s/^temperature = .*/temperature = $t/; s/^mutempco = .*/mutemptable = 2 20 1 120 0.7/" \
		"$SRC/temp_mur.ofe" > "$WORK/ttt$t/x.ofe"
	tt_run "ttl$t"; tt_run "ttt$t"
	compare "mutemptable (a straight line) == mutempco, T = $t [H/m]" "$(tt_val ttt$t Ldc)" "$(tt_val ttl$t Ldc)" 1e-12
done
#   (b) キュリー点 : 100 degC まで 1、110 degC で 0.25 (μr0 = 4 なので μr = 1)。
#       Tc を越えると**真空と同じ問題**になる : μr = 4 x 0.25 は浮動小数点でも厳密に 1
#       なので、mur = 1 の入力と ofe.out がバイト一致する。表の外 (300 degC) は端の値
#       のまま (バイト一致)、表の内側 (105 degC) は f = 0.625 -> μr = 2.5 の閉形式
tt_curie() {	# tt_curie <dir> <T>
	mkdir -p "$WORK/$1"
	sed "s/^temperature = .*/temperature = $2/; s/^mutempco = .*/mutemptable = 2 20 1 100 1 110 0.25/" \
		"$SRC/temp_mur.ofe" > "$WORK/$1/x.ofe"
	tt_run "$1"
}
tt_curie ttc150 150
tt_curie ttc300 300
tt_curie ttc105 105
tt_curie ttc0 0
mkdir -p "$WORK/ttvac" "$WORK/tt20"
sed "s/^temperature = .*/temperature = 150/; s/^mur = .*/mur = 2 1.0/; /^mutempco/d" "$SRC/temp_mur.ofe" > "$WORK/ttvac/x.ofe"
sed "s/^temperature = .*/temperature = 20/; /^mutempco/d" "$SRC/temp_mur.ofe" > "$WORK/tt20/x.ofe"
tt_run ttvac; tt_run tt20
tt_same() {	# tt_same <ラベル> <dir1> <dir2>
	if [ -s "$WORK/$2/ofe.out" ] && cmp -s "$WORK/$2/ofe.out" "$WORK/$3/ofe.out"; then echo "  $1 : byte-identical -> OK"
	else echo "  $1 : differ -> NG" >&2; status=1; fi
}
tt_same "above Tc (150 degC) == mur = 1 (vacuum)" ttc150 ttvac
tt_same "beyond the table (300 degC) == its last point (150 degC)" ttc300 ttc150
tt_same "below the table (0 degC) == its first point (no temperature factor)" ttc0 tt20
compare "inside the drop (105 degC, mur = 2.5) [H/m]" "$(tt_val ttc105 Ldc)" \
	"$(awk 'BEGIN{ pi = 3.14159265358979324; mu0 = 4 * pi * 1e-7
		printf "%.10e", mu0 * ((2.5 * 0.2e-3) + (2 * 0.05e-3 / 3)) / 1e-3 }')" 0.002
#   (c) εr の表 (2 区間)。C は εr に厳密に比例する。T = 170 は折れ目の先の区間で
#       f = 1.5 + 0.5 (1.2 - 1.5) = 1.35。T = 0 は表の手前で最初の点の値 1.1
#       (**最初の点を 1 以外にしておかないと、手前を 1 にする誤りが等価になる**)
for pair in "170 1.35" "0 1.1"; do
	set -- $pair
	mkdir -p "$WORK/tte$1"
	sed "s/^temperature = .*/temperature = $1/; s/^epstempco = .*/epstemptable = 2 20 1.1 120 1.5 220 1.2/" \
		"$SRC/temp_material.ofe" > "$WORK/tte$1/x.ofe"
	tt_run "tte$1"
	compare "epstemptable : C($1 degC) = $2 C0 [F]" "$(tt_val tte$1 C)" \
		"$(awk -v f="$2" 'BEGIN { printf "%.10e", 1.7708375600e-13 * f }')" 1e-6
done
#   (d) B-H の表 : B 軸を 0.8 倍すると (L - L_int) がちょうど 0.8 倍
for t in 20 120; do
	mkdir -p "$WORK/ttb$t"
	sed 's/^analysis = M/bhtemptable = 2 20 1 120 0.8\ntemperature = '"$t"'\nanalysis = M/' \
	    "$SRC/plate_line_bh.ofe" > "$WORK/ttb$t/x.ofe"
	tt_run "ttb$t"
done
compare "bhtemptable : (L - L_int) ratio at T = 120 vs 20" \
	"$(awk -v a="$(tt_val ttb120 Ldc)" -v b="$(tt_val ttb20 Ldc)" 'BEGIN{ pi = 3.14159265358979324; mu0 = 4 * pi * 1e-7
		li = mu0 * 2 * 0.05e-3 / (3 * 1e-3)
		printf "%.10e", (a - li) / (b - li) }')" 0.8 1e-5
#   (d2) 抵抗率の表 (tempcotable、σ(T) = σ0 / f(T))。R は抵抗率に厳密に比例する。
#       temp_resistor.ofe (85 degC) : 直線の表は tempco と一致、折れ線の表は
#       f(85) = 1.2 + (35/50) 0.8 = 1.76 で R = 1.76 R(20 degC, 補正なし)
for c in lin tab ref; do mkdir -p "$WORK/tts_$c"; done
cp "$SRC/temp_resistor.ofe" "$WORK/tts_lin/x.ofe"
sed "s/^tempco = .*/tempcotable = 2 20 1 50 1.2 100 2.0/" "$SRC/temp_resistor.ofe" > "$WORK/tts_tab/x.ofe"
sed "/^tempco/d; s/^temperature = .*/temperature = 20/" "$SRC/temp_resistor.ofe" > "$WORK/tts_ref/x.ofe"
for c in lin tab ref; do tt_run "tts_$c"; done
mkdir -p "$WORK/tts_line"
sed "s/^tempco = .*/tempcotable = 2 20 1 220 1.786/" "$SRC/temp_resistor.ofe" > "$WORK/tts_line/x.ofe"
tt_run tts_line
compare "tempcotable (a straight line) == tempco, R(85 degC) [ohm]" "$(tt_val tts_line R)" "$(tt_val tts_lin R)" 1e-12
compare "tempcotable : R(85 degC) = 1.76 R(20 degC) [ohm]" "$(tt_val tts_tab R)" \
	"$(awk -v r="$(tt_val tts_ref R)" 'BEGIN { printf "%.12e", 1.76 * r }')" 1e-7
#   (d4) 非線形 (B-H) の M で、同居する**線形**磁性体の μr にも表が効くこと。非線形の
#       経路は要素ごとの ν を自前で作るので、そこで 1 次係数を直接書くと表が黙って
#       無視される (実際にあった)。板間を非線形層と線形層 (μr0 = 4、100 degC で 0.5 倍)
#       に分け、μr = 2 を直接書いたものと ofe.out がバイト一致することを見る
for c in tab ref; do mkdir -p "$WORK/ttnl_$c"; done
awk '/^geometry/{print "geometry = 2 1 0 1e-3 0.05e-3 0.15e-3 0 1e-4"
		print "material = 1.0 0"; print "mur = 3 4.0"; print "mutemptable = 3 20 1 100 0.5"
		print "geometry = 3 1 0 1e-3 0.15e-3 0.25e-3 0 1e-4"; next}
	/^analysis/{print "temperature = 100"} {print}' "$SRC/plate_line_bh.ofe" > "$WORK/ttnl_tab/x.ofe"
awk '/^geometry/{print "geometry = 2 1 0 1e-3 0.05e-3 0.15e-3 0 1e-4"
		print "material = 1.0 0"; print "mur = 3 2.0"
		print "geometry = 3 1 0 1e-3 0.15e-3 0.25e-3 0 1e-4"; next}
	/^analysis/{print "temperature = 100"} {print}' "$SRC/plate_line_bh.ofe" > "$WORK/ttnl_ref/x.ofe"
mkdir -p "$WORK/ttnl_lin"
awk '/^mutemptable/{print "mutempco = 3 -5e-3 0"; next} {print}' "$WORK/ttnl_tab/x.ofe" > "$WORK/ttnl_lin/x.ofe"
tt_run ttnl_tab; tt_run ttnl_ref; tt_run ttnl_lin
tt_same "mutemptable on a linear layer inside a nonlinear (B-H) M problem == mur written directly" ttnl_tab ttnl_ref
#       1 次係数 (mutempco) も同じ。**Newton の反復の中で線形材料の ν を作る関数が温度を
#       見ておらず、1 次係数でも黙って無視されていた** (この検査を足して見つかった)
tt_same "mutempco on a linear layer inside a nonlinear (B-H) M problem == mur written directly" ttnl_lin ttnl_ref

#   (e) 弾く入力 (理由まで見る)
tt_reject() {	# tt_reject <ラベル> <sed の置換> <理由>
	sed "$2" "$SRC/temp_mur.ofe" > "$WORK/ttx.ofe"
	mesh_reject "$1" ttx.ofe
	if ! (cd "$WORK" && "$OFE" -n 2 ttx.ofe 2>&1 | grep -q "$3"); then
		echo "  *** it was rejected for a different reason than: $3" >&2
		status=1
	fi
}
tt_reject "mutemptable with temperatures not increasing" "s/^mutempco = .*/mutemptable = 2 20 1 20 0.5/" "invalid mutemptable data"
tt_reject "mutemptable with a non-positive factor" "s/^mutempco = .*/mutemptable = 2 20 1 120 0/" "invalid mutemptable data"
tt_reject "mutemptable with an odd number of values" "s/^mutempco = .*/mutemptable = 2 20 1 120/" "invalid mutemptable data"
tt_reject "mutempco and mutemptable together" "s/^mutempco = .*/mutempco = 2 -3e-3 20\nmutemptable = 2 20 1 120 0.7/" "both mutempco and mutemptable"
sed "s/^tempco = .*/tempco = 2 3.93e-3 20\ntempcotable = 2 20 1 100 1.3/" "$SRC/temp_resistor.ofe" > "$WORK/ttx.ofe"
mesh_reject "tempco and tempcotable together" ttx.ofe
if ! (cd "$WORK" && "$OFE" -n 2 ttx.ofe 2>&1 | grep -q "both tempco and tempcotable"); then
	echo "  *** it was rejected for a different reason than: both tempco and tempcotable" >&2
	status=1
fi

# 非線形 (B-H) と渦電流 (F) の同時解析 (実効磁気抵抗率法)。
echo "[plate_line_bh_ac] nonlinear B-H with eddy currents (effective reluctivity)"
nf_run() {	# nf_run <dir> : $WORK/<dir>/x.ofe を解いて ofe_post まで。非 0 なら NG
	if ! (cd "$WORK/$1" && "$OFE" -n 2 x.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1); then
		echo "  $1 : the run failed -> NG" >&2; status=1
	fi
}
nf_val() { awk -F, -v name="$2" '$1 == name && NF == 2 { f = 1; next } f && $1 == "1" { print $2; exit }' "$WORK/$1/rlc.csv"; }
nf_iters() { awk '/nl-it/ { f = 1; next } f && ($1 ~ /^[0-9]+$/) && (NF == 4) { n = $1 } END { print n + 0 }' "$WORK/$1/ofe.log"; }
#   (a) 1 次元の板間 (data/sample/plate_line_bh_ac.ofe)。板間の H = I/W は周波数にも
#       μ にも依らないので、空気の板間との差は磁性体の磁束だけ:
#         R_nl = R_air、L_nl = L_air + (B(I/W) - mu0 I/W) d / I (どの周波数でも厳密)
#       H が固定されているので H を独立変数にした更新は 2 回で収束する
#       (**B を独立変数にした逐次代入は 2 値を往復して収束しない**、実測)
for f in 1e3 1e5; do
	mkdir -p "$WORK/nfair$f"
	awk -v f="$f" '/^material|^bh|^geometry|^current|^nlsolver/{next} /^frequency/{print "frequency = " f; next} {print}' \
		"$SRC/plate_line_bh_ac.ofe" > "$WORK/nfair$f/x.ofe"
	nf_run "nfair$f"
	# 1e3 Hz では B-H 曲線の最初の点の手前 (0.05 A、初期透磁率の直線) と最後の点の先
	# (20 A、最終区間の傾きで外挿) も通す (逆関数 B(H) の 3 つの枝を全部見る)
	pairs="1.0 1.0|0.5 0.72222222222222|5.0 1.0888888888889"
	[ "$f" = 1e3 ] && pairs="$pairs|0.05 0.25|20.0 1.4222222222222"
	old_ifs=$IFS; IFS='|'; set -- $pairs; IFS=$old_ifs
	for pair in "$@"; do
		set -- $pair
		d="nf$f-$1"
		mkdir -p "$WORK/$d"
		awk -v f="$f" -v cur="$1" '/^frequency/{print "frequency = " f; next} /^current/{print "current = " cur; next} {print}' \
			"$SRC/plate_line_bh_ac.ofe" > "$WORK/$d/x.ofe"
		nf_run "$d"
		compare "f = $f, I = $1 : L = L_air + (B(I/W) - mu0 I/W) d/I [H/m]" "$(nf_val $d Lf)" \
			"$(awk -v la="$(nf_val nfair$f Lf)" -v I="$1" -v B="$2" 'BEGIN { pi = 3.14159265358979324; mu0 = 4 * pi * 1e-7
				printf "%.12e", la + (B - mu0 * I / 1e-3) * 0.2e-3 / I }')" 1e-7
		compare "f = $f, I = $1 : R == R_air (the core carries no current) [ohm/m]" "$(nf_val $d Rf)" "$(nf_val nfair$f Rf)" 1e-8
		it=$(nf_iters "$d")
		if [ "$it" -ge 1 ] && [ "$it" -le 3 ]; then echo "  f = $f, I = $1 : converged in $it nonlinear iterations (<= 3) -> OK"
		else echo "  f = $f, I = $1 : $it nonlinear iterations (expected <= 3) -> NG" >&2; status=1; fi
	done
done
#   (b) 直線の B-H 曲線 (初期透磁率と同じ傾きの 2 点) は線形の μr = 3979 と一致する
mkdir -p "$WORK/nflin" "$WORK/nfmur"
awk '/^bh = 2 10000 /{next} /^bh = 2 1000 /{print "bh = 2 1000 5.0"; next} {print}' "$SRC/plate_line_bh_ac.ofe" > "$WORK/nflin/x.ofe"
awk '/^bh = 2 100 /{print "mur = 2 3978.873577297384"; next} /^bh/{next} /^current|^nlsolver/{next} {print}' "$SRC/plate_line_bh_ac.ofe" > "$WORK/nfmur/x.ofe"
nf_run nflin; nf_run nfmur
compare "a straight-line B-H curve == linear mur, L(f) [H/m]" "$(nf_val nflin Lf)" "$(nf_val nfmur Lf)" 1e-9
#   R は μr ~ 4000 の板間で系が悪条件になり 1e-9 程度で頭打ち (COCG の判定 1e-11)
compare "a straight-line B-H curve == linear mur, R(f) [ohm/m]" "$(nf_val nflin Rf)" "$(nf_val nfmur Rf)" 1e-8
#   (c) 場の出力は与えた電流での物理的な場 : 板間の |B| = B(I/W) = 1.0 T (I = 1 A)
mkdir -p "$WORK/nfvtk"
awk '/^analysis/{print "fieldout = 1"} {print}' "$SRC/plate_line_bh_ac.ofe" > "$WORK/nfvtk/x.ofe"
nf_run nfvtk
bcore=$(awk '/^VECTORS B_F_re_port1 / { f = 1; g = 0; c = 0; next } /^VECTORS B_F_im_port1 / { g = 1; f = 0; c = 0; next }
	f && (NF == 3) && ($1 ~ /^-?[0-9]/) { br[c++] = ($1 * $1) + ($2 * $2) + ($3 * $3); next } f { f = 0 }
	g && (NF == 3) && ($1 ~ /^-?[0-9]/) { bi[c++] = ($1 * $1) + ($2 * $2) + ($3 * $3); next } g { g = 0 }
	END { for (k in br) { b = sqrt(br[k] + bi[k]); if (b > m) m = b } printf "%.9e", m }' "$WORK/nfvtk/ofe_field.vtk")
compare "field output : max |B| in the core == B(I/W) = 1.0 T" "$bcore" 1.0 1e-6
#       電流密度 J_z = σ (s V' - jω Az) も与えた電流での値 : 導体 1 のセルの J_z x 面積の和が
#       current (実部 1 A、虚部 0)。V' に s を掛け忘れると実部がずれる
jsum=$(awk '/^SCALARS conductor / { sc = 1; c = 0; next } sc && /^LOOKUP/ { next }
	sc && (NF == 1) && ($1 ~ /^-?[0-9]/) { cd[c++] = $1; next } sc { sc = 0 }
	/^VECTORS J_F_re_port1 / { fr = 1; fi = 0; c = 0; next } /^VECTORS J_F_im_port1 / { fi = 1; fr = 0; c = 0; next }
	fr && (NF == 3) && ($1 ~ /^-?[0-9]/) { jr[c++] = $3; next } fr { fr = 0 }
	fi && (NF == 3) && ($1 ~ /^-?[0-9]/) { ji[c++] = $3; next } fi { fi = 0 }
	END { a = 0.25e-3 * 2.5e-6; for (k in cd) if (cd[k] == 1) { sr += jr[k] * a; si += ji[k] * a }
	      printf "%.9e %.9e", sr, si }' "$WORK/nfvtk/ofe_field.vtk")
set -- $jsum
compare "field output : sum of J_z over conductor 1 == current (real part) [A]" "$1" 1.0 1e-6
res=$(awk -v v="$2" 'BEGIN { a = (v < 0) ? -v : v; printf "%s (%.1e)", (a < 1e-9) ? "OK" : "NG", v }')
echo "  field output : sum of J_z over conductor 1, imaginary part == 0 : $res"
case "$res" in NG*) status=1 ;; esac
#   (d) 2 次元 (板が磁性体より狭く、縁で磁界が広がる)。F (100 Hz、渦電流はほぼ効かない)
#       は M (Newton) の非線形解と一致する。線形域 (0.2 A) は 1e-5 で一致。飽和域 (5 A)
#       では ν の評価点が違う (F は要素中心、M は Gauss 点) ので h² の差が残り、格子を
#       半分にすると差が 1/4 になる (実測 3.89)。低い周波数ほど F の系は悪条件になる
#       (Az の定数分が 1/ω で大きくなる) ので 1 Hz では比べない (実測 0.18% ずれる)
nf_2d() {	# nf_2d <dir> <F|M> <current> <xdiv> <ydiv-per-region>
	mkdir -p "$WORK/$1"
	awk -v an="$2" -v cur="$3" -v nx="$4" -v ny="$5" '
		/^xmesh/{print "xmesh = 0 " nx " 1e-3"; next}
		/^ymesh/{print "ymesh = 0 " ny " 0.05e-3 " ny " 0.25e-3 " ny " 0.30e-3"; next}
		/^conductor = 0/{print "conductor = 0 1 0.25e-3 0.75e-3 0 0.05e-3 0 1e-4"; next}
		/^conductor = 1/{print "conductor = 1 1 0.25e-3 0.75e-3 0.25e-3 0.30e-3 0 1e-4"; next}
		/^frequency/{ if (an == "F") print "frequency = 100"; next}
		/^current/{print "current = " cur; next}
		/^nlsolver/{print "nlsolver = 100 1e-10 1.0"; next}
		/^analysis/{print "analysis = " an; next} {print}' "$SRC/plate_line_bh_ac.ofe" > "$WORK/$1/x.ofe"
	nf_run "$1"
}
nf_2d nf2F02 F 0.2 8 20; nf_2d nf2M02 M 0.2 8 20
compare "2-D, 0.2 A (initial permeability) : F(100 Hz) == M [H/m]" "$(nf_val nf2F02 Lf)" "$(nf_val nf2M02 Ldc)" 1e-5
nf_2d nf2F5 F 5 8 20; nf_2d nf2M5 M 5 8 20
nf_2d nf2F5h F 5 16 40; nf_2d nf2M5h M 5 16 40
res=$(awk -v a="$(nf_val nf2F5 Lf)" -v b="$(nf_val nf2M5 Ldc)" -v c="$(nf_val nf2F5h Lf)" -v d="$(nf_val nf2M5h Ldc)" 'BEGIN {
	e1 = (a - b) / b; if (e1 < 0) e1 = -e1
	e2 = (c - d) / d; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (diff %.2e -> %.2e, ratio %.2f)", ((e1 < 1e-3) && (q > 3.5) && (q < 4.5)) ? "OK" : "NG", e1, e2, q }')
echo "  2-D, 5 A (saturated) : F(100 Hz) - M shrinks as h^2 : $res"
case "$res" in NG*) status=1 ;; esac
#      Anderson 加速 : 2 次元の飽和域でも 40 回以内 (加速なしの固定点は 0.2 A で 77 回、
#      1 A / 5 A は 100 回で収束しなかった、実測)
it=$(nf_iters nf2F5)
if [ "$it" -ge 1 ] && [ "$it" -le 40 ]; then echo "  2-D, 5 A : converged in $it nonlinear iterations (<= 40, Anderson) -> OK"
else echo "  2-D, 5 A : $it nonlinear iterations (expected <= 40) -> NG" >&2; status=1; fi
#   (e) 弾く入力
nf_reject() {	# nf_reject <ラベル> <awk の置換> <理由>
	awk "$2" "$SRC/plate_line_bh_ac.ofe" > "$WORK/nfx.ofe"
	mesh_reject "$1" nfx.ofe
	if ! (cd "$WORK" && "$OFE" -n 2 nfx.ofe 2>&1 | grep -q "$3"); then
		echo "  *** it was rejected for a different reason than: $3" >&2; status=1
	fi
}
nf_reject "nonlinear F with two ports" '/^conductor = 1/{print; print "conductor = 2 1 0.4e-3 0.6e-3 0.12e-3 0.18e-3 0 1e-4"; print "conductorsigma = 2 5.8e7"; next} {print}' "exactly one port"
nf_reject "nonlinear F with a frequency sweep" '/^frequency/{print "frequencysweep = 1e3 1e4"; next} {print}' "cannot be combined with frequencysweep"
nf_reject "nonlinear F with per-axis B-H curves" '/^bh = 2 100 /{print "bh = 2 x 100 0.5"; print "bh = 2 y 100 0.5"; print "bh = 2 z 100 0.5"; next} /^bh/{next} {print}' "isotropic bh curves only"

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

# 非線形 (B-H) と 3 次元渦電流 (A) : bar_air の層を B-H 曲線の磁性体に替える
# (data/sample/bar_air_bh.ofe)。層の H が一様なので、界面の Robin 条件の κ = g μr を
# 自己無撞着に解いた 1 次元の閉形式 (data/sample/bar_layer.awk) と比べる
echo "[bar_air_bh] nonlinear B-H with 3-D eddy currents (effective reluctivity)"
na_val() { awk -F, -v name="$2" '$1 == name && NF == 2 { f = 1; next } f && $1 == "1" { print $2; exit }' "$WORK/$1/rlc.csv"; }
na_run() {	# na_run <dir> : $WORK/<dir>/x.ofe を解いて ofe_post まで
	cp "$SRC/bar_air.msh" "$WORK/$1/"
	if ! (cd "$WORK/$1" && "$OFE" -n 2 x.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1); then
		echo "  $1 : the run failed -> NG" >&2; status=1
	fi
}
#   (a) 3 つの電流で曲線の 2 区間目 (150 A)・3 区間目 (200 A)・外挿の区間 (250 A) を通す。
#       許容 0.1% : 線形の層 (μr = 50) と同じ格子の離散化誤差が L -0.028% (実測)
#       磁束で駆動される構成 (層の μ を上げると H が下がる) なので収束は遅く、
#       上限 100 回 (実測 16〜38 回、許容 1e-7)。上限を越えたら入力エラーになる。
#       許容を内側の COCG と同じ 1e-9 にすると回数が丸めの雑音で決まり、CI で上限を
#       越えた (手元 54 回)
for I in 150 200 250; do
	mkdir -p "$WORK/na$I"
	sed "s/^current = .*/current = $I/" "$SRC/bar_air_bh.ofe" > "$WORK/na$I/x.ofe"
	na_run "na$I"
	set -- $(awk -v f=1e2 -v mode=nl -v I="$I" -f "$SRC/bar_layer.awk")
	compare "I = $I A (layer H = $(printf '%.0f' "$3") A/m, mur = $(printf '%.0f' "$4")) : R [ohm]" "$(na_val na$I Rf)" "$1" 1e-3
	compare "I = $I A : L [H]" "$(na_val na$I Lf)" "$2" 1e-3
done
#   非線形は電流振幅で決まるので、励振の電圧 (voltage) に依らない。既定の 1 V だと
#   端子電流に電圧を掛け忘れる誤りが等価になるので、0.5 V でも解いて一致を見る
mkdir -p "$WORK/na200v"
awk '/^frequency/{print; print "voltage = 0.5"; next} {print}' "$WORK/na200/x.ofe" > "$WORK/na200v/x.ofe"
na_run na200v
compare "nonlinear A does not depend on the drive voltage (0.5 V vs 1 V) : L [H]" "$(na_val na200v Lf)" "$(na_val na200 Lf)" 1e-7
#   線形の層の閉形式 (κ = g μr) そのものも見ておく (μr = 50、同じ許容)
mkdir -p "$WORK/nal50"
awk '/^bh/{next} /^current|^nlsolver/{next} /^region = 2 3/{print; print "mur = 3 50"; next} {print}' \
	"$SRC/bar_air_bh.ofe" > "$WORK/nal50/x.ofe"
na_run nal50
set -- $(awk -v f=1e2 -v mode=lin -v mr=50 -v I=1 -f "$SRC/bar_layer.awk")
compare "linear layer, mur = 50 (kappa = g mur) : L [H]" "$(na_val nal50 Lf)" "$2" 1e-3
#   (b) 直線の B-H 曲線 (初期透磁率と同じ傾き) は線形の μr = 3979 の層と一致する
mkdir -p "$WORK/nalin" "$WORK/namur"
awk '/^bh = 3 10000 /{next} /^bh = 3 1000 /{print "bh = 3 1000 5.0"; next} {print}' "$SRC/bar_air_bh.ofe" > "$WORK/nalin/x.ofe"
awk '/^bh/{next} /^current|^nlsolver/{next} /^region = 2 3/{print; print "mur = 3 3978.873577297384"; next} {print}' \
	"$SRC/bar_air_bh.ofe" > "$WORK/namur/x.ofe"
na_run nalin; na_run namur
compare "a straight-line B-H layer == linear mur, L [H]" "$(na_val nalin Lf)" "$(na_val namur Lf)" 1e-7
compare "a straight-line B-H layer == linear mur, R [ohm]" "$(na_val nalin Rf)" "$(na_val namur Rf)" 1e-7
#   (c) 弾く入力。**A の検査は非構造格子の分岐の中に置く** (構造格子向けの bh の検査は
#       非構造格子の早期 return の後ろにあって A には効かない)
na_reject() {	# na_reject <ラベル> <awk の置換> <理由>
	awk "$2" "$SRC/bar_air_bh.ofe" > "$WORK/nax.ofe"
	mesh_reject "$1" nax.ofe
	if ! (cd "$WORK" && "$OFE" -n 2 nax.ofe 2>&1 | grep -q "$3"); then
		echo "  *** it was rejected for a different reason than: $3" >&2; status=1
	fi
}
na_reject "nonlinear A with a frequency sweep" '/^frequency/{print "frequencysweep = 1e2 1e3"; next} {print}' "cannot be combined with frequencysweep"
na_reject "nonlinear A with per-axis B-H curves" '/^bh = 3 100 /{print "bh = 3 x 100 0.5"; print "bh = 3 y 100 0.5"; print "bh = 3 z 100 0.5"; next} /^bh/{next} {print}' "isotropic bh curves only"

# ---- 電気-熱連成 (thermal) ----
#
# 電流が一様な棒では熱の式が T について線形になり、閉形式が残る
# (data/sample/thermal_bar.ofe のコメント参照):
#   R = (ρ_s L / A) tan(mL/2) / (mL/2),  T_max = T_s + (1/α')(1/cos(mL/2) - 1)
# α < 0 では tanh / cosh、α = 0 では放物線。mL = π で定常解が無くなる (熱暴走)。
echo "[thermal] electro-thermal coupling : a bar cooled at both ends (closed form)"
th_r() { awk '/^  R = / { print $3; exit }' "$WORK/ofe.log"; }
th_t() { awk '/^  T_max = / { print $3; exit }' "$WORK/ofe.log"; }
# th_exact <I> <alpha> [one] -> "R T_max" (thermal_bar.ofe の形状・材料)
#   one を付けると片側 (x = 0) だけヒートシンクで反対側は断熱。長さ 2L の両側
#   ヒートシンクの棒の半分と同じなので、式の mL/2 が mL になる (熱暴走は mL = π/2)
th_exact() {
	awk -v I="$1" -v a="$2" -v one="$3" 'BEGIN {
		r0 = 1e-6; k = 10; L = 0.01; A = 1e-6; T0 = 0; Ts = 20
		rs = r0 * (1 + a * (Ts - T0))
		J = I / A
		h = ((one == "one") ? L : L / 2)
		if (a == 0) { printf "%.12e %.12e", rs * L / A, Ts + J * J * rs * h * h / (2 * k); exit }
		ap = a / (1 + a * (Ts - T0))
		m2 = J * J * r0 * a / k
		if (m2 > 0) { x = sqrt(m2) * h; t = sin(x) / cos(x); c = cos(x) }
		else { x = sqrt(-m2) * h; ex = exp(x); t = (ex - 1 / ex) / (ex + 1 / ex); c = (ex + 1 / ex) / 2 }
		printf "%.12e %.12e", rs * L / A * t / x, Ts + (1 / ap) * (1 / c - 1) }'
}
# 発熱の恒等式 : I^2 R (端子) = Σ f (場の積分) = ヒートシンク + 対流 + 輻射で出ていく熱
th_energy() {	# th_energy <ラベル>
	res=$(awk '/joule power/ {
		for (i = 1; i <= NF; i++) {
			if ($i == "R") a = $(i+2)
			if ($i == "integral") b = $(i+2)
			if ($i == "sinks") c = $(i+2) + 0
			if ($i == "convection") e = $(i+2) + 0
			if ($i == "radiation") g = $(i+2) + 0
		}
	}
	END { if ((a == "") || (e == "") || (g == "")) { printf "NG (no joule line)"; exit }
	      d1 = (b - a) / a; if (d1 < 0) d1 = -d1
	      d2 = (c + e + g - a) / a; if (d2 < 0) d2 = -d2
	      printf "%s (field %.1e, out %.1e)", ((d1 < 1e-9) && (d2 < 1e-9)) ? "OK" : "NG", d1, d2 }' "$WORK/ofe.log")
	echo "  $1 : I^2 R = field integral = heat out (sinks + convection + radiation) : $res"
	case "$res" in NG*) status=1 ;; esac
}
th_run() {	# th_run <ofe> : 実行して非 0 なら NG
	if ! (cd "$WORK" && "$OFE" -n 2 "$1" > /dev/null 2>&1); then
		echo "  $1 : the run failed -> NG" >&2
		status=1
	fi
}
cp "$SRC/thermal_bar.ofe" "$SRC/thermal_neck.ofe" "$WORK/"
# (1) 閉形式 (I = 7.5 A, mL = 1.5)。許容 0.1% の根拠 : 要素内で σ を一定にするので
#     離散化誤差は h² で、40 分割で R -0.018% / T_max -0.025% (実測)。
#     次の (2) で誤差の次数そのものを見るので、許容を緩めても退行は隠れない
th_run thermal_bar.ofe
set -- $(th_exact 7.5 4e-3)
compare "R(I = 7.5 A, mL = 1.5) [ohm]" "$(th_r)" "$1" 0.001
compare "T_max(I = 7.5 A) [degC]" "$(th_t)" "$2" 0.001
th_energy "thermal_bar"
r40=$(th_r); rex=$1
# (2) 誤差の次数 : 20 分割と 40 分割の誤差比が 4 (h²)。要素の温度を 8 節点の平均で
#     なく 1 節点で取るような誤りは h の 1 次になって比が 2 に落ちる
awk '/^xmesh/{print "xmesh = 0 20 0.01"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thb20.ofe"
th_run thb20.ofe
r20=$(th_r)
res=$(awk -v a="$r20" -v b="$r40" -v e="$rex" 'BEGIN {
	e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.2)) ? "OK" : "NG", e1, e2, q }')
echo "  second-order convergence of R (20 -> 40 cells) : $res"
case "$res" in NG*) status=1 ;; esac
# (2b) 片側だけヒートシンク (反対側は断熱)。**両側の棒は左右対称なので、要素の
#     温度を 8 節点の平均でなく 1 節点で取る誤り (評価点が h/2 ずれる) が
#     左半分と右半分で逆符号の 1 次誤差になって打ち消し、h² のまま通ってしまう**
#     (実測: 両側の棒では誤差比 4.0 のまま)。温度が単調な片側の棒では打ち消さず、
#     1 次に落ちる (実測: 誤差比 2.0、40 分割で -0.44%)
th_one() {	# th_one <nx> : 片側ヒートシンク、I = 3.75 A (mL = 0.75 < π/2)
	awk -v nx="$1" '/^xmesh/{print "xmesh = 0 " nx " 0.01"; next}
		/^heatsink = 1/{next} /^thermal/{print "thermal = 3.75"; next} {print}' \
		"$SRC/thermal_bar.ofe" > "$WORK/thone$1.ofe"
	th_run "thone$1.ofe"
}
th_one 40
set -- $(th_exact 3.75 4e-3 one)
compare "R(one-sided heat sink, mL = 0.75) [ohm]" "$(th_r)" "$1" 0.001
compare "T_max(one-sided heat sink) [degC]" "$(th_t)" "$2" 0.001
r40=$(th_r); rex=$1
th_one 20
r20=$(th_r)
res=$(awk -v a="$r20" -v b="$r40" -v e="$rex" 'BEGIN {
	e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.2)) ? "OK" : "NG", e1, e2, q }')
echo "  second-order convergence, one-sided (20 -> 40 cells) : $res"
case "$res" in NG*) status=1 ;; esac
# (3) α = 0 : σ が温度に依らないので R は冷えた棒の値に厳密に一致し、一様断面の
#     棒の 1 次要素解は節点で厳密 (定数の発熱源の 1 次元 Poisson) なので
#     T_max も放物線の頂点に全桁一致する
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/tha0.ofe"
th_run tha0.ofe
set -- $(th_exact 7.5 0)
compare "R(alpha = 0) == the uncoupled value [ohm]" "$(th_r)" "$1" 1e-10
compare "T_max(alpha = 0) == parabola, nodally exact [degC]" "$(th_t)" "$2" 1e-10
# (4) α < 0 (NTC) : tanh / cosh の枝
awk '/^tempco/{print "tempco = 2 -4e-3 0"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thntc.ofe"
th_run thntc.ofe
set -- $(th_exact 7.5 -4e-3)
compare "R(alpha < 0, NTC) [ohm]" "$(th_r)" "$1" 0.001
compare "T_max(alpha < 0, NTC) [degC]" "$(th_t)" "$2" 0.001
# (5) 熱暴走の境目 (臨界電流 15.708 A, mL = π) の両側。
#     手前 (mL = 2.9) は収束して閉形式に合う。臨界に近いほど解が敏感になる
#     (tan の極) ので許容は 1% (実測 -0.53%)。Picard 反復の収縮率はほぼ
#     (mL/π)² なので、ここでは 130 回程度かかる
awk '/^thermal/{print "thermal = 14.5"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thnear.ofe"
th_run thnear.ofe
set -- $(th_exact 14.5 4e-3)
compare "R(I = 14.5 A, mL = 2.9, near runaway) [ohm]" "$(th_r)" "$1" 0.01
#     向こう側 (mL = 3.3) は定常解が無いので、黙って値を出さずに弾くこと
awk '/^thermal/{print "thermal = 16.5"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thrun.ofe"
mesh_reject "I = 16.5 A (mL = 3.3 > pi, no steady state)" thrun.ofe
if ! grep -q "thermal runaway" "$WORK/ofe.log"; then
	echo "  *** it was rejected for a different reason than thermal runaway" >&2
	status=1
fi
# (6) 首のある導体 : ∇φ が要素の中で変わるので、発熱の積分則を検査できる
#     (1 次元の棒では要素中心の 1 点則でも厳密に通ってしまう)
th_run thermal_neck.ofe
th_energy "thermal_neck (current crowding)"
#     場の出力の温度の最大 = ログの T_max (VTK は %.9e、ログは %.10e)
tvtk=$(awk '/^SCALARS T / { f = 1; next } f && /^LOOKUP/ { next }
	f && (NF == 1) { if (($1 + 0) > m) m = $1 + 0; next } f { f = 0 }
	END { printf "%.9e", m }' "$WORK/ofe_field.vtk")
compare "max of the T field == T_max in the log [degC]" "$tvtk" "$(th_t)" 1e-8
# (8) 対流 (convection)。data/sample/thermal_conv.ofe : 側面だけから対流で冷やす棒。
#     両端が断熱で発熱も冷却も一様なので温度は一様、断面を 1 x 1 セルにしてあるので
#     離散解も厳密に一様になる (どの節点でも受け持つ体積と側面積の比が A/P)。
#     許容 1e-9 は Picard 反復の打ち切り (1e-10、収縮率 (I/I_c)² = 0.48) の分
echo "[thermal] convection : a bar cooled from its sides only (uniform, exact)"
th_conv() {	# th_conv <I> -> "R T_max" (thermal_conv.ofe の形状・材料)
	awk -v I="$1" 'BEGIN {
		r0 = 1e-6; a = 4e-3; T0 = 0; Ta = 30; h = 200; L = 0.01; A = 5e-7; P = 3e-3
		J = I / A; ra = r0 * (1 + a * (Ta - T0)); aa = a / (1 + a * (Ta - T0))
		th = J * J * ra * A / (h * P - J * J * r0 * a * A)
		printf "%.12e %.12e", ra * (1 + aa * th) * L / A, Ta + th }'
}
th_faces() { awk '/exposed faces/ { for (i = 1; i < NF; i++) if ($(i+1) == "exposed") print $i; exit }' "$WORK/ofe.log"; }
th_nface() {	# th_nface <ラベル> <期待値> : 対流を課した面の数
	if [ "$(th_faces)" = "$2" ]; then echo "  $1 : $2 exposed faces -> OK"
	else echo "  $1 : $(th_faces) exposed faces, expected $2 -> NG" >&2; status=1; fi
}
cp "$SRC/thermal_conv.ofe" "$WORK/"
th_run thermal_conv.ofe
set -- $(th_conv 6)
compare "R(convection only, I = 6 A) [ohm]" "$(th_r)" "$1" 1e-9
compare "T_max(convection only, uniform) [degC]" "$(th_t)" "$2" 1e-9
th_energy "thermal_conv"
#     電極の面 (棒の両端) は除くので側面の 40 x 4 面だけ。両端に課すと温度が
#     一様でなくなるうえ数が 162 になる
th_nface "thermal_conv (sides only, not the electrode faces)" 160
#     棒の向きを y / z に回す (座標の巡回置換)。x 向きの棒では x 方向の面が
#     電極の面として除かれるだけで、面の質量の組み立てを一度も通らない
th_rot() {	# th_rot <y|z> : thermal_conv.ofe の座標を巡回置換した棒
	awk -v ax="$1" '
		function perm(a, b, c) { if (ax == "y") return c " " a " " b; return b " " c " " a }
		/^xmesh/ { X = $3 " " $4 " " $5; next }
		/^ymesh/ { Y = $3 " " $4 " " $5; next }
		/^zmesh/ { Z = $3 " " $4 " " $5
			if (ax == "y") { print "xmesh = " Z; print "ymesh = " X; print "zmesh = " Y }
			else { print "xmesh = " Y; print "ymesh = " Z; print "zmesh = " X }
			next }
		/^geometry|^conductor/ { print $1, $2, $3, $4, perm($5 " " $6, $7 " " $8, $9 " " $10); next }
		{ print }' "$SRC/thermal_conv.ofe" > "$WORK/thconv_$1.ofe"
	th_run "thconv_$1.ofe"
}
for ax in y z; do
	th_rot $ax
	set -- $(th_conv 6)
	compare "R(convection only, bar along $ax) [ohm]" "$(th_r)" "$1" 1e-9
	th_nface "bar along $ax" 160
done
#     熱を通さないセル (空気) に囲まれた棒。露出面は格子の外だけでなく「隣が
#     熱を通さないセル」の面も含む。格子の外だけを見る誤りはこれでだけ落ちる
awk '/^ymesh/{print "ymesh = -1e-3 1 0 1 1e-3 1 2e-3"; next}
	/^zmesh/{print "zmesh = -1e-3 1 0 1 0.5e-3 1 1.5e-3"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thconv_air.ofe"
th_run thconv_air.ofe
set -- $(th_conv 6)
compare "R(convection only, bar embedded in air) [ohm]" "$(th_r)" "$1" 1e-9
th_nface "bar embedded in air" 160
#     熱暴走の境目 (I_c = 8.660 A) の両側。手前 (I = 8、収縮率 0.85) は
#     130 回程度で収束して厳密値に合い、向こう側 (I = 9) は弾く
awk '/^thermal/{print "thermal = 8"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thconv8.ofe"
th_run thconv8.ofe
set -- $(th_conv 8)
compare "R(convection only, I = 8 A, near runaway) [ohm]" "$(th_r)" "$1" 1e-8
awk '/^thermal/{print "thermal = 9"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thconv9.ofe"
mesh_reject "convection only, I = 9 A (> I_c = 8.66 A, no steady state)" thconv9.ofe
if ! grep -q "thermal runaway" "$WORK/ofe.log"; then
	echo "  *** it was rejected for a different reason than thermal runaway" >&2
	status=1
fi
# (9) ヒートシンク + 側面の対流 (フィン)。θ = T - T_s として
#       θ'' + μ² θ = -c,  μ² = (J² ρ0 α A - h P) / (k A),  c = (J² ρ_s A - h P (T_s - T_a)) / (k A)
#       R = (ρ_s L + ρ0 α ∫θ) / A,  ∫θ = (c/μ²) ((2/μ) tan(μL/2) - L)   (μ² < 0 なら tanh)
#     h = 20 で μ² > 0 (tan)、h = 200 で μ² < 0 (tanh)。T_a = 30 ≠ T_s = 20 なので
#     対流が熱を持ち込む項 (c の第 2 項) も効く。ヒートシンクの節点は側面にも
#     接しているので、「ヒートシンクへ出ていく熱」に対流の右辺を入れ忘れると
#     エネルギーの恒等式が落ちる
th_fin() {	# th_fin <I> <h> -> "R T_max" (thermal_bar.ofe + convection = h 30)
	awk -v I="$1" -v h="$2" 'BEGIN {
		r0 = 1e-6; k = 10; L = 0.01; A = 1e-6; P = 4e-3; T0 = 0; Ts = 20; Ta = 30; a = 4e-3
		J = I / A; rs = r0 * (1 + a * (Ts - T0))
		mu2 = (J * J * r0 * a * A - h * P) / (k * A)
		c = (J * J * rs * A - h * P * (Ts - Ta)) / (k * A)
		if (mu2 > 0) { x = sqrt(mu2); t = sin(x * L / 2) / cos(x * L / 2); cc = cos(x * L / 2) }
		else { x = sqrt(-mu2); ex = exp(x * L / 2); t = (ex - 1 / ex) / (ex + 1 / ex); cc = (ex + 1 / ex) / 2 }
		i1 = (c / mu2) * ((2 / x) * t - L)
		printf "%.12e %.12e", (rs * L + r0 * a * i1) / A, Ts + (c / mu2) * (1 / cc - 1) }'
}
for h in 20 200; do
	for nx in 20 40; do
		awk -v nx="$nx" -v h="$h" '/^xmesh/{print "xmesh = 0 " nx " 0.01"; next}
			/^thermal/{print "convection = " h " 30"} {print}' \
			"$SRC/thermal_bar.ofe" > "$WORK/thfin$h-$nx.ofe"
		th_run "thfin$h-$nx.ofe"
		eval "r$nx=\$(th_r)"
	done
	set -- $(th_fin 7.5 "$h")
	#   許容 0.1% : h² の離散化誤差 (40 分割で -0.015% / -0.007%、実測)
	compare "R(fin, h = $h) [ohm]" "$r40" "$1" 0.001
	compare "T_max(fin, h = $h) [degC]" "$(th_t)" "$2" 0.001
	th_energy "fin, h = $h"
	res=$(awk -v a="$r20" -v b="$r40" -v e="$1" 'BEGIN {
		e1 = (a - e) / e; if (e1 < 0) e1 = -e1
		e2 = (b - e) / e; if (e2 < 0) e2 = -e2
		q = (e2 > 0) ? e1 / e2 : 0
		printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.2)) ? "OK" : "NG", e1, e2, q }')
	echo "  second-order convergence, fin h = $h (20 -> 40 cells) : $res"
	case "$res" in NG*) status=1 ;; esac
done
# (10) 首のある導体を対流だけで冷やす。x / y 方向の面 (パッドの肩と首の側面) と
#      空気に接する面が全部出てくる 3 次元の形で、発熱 = 対流で出ていく熱
awk '!/^heatsink/ && /^thermal/{print "convection = 500 30"} !/^heatsink/{print}' \
	"$SRC/thermal_neck.ofe" > "$WORK/thneck_conv.ofe"
th_run thneck_conv.ofe
th_energy "thermal_neck cooled by convection only"
#      エネルギーの恒等式は「どの面に課したか」を見ない (内部の面に課しても
#      つじつまは合う) ので、面の数を別に数える:
#        上下 (z) : 導体の柱 (パッド 6x16x2 + 首 6x4) x 2 = 432
#        y 方向   : パッドの両脇 6x2x2x2 = 48 + 首の両脇 6x2x2 = 24
#        x 方向   : パッドの肩 (首の外の 12x2) x 2 = 48 (両端は電極の面で除く)
th_nface "thermal_neck (432 + 72 + 48)" 552

# (11) 輻射 (radiation)。-k dT/dn = ε σ_SB (T_K^4 - T_r,K^4)、絶対温度。
#      側面だけから冷やす棒 (thermal_conv.ofe の形状) は温度が一様なので、釣り合い
#        J² ρ0 (1 + α (T - T0)) A = P [h (T - T_c) + ε σ_SB (T_K^4 - T_r,K^4)]
#      を二分法で解いたものと機械精度で一致する (断面 1 x 1 セルなので離散解も一様)
echo "[thermal] radiation : uniform bars (algebraic), a bar with heat sinks (shooting)"
th_radu() {	# th_radu <I> <h> <T_c> <ε> <T_r> [<α>] -> "R T"
	awk -v I="$1" -v h="$2" -v Tc="$3" -v e="$4" -v Tr="$5" -v a="${6:-4e-3}" 'BEGIN {
		r0 = 1e-6; T0 = 0; L = 0.01; A = 5e-7; P = 3e-3; sb = 5.670374419e-8; K = 273.15
		J = I / A; lo = -200; hi = 1e5
		for (it = 0; it < 300; it++) {
			T = (lo + hi) / 2
			F = J * J * r0 * (1 + a * (T - T0)) * A - P * (h * (T - Tc) + e * sb * ((T + K) ^ 4 - (Tr + K) ^ 4))
			if (F > 0) lo = T; else hi = T
		}
		printf "%.12e %.12e", r0 * (1 + a * (T - T0)) * L / A, T }'
}
#      輻射だけ (ε = 0.8, T_r = 10, I = 2 A : 約 345 degC)。T_r を他のどの温度とも
#      違えてあるのは取り違えを落とすため。摂氏のまま 4 乗する誤りもここで落ちる
awk '/^convection/{print "radiation = 0.8 10"; next} /^thermal/{print "thermal = 2"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thrad_u.ofe"
th_run thrad_u.ofe
set -- $(th_radu 2 0 30 0.8 10)
compare "R(radiation only, uniform) [ohm]" "$(th_r)" "$1" 1e-9
compare "T_max(radiation only, uniform) [degC]" "$(th_t)" "$2" 1e-9
th_energy "radiation only"
th_nface "radiation only (sides, not the electrode faces)" 160
#      対流 + 輻射 (I = 6 A、輻射が熱の約 7%)
awk '/^convection/{print; print "radiation = 0.8 10"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thrad_cu.ofe"
th_run thrad_cu.ofe
set -- $(th_radu 6 200 30 0.8 10)
compare "R(convection + radiation, uniform) [ohm]" "$(th_r)" "$1" 1e-9
th_energy "convection + radiation"
#      α = 0 なら非線形は輻射だけで、Newton 線形化なら 2 次収束して数回で止まる
#      (実測 7 回)。**線形化を誤っても (接線を落とす・割線にする) 同じ不動点に
#      収束し得るので値の比較では捕まらない**。回数そのものを assert する
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$WORK/thrad_u.ofe" > "$WORK/thrad_a0.ofe"
th_run thrad_a0.ofe
set -- $(th_radu 2 0 30 0.8 10 0)
compare "T_max(radiation only, alpha = 0) [degC]" "$(th_t)" "$2" 1e-9
nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
if [ -n "$nit" ] && [ "$nit" -le 8 ]; then
	echo "  Newton on radiation (alpha = 0) : converged in $nit iterations (<= 8) -> OK"
else
	echo "  Newton on radiation (alpha = 0) : ${nit:-?} iterations (> 8, not quadratic) -> NG" >&2
	status=1
fi
#      ヒートシンク + 輻射 (data/sample/thermal_rad.ofe)。1 次元の非線形 2 点境界値問題を
#      中央 (T' = 0) から RK4 で撃ち、端で T_s になる中央の温度を二分法で求める
th_radfin() {	# th_radfin -> "R T_max" (thermal_rad.ofe の形状・材料)
	awk 'function f(T) { return (-J * J * r0 * (1 + a * (T - T0)) + (P / A) * e * sb * ((T + K) ^ 4 - (Tr + K) ^ 4)) / k }
	function rho(T) { return r0 * (1 + a * (T - T0)) }
	function shoot(Tc,   i, hs, T, D, k1, k2, k3, k4, l1, l2, l3, l4) {
		hs = (L / 2) / n; T = Tc; D = 0; S = 0
		for (i = 0; i < n; i++) {
			k1 = D;               l1 = f(T)
			k2 = D + hs * l1 / 2; l2 = f(T + hs * k1 / 2)
			k3 = D + hs * l2 / 2; l3 = f(T + hs * k2 / 2)
			k4 = D + hs * l3;     l4 = f(T + hs * k3)
			S += hs * (rho(T) + 2 * rho(T + hs * k1 / 2) + 2 * rho(T + hs * k2 / 2) + rho(T + hs * k3)) / 6
			T += hs * (k1 + 2 * k2 + 2 * k3 + k4) / 6
			D += hs * (l1 + 2 * l2 + 2 * l3 + l4) / 6
		}
		return T
	}
	BEGIN {
		r0 = 1e-6; a = 4e-3; T0 = 0; L = 0.01; A = 1e-6; P = 4e-3; k = 1; e = 1; Tr = 10
		sb = 5.670374419e-8; K = 273.15; Ts = 20; J = 2.5 / A; n = 4000
		lo = Ts; hi = Ts + 500
		for (it = 0; it < 80; it++) {
			Tc = (lo + hi) / 2
			if (shoot(Tc) > Ts) hi = Tc; else lo = Tc
		}
		shoot(Tc)
		printf "%.12e %.12e", 2 * S / A, Tc }'
}
for nx in 20 40; do
	awk -v nx="$nx" '/^xmesh/{print "xmesh = 0 " nx " 0.01"; next} {print}' \
		"$SRC/thermal_rad.ofe" > "$WORK/thradf$nx.ofe"
	th_run "thradf$nx.ofe"
	eval "r$nx=\$(th_r)"
done
set -- $(th_radfin)
#   許容 0.1% : h² の離散化誤差 (40 分割で -0.011%、実測)
compare "R(heat sinks + radiation, 40 cells) [ohm]" "$r40" "$1" 0.001
compare "T_max(heat sinks + radiation) [degC]" "$(th_t)" "$2" 0.001
th_energy "heat sinks + radiation"
res=$(awk -v a="$r20" -v b="$r40" -v e="$1" 'BEGIN {
	e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.2)) ? "OK" : "NG", e1, e2, q }')
echo "  second-order convergence, radiation (20 -> 40 cells) : $res"
case "$res" in NG*) status=1 ;; esac
#      首のある導体を対流 + 輻射で冷やす (ヒートシンクも残す) : 3 経路の収支
awk '/^thermal/{print "convection = 100 30"; print "radiation = 0.9 10"} {print}' \
	"$SRC/thermal_neck.ofe" > "$WORK/thneck_rad.ofe"
th_run thneck_rad.ofe
th_energy "thermal_neck : sinks + convection + radiation"
th_nface "thermal_neck with convection and radiation" 552

# (12) 過渡解析 (transient = <t_end> <nstep> [<θ>]、heatcap = ρc)。θ 法で
#        M (T^{n+1} - T^n)/Δt + θ R(T^{n+1}) + (1-θ) R(T^n) = 0
#      を各刻みで (σ(T) と輻射を含めて) 完全陰的に解く。
echo "[thermal] transient : uniform bars (exact recurrence), Fourier series, steady limit"
th_balance() {	# th_balance <ラベル> : 発熱 = 蓄熱 + ヒートシンク + 対流 + 輻射 [J]
	res=$(awk '/heat balance/ {
		for (i = 1; i <= NF; i++) {
			if ($i == "joule") a = $(i+2) + 0
			if ($i == "stored") b = $(i+2) + 0
			if ($i == "sinks") c = $(i+2) + 0
			if ($i == "convection") e = $(i+2) + 0
			if ($i == "radiation") g = $(i+2) + 0
		}
		got = 1
	}
	END { if (!got || (a == 0)) { printf "NG (no heat balance line)"; exit }
	      d = (b + c + e + g - a) / a; if (d < 0) d = -d
	      printf "%s (%.1e)", (d < 1e-9) ? "OK" : "NG", d }' "$WORK/ofe.log")
	echo "  $1 : joule = stored + sinks + convection + radiation : $res"
	case "$res" in NG*) status=1 ;; esac
}
#      一様な棒では温度が一様のまま (断面 1 x 1 セル、電極の面は断熱) なので、
#      θ 法の離散解が 1 次の漸化式で**厳密に**書ける。y = θ - θ_∞ が y' = λ y に
#      従うとき y_N = y_0 ((1 + (1-θ) λ Δt) / (1 - θ λ Δt))^N (連続解は y_0 e^{λt})
th_ode() {	# th_ode <λ> <y_0> <T_∞> <t_end> <N> <θ> <ρ0 (1+α(T-T0)) を R にする L/A> <α> <T0>
	awk -v lam="$1" -v y0="$2" -v Tf="$3" -v te="$4" -v ns="$5" -v th="$6" -v la="$7" -v a="$8" -v T0="$9" 'BEGIN {
		dt = te / ns
		yd = y0 * ((1 + (1 - th) * lam * dt) / (1 - th * lam * dt)) ^ ns
		yc = y0 * exp(lam * te)
		printf "%.12e %.12e %.12e %.12e", 1e-6 * (1 + a * (Tf + yd - T0)) * la, Tf + yd,
			1e-6 * (1 + a * (Tf + yc - T0)) * la, Tf + yc }'
}
#      (a) 断熱の棒 (thermal_bar.ofe からヒートシンクを外す)。
#          ρc dθ/dt = J² ρ_i (1 + α_i θ) なので λ = J² ρ0 α / ρc、y = θ + 1/α_i
th_adia() {	# th_adia <N> <θ> -> 離散解の R, T と連続解の R, T
	awk -v I=7.5 'BEGIN { J = I / 1e-6; ai = 4e-3 / (1 + 4e-3 * 50)
		printf "%.15e %.15e %.15e", J * J * 1e-6 * 4e-3 / 3.4e6, 1 / ai, 50 - 1 / ai }' > "$WORK/thode.txt"
	set -- $(cat "$WORK/thode.txt") "$1" "$2"
	th_ode "$1" "$2" "$3" 10 "$4" "$5" 1e4 4e-3 0
}
for th in 1 0.5; do
	for ns in 50 100; do
		awk -v th="$th" -v ns="$ns" '/^heatsink/{next}
			/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 10 " ns " " th} {print}' \
			"$SRC/thermal_bar.ofe" > "$WORK/thad$ns.ofe"
		th_run "thad$ns.ofe"
		eval "r$ns=\$(th_r)"
	done
	set -- $(th_adia 100 "$th")
	#   許容 1e-9 : 内側の反復の打ち切り (1e-10) の分。漸化式は厳密
	compare "adiabatic bar, theta = $th : R == the discrete recurrence [ohm]" "$r100" "$1" 1e-9
	compare "adiabatic bar, theta = $th : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
	th_balance "adiabatic bar, theta = $th"
	#   時間の次数 : 連続解に対する 50 -> 100 刻みの誤差比が 2 (後退 Euler) / 4 (CN)
	rc=$3
	set -- $(th_adia 50 "$th")
	res=$(awk -v a="$r50" -v b="$r100" -v e="$rc" -v th="$th" 'BEGIN {
		e1 = (a - e) / e; if (e1 < 0) e1 = -e1
		e2 = (b - e) / e; if (e2 < 0) e2 = -e2
		q = (e2 > 0) ? e1 / e2 : 0; w = ((th == 1) ? 2 : 4)
		printf "%s (err %.2e -> %.2e, ratio %.3f, expected %d)", ((q > 0.95 * w) && (q < 1.05 * w)) ? "OK" : "NG", e1, e2, q, w }')
	echo "  time order, theta = $th (50 -> 100 steps) : $res"
	case "$res" in NG*) status=1 ;; esac
done
#      csv : 行数 = 刻み + 1 (t = 0 を含む)、最初の行は初期温度の R、最後の行はログの R
nrow=$(awk 'NR > 1' "$WORK/ofe_transient.csv" | wc -l | tr -d ' ')
if [ "$nrow" = 101 ]; then echo "  ofe_transient.csv : 101 rows (100 steps + t = 0) -> OK"
else echo "  ofe_transient.csv : $nrow rows, expected 101 -> NG" >&2; status=1; fi
compare "ofe_transient.csv : R at t = 0 == rho(50 degC) L / A [ohm]" \
	"$(awk -F, 'NR == 2 { print $2 }' "$WORK/ofe_transient.csv")" 1.2e-2 1e-8
compare "ofe_transient.csv : last R == the log [ohm]" \
	"$(awk -F, 'END { print $2 }' "$WORK/ofe_transient.csv")" "$(th_r)" 1e-8
#      (b) 対流で冷やす一様な棒 (thermal_conv.ofe)。θ = T - T_a は
#          dθ/dt = a - b θ、a = J² ρ_a / ρc、b = (hP/A - J² ρ0 α) / ρc、θ(0) = 20
set -- $(awk 'BEGIN { J = 6 / 5e-7; ra = 1e-6 * 1.12; a = J * J * ra / 3.4e6
	b = (200 * 3e-3 / 5e-7 - J * J * 1e-6 * 4e-3) / 3.4e6
	printf "%.15e %.15e %.15e", -b, 20 - a / b, 30 + a / b }')
set -- $(th_ode "$1" "$2" "$3" 5 50 1 2e4 4e-3 0)
awk '/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 5 50"} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thconv_tr.ofe"
th_run thconv_tr.ofe
compare "convection-cooled bar, transient : R == the discrete recurrence [ohm]" "$(th_r)" "$1" 1e-9
compare "convection-cooled bar, transient : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
th_balance "convection-cooled bar, transient"
#      (c) ヒートシンクの棒 (data/sample/thermal_transient.ofe) : Fourier 級数解。
#          比べるのは中央の温度上昇 θ_max = T_max - 20 (T_max そのものだと 20 の分だけ
#          誤差が薄まる)
th_four() {	# th_four -> t = 2 s の中央の温度 [degC]
	awk 'BEGIN { PI = atan2(0, -1); k = 10; c = 3.4e6; L = 0.01
		s = (7.5e6 * 7.5e6 * 1e-6) / c; kap = k / c; th = (s / (2 * kap)) * (L / 2) * (L / 2)
		for (nn = 1; nn < 4000; nn += 2) th += -(s / kap) * 4 * L * L / (nn * PI) ^ 3 * sin(nn * PI / 2) * exp(-kap * (nn * PI / L) ^ 2 * 2)
		printf "%.12e", 20 + th }'
}
cp "$SRC/thermal_transient.ofe" "$WORK/"
th_run thermal_transient.ofe
tf=$(th_four)
t40=$(th_t)
#   許容 0.1% : 空間の h² (40 分割で +0.040%、実測)。時間の誤差は 1e-5 以下
compare "heat-sink bar, transient : T_max - 20 vs Fourier series (t = 2 s) [K]" \
	"$(awk -v t="$t40" 'BEGIN { printf "%.12e", t - 20 }')" "$(awk -v t="$tf" 'BEGIN { printf "%.12e", t - 20 }')" 0.001
th_balance "thermal_transient"
awk '/^xmesh/{print "xmesh = 0 20 0.01"; next} {print}' "$SRC/thermal_transient.ofe" > "$WORK/thtr20.ofe"
th_run thtr20.ofe
t20=$(th_t)
res=$(awk -v a="$t20" -v b="$t40" -v e="$tf" 'BEGIN {
	e1 = (a - e) / (e - 20); if (e1 < 0) e1 = -e1
	e2 = (b - e) / (e - 20); if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.2)) ? "OK" : "NG", e1, e2, q }')
echo "  second-order convergence in space, transient (20 -> 40 cells) : $res"
case "$res" in NG*) status=1 ;; esac
#      (d) 長い時間で定常解に行き着く (thermal_bar.ofe を 400 s = 時定数の約 90 倍)。
#          後退 Euler の刻みの式は ΔT -> 0 で定常の式そのものになるので、定常の
#          ソルバーの R と内側の反復の打ち切りの程度で一致しなければならない
th_run thermal_bar.ofe
rss=$(th_r)
awk '/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 400 200"} {print}' \
	"$SRC/thermal_bar.ofe" > "$WORK/thlong.ofe"
th_run thlong.ofe
compare "transient at t = 400 s == the steady solver [ohm]" "$(th_r)" "$rss" 1e-8
#          初期温度 (50) がヒートシンク (20) と違う唯一のケースなので、蓄熱を
#          「初期の温度場」でなく「temperature の値」から測る誤りはここでだけ落ちる
th_balance "thermal_bar, transient to steady (initial 50 degC, sinks 20 degC)"
#      (e) 首のある導体を ヒートシンク + 対流 + 輻射で、CN (θ = 0.5) : 4 経路の収支。
#          (1-θ) 側の項 R(T^n) を取り違えると後退 Euler では見えないので CN で見る
awk '/^thermal/{print "convection = 100 30"; print "radiation = 0.9 10"; print "heatcap = 2 3.4e6"; print "transient = 1 50 0.5"} {print}' \
	"$SRC/thermal_neck.ofe" > "$WORK/thneck_tr.ofe"
th_run thneck_tr.ofe
th_balance "thermal_neck, transient (CN, sinks + convection + radiation)"

# (13) 電流の波形 (waveform = <t0> <f0> <t1> <f1> ...、I(t) = I x f(t) の折れ線)。
#      刻み [t_n, t_{n+1}] の中では f(t_n^+) と f(t_{n+1}^-) を使う
echo "[thermal] transient current waveform : a pulse (exact recurrence, CN order), a ramp"
#      (a) 一定の波形 (waveform = 0 1) は波形なしとバイト一致 (恒等式)
mkdir -p "$WORK/thw1" "$WORK/thw0"
cp "$SRC/thermal_transient.ofe" "$WORK/thw0/"
awk '/^transient/{print; print "waveform = 0 1"; next} {print}' "$SRC/thermal_transient.ofe" \
	> "$WORK/thw1/thermal_transient.ofe"
(cd "$WORK/thw0" && "$OFE" -n 2 thermal_transient.ofe > /dev/null 2>&1)
(cd "$WORK/thw1" && "$OFE" -n 2 thermal_transient.ofe > /dev/null 2>&1)
if [ -s "$WORK/thw0/ofe_transient.csv" ] && cmp -s "$WORK/thw0/ofe.out" "$WORK/thw1/ofe.out" \
	&& cmp -s "$WORK/thw0/ofe_transient.csv" "$WORK/thw1/ofe_transient.csv"; then
	echo "  waveform = 0 1 == no waveform (ofe.out, ofe_transient.csv byte-identical) -> OK"
else
	echo "  waveform = 0 1 differs from no waveform -> NG" >&2; status=1
fi
#      (b) 矩形パルス (data/sample/thermal_pulse.ofe)。一様な棒なので θ 法の離散解は
#          1 次の漸化式で厳密 (段差の両側の電流を刻みの両端で取り分ける)
th_pulse() {	# th_pulse <N> <θ> -> 離散解の R T、連続解の R T
	awk -v ns="$1" -v th="$2" 'function curr(t) { return (t < 2) ? 6 : 0 }
	function curl(t) { return (t <= 2) ? 6 : 0 }
	BEGIN {
		r0 = 1e-6; al = 4e-3; T0 = 0; L = 0.01; A = 5e-7; P = 3e-3; h = 200; Ta = 30; c = 3.4e6; te = 4
		ra = r0 * (1 + al * (Ta - T0)); aa = al / (1 + al * (Ta - T0)); hp = h * P / A
		dt = te / ns; x = 50 - Ta
		for (k = 0; k < ns; k++) {
			J0 = curr(k * dt) / A; J1 = curl((k + 1) * dt) / A
			G0 = J0 * J0 * ra * (1 + aa * x) - hp * x
			x = (x + (dt / c) * ((1 - th) * G0 + th * J1 * J1 * ra)) / (1 - (dt / c) * th * (J1 * J1 * ra * aa - hp))
		}
		J = 6 / A; a1 = J * J * ra / c; b1 = (hp - J * J * r0 * al) / c; b2 = hp / c
		thi = a1 / b1; x2 = thi + (20 - thi) * exp(-b1 * 2); xc = x2 * exp(-b2 * (te - 2))
		printf "%.12e %.12e %.12e %.12e", ra * (1 + aa * x) * L / A, Ta + x, ra * (1 + aa * xc) * L / A, Ta + xc }'
}
for th in 1 0.5; do
	for ns in 64 128; do
		awk -v th="$th" -v ns="$ns" '/^transient/{print "transient = 4 " ns " " th; next} {print}' \
			"$SRC/thermal_pulse.ofe" > "$WORK/thpu$ns.ofe"
		th_run "thpu$ns.ofe"
		eval "r$ns=\$(th_r)"
	done
	set -- $(th_pulse 128 "$th")
	compare "pulse, theta = $th : R == the discrete recurrence [ohm]" "$r128" "$1" 1e-9
	compare "pulse, theta = $th : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
	th_balance "pulse, theta = $th"
	#   段差が刻みの境目にあっても次数が落ちない : 後退 Euler 2、CN 4。
	#   **段差の両側を取り分けずに f(t_{n+1}) をそのまま使うと CN は 1 次に落ちる**
	res=$(awk -v a="$r64" -v b="$r128" -v e="$3" -v th="$th" 'BEGIN {
		e1 = (a - e) / e; if (e1 < 0) e1 = -e1
		e2 = (b - e) / e; if (e2 < 0) e2 = -e2
		q = (e2 > 0) ? e1 / e2 : 0; w = ((th == 1) ? 2 : 4)
		printf "%s (err %.2e -> %.2e, ratio %.3f, expected %d)", ((q > 0.95 * w) && (q < 1.05 * w)) ? "OK" : "NG", e1, e2, q, w }')
	echo "  time order across the jump, theta = $th (64 -> 128 steps) : $res"
	case "$res" in NG*) status=1 ;; esac
done
#          csv の電流の列 : 段差の時刻 t = 2 の行は右極限 (0)、その 1 刻み前は 6
res=$(awk -F, '($1 + 0 == 2) { a = $6 } (($1 + 0) > 1.96) && (($1 + 0) < 1.97) { b = $6 }
	END { printf "%s (I(1.96875) = %s, I(2) = %s)", ((b + 0 == 6) && (a != "") && (a + 0 == 0)) ? "OK" : "NG", b, a }' "$WORK/ofe_transient.csv")
echo "  ofe_transient.csv current column at the jump : $res"
case "$res" in NG*) status=1 ;; esac
#      (c) 直線のランプ (waveform = 0 0 10 1) を断熱の棒 (α = 0) に。
#          ρc dT/dt = J0² ρ0 (t/10)² なので T = T_0 + c0 t³ / 300 (c0 = J0² ρ0 / ρc)。
#          離散解は和で厳密に書ける (後退 Euler : Δt Σ f(t_n)、CN : 台形則)。
#          電流を刻みの始めの時刻で取り違えると後退 Euler の和が 1 項ずれて落ちる
th_ramp() {	# th_ramp <N> <θ> -> 離散解の T、連続解の T
	awk -v ns="$1" -v th="$2" 'BEGIN {
		c0 = 7.5e6 * 7.5e6 * 1e-6 / 3.4e6; dt = 10 / ns; T = 50
		for (k = 1; k <= ns; k++) T += dt * c0 * (th * (k * dt / 10) ^ 2 + (1 - th) * ((k - 1) * dt / 10) ^ 2)
		printf "%.12e %.12e", T, 50 + c0 * 1000 / 300 }'
}
for th in 1 0.5; do
	for ns in 50 100; do
		awk -v th="$th" -v ns="$ns" '/^heatsink/{next} /^tempco/{print "tempco = 2 0 0"; next}
			/^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 10 " ns " " th; print "waveform = 0 0 10 1"; next} {print}' \
			"$SRC/thermal_bar.ofe" > "$WORK/thramp$ns.ofe"
		th_run "thramp$ns.ofe"
		eval "t$ns=\$(th_t)"
	done
	set -- $(th_ramp 100 "$th")
	compare "ramp, theta = $th : T == the discrete sum [degC]" "$t100" "$1" 1e-9
	th_balance "ramp, theta = $th"
	res=$(awk -v a="$t50" -v b="$t100" -v e="$2" -v th="$th" 'BEGIN {
		e1 = (a - e) / (e - 50); if (e1 < 0) e1 = -e1
		e2 = (b - e) / (e - 50); if (e2 < 0) e2 = -e2
		q = (e2 > 0) ? e1 / e2 : 0; w = ((th == 1) ? 2 : 4)
		printf "%s (err %.2e -> %.2e, ratio %.3f, expected %d)", ((q > 0.95 * w) && (q < 1.05 * w)) ? "OK" : "NG", e1, e2, q, w }')
	echo "  time order of the ramp, theta = $th (50 -> 100 steps) : $res"
	case "$res" in NG*) status=1 ;; esac
done
#      (d) 過渡解析の場の出力。電気は単位電流で解いて最後に終了時刻の電流を掛けるので、
#          φ の最大 (ポートの電位) = ログの V = I R、T の最大 = ログの T_max を見る。
#          終了時刻の電流が 1 でない波形 (0.5 倍) にして、掛け忘れを落とす
awk '/^waveform/{print "waveform = 0 1 3 0.5"; print "fieldout = 1"; next} {print}' \
	"$SRC/thermal_pulse.ofe" > "$WORK/thpu_vtk.ofe"
th_run thpu_vtk.ofe
vmax=$(awk '/^SCALARS phi_R_port1 / { f = 1; next } f && /^LOOKUP/ { next }
	f && (NF == 1) { if (($1 + 0) > m) m = $1 + 0; next } f { f = 0 }
	END { printf "%.9e", m }' "$WORK/ofe_field.vtk")
vlog=$(awk '/^  R = / { gsub(/\(V = /, ""); for (i = 1; i <= NF; i++) if ($i == "[V]),") print $(i-1); exit }' "$WORK/ofe.log")
compare "transient field : max phi == V = I R at t_end [V]" "$vmax" "$vlog" 1e-8
tvtk=$(awk '/^SCALARS T / { f = 1; next } f && /^LOOKUP/ { next }
	f && (NF == 1) { if (($1 + 0) > m) m = $1 + 0; next } f { f = 0 }
	END { printf "%.9e", m }' "$WORK/ofe_field.vtk")
compare "transient field : max T == T_max in the log [degC]" "$tvtk" "$(th_t)" 1e-8

# (14) 領域ごとの対流・輻射 (convection / radiation に直方体を付ける。面は
#      あとに書いた行ほど優先)
echo "[thermal] per-face convection / radiation : closed form, ordering, mirror symmetry"
#      (a) 下面だけ冷やす板 (data/sample/thermal_plate.ofe)。α = 0 なので温度は z だけの
#          2 次式で、1 次要素は節点で厳密 : T_max = T_a + qd/h + qd²/(2k) = 145.3125
cp "$SRC/thermal_plate.ofe" "$WORK/"
th_run thermal_plate.ofe
compare "plate cooled from the bottom only : T_max == the closed form (nodally exact) [degC]" \
	"$(th_t)" 145.3125 1e-10
th_energy "plate cooled from the bottom only"
th_nface "plate cooled from the bottom only (10 bottom faces)" 10
#      (b) 同じ板を逆の書き方で : 既定 h = 500、側面 2 枚と上面を h = 0 で上書き。
#          面ごとの係数が同じなら同じ演算になるので ofe.out はバイト一致する。
#          **「最初に当たった行」を取る誤りは、どちらか片方の書き方で必ず落ちる**
mkdir -p "$WORK/thpl_a" "$WORK/thpl_b"
cp "$SRC/thermal_plate.ofe" "$WORK/thpl_a/"
awk '/^convection = 0 30$/{print "convection = 500 30"
		print "convection = 0 30  0 0.01  0 0  0 1e-3"
		print "convection = 0 30  0 0.01  1e-3 1e-3  0 1e-3"
		print "convection = 0 30  0 0.01  0 1e-3  1e-3 1e-3"; next}
	/^convection = 500/{next} {print}' "$SRC/thermal_plate.ofe" > "$WORK/thpl_b/thermal_plate.ofe"
(cd "$WORK/thpl_a" && "$OFE" -n 2 thermal_plate.ofe > /dev/null 2>&1)
(cd "$WORK/thpl_b" && "$OFE" -n 2 thermal_plate.ofe > /dev/null 2>&1)
if [ -s "$WORK/thpl_a/ofe.out" ] && cmp -s "$WORK/thpl_a/ofe.out" "$WORK/thpl_b/ofe.out"; then
	echo "  the same faces written as 'default 0 + bottom 500' and 'default 500 + three zeros' : byte-identical -> OK"
else
	echo "  the two ways of writing the same faces differ -> NG" >&2; status=1
fi
#      (b2) 直方体の許容幅。zmesh = 0 1 1e-4 5 8e-4 の最上面は丸めで
#          0.00079999999999999993 になり、利用者の書く 8e-4 (0.00080000000000000004) と
#          1 ulp ずれる。上面を z = 8e-4 の直方体で断熱にしても、許容幅が無いと
#          上面が冷えたまま残る。厚さ d = 0.8mm の閉形式 T_max = 30 + 140.625 + 2.8125
awk '/^zmesh/{print "zmesh = 0 1 1e-4 5 8e-4"; next}
	/^geometry/{print "geometry = 2 1 0 0.01 0 1e-3 0 8e-4"; next}
	/^conductor = 0/{print "conductor = 0 1 0 0 0 1e-3 0 8e-4"; next}
	/^conductor = 1/{print "conductor = 1 1 0.01 0.01 0 1e-3 0 8e-4"; next}
	/^convection = 0 30$/{print "convection = 500 30"
		print "convection = 0 30  0 0.01  0 0  0 8e-4"
		print "convection = 0 30  0 0.01  1e-3 1e-3  0 8e-4"
		print "convection = 0 30  0 0.01  0 1e-3  8e-4 8e-4"; next}
	/^convection = 500/{next} {print}' "$SRC/thermal_plate.ofe" > "$WORK/thpl_eps.ofe"
th_run thpl_eps.ofe
compare "plate on a grid whose top plane is 1 ulp off the box : T_max == the closed form [degC]" \
	"$(th_t)" 173.4375 1e-10
th_nface "plate on a grid whose top plane is 1 ulp off the box (10 bottom faces)" 10
#      (b3) 直方体が棒の途中で切れる場合。下面の左半分 (x <= 5mm) だけ冷やすものと
#          右半分 (x >= 5mm) だけ冷やすものは x = L/2 の鏡映で合同 (電極は両端、
#          断熱) なので R・T_max が一致し、どちらも 5 面 (α = 4e-3 にして R も温度に
#          依らせる。半分しか冷やさないので 7.5 A では熱暴走する。3 A にする)。**面の係数を中心でなく
#          角の座標で引く誤り**は、境目 x = 5mm の面を片側にだけ入れて落ちる
#          (ここまでの直方体は棒の全長を覆うので、x の位置を見ていなかった)
for side in left right; do
	if [ "$side" = left ]; then box="0 5e-3"; else box="5e-3 0.01"; fi
	awk -v box="$box" '/^convection = 500/{print "convection = 500 30  " box "  0 1e-3  0 0"; next}
		/^tempco/{print "tempco = 2 4e-3 0"; next} /^thermal/{print "thermal = 3"; next} {print}' \
		"$SRC/thermal_plate.ofe" > "$WORK/thpl_$side.ofe"
	th_run "thpl_$side.ofe"
	eval "rr_$side=\$(th_r); tt_$side=\$(th_t)"
	th_nface "plate, bottom cooled on the $side half only" 5
done
compare "plate, left half vs right half cooled (mirror) : R [ohm]" "$rr_right" "$rr_left" 1e-9
compare "plate, left half vs right half cooled (mirror) : T_max [degC]" "$tt_right" "$tt_left" 1e-9
#      (c) 鏡映の恒等式 : 厚さ 1mm (z 2 セル) の棒を全面から冷やしたものと、その半分
#          (厚さ 0.5mm、1 セル、電流も半分) の上面 (鏡面) だけを直方体で h = 0 / ε = 0
#          にしたもの。離散問題が合同なので T_max が一致し、R は半分の棒が 2 倍になる。
#          対流と輻射の両方で見る (輻射は直方体の係数で放射率を引く別の経路)
for kind in convection radiation; do
	if [ "$kind" = convection ]; then bcl="convection = 200 30"; bcm="convection = 0 30  0 0.01  0 1e-3  0.5e-3 0.5e-3"
	else bcl="radiation = 0.8 10"; bcm="radiation = 0 10  0 0.01  0 1e-3  0.5e-3 0.5e-3"; fi
	awk -v bcl="$bcl" '/^zmesh/{print "zmesh = 0 2 1e-3"; next}
		/^geometry/{print "geometry = 2 1 0 0.01 0 1e-3 0 1e-3"; next}
		/^conductor = 0/{print "conductor = 0 1 0 0 0 1e-3 0 1e-3"; next}
		/^conductor = 1/{print "conductor = 1 1 0.01 0.01 0 1e-3 0 1e-3"; next}
		/^convection/{print bcl; next} /^thermal/{print "thermal = 12"; next} {print}' \
		"$SRC/thermal_conv.ofe" > "$WORK/thmir_full.ofe"
	awk -v bcl="$bcl" -v bcm="$bcm" '/^convection/{print bcl; print bcm; next} {print}' \
		"$SRC/thermal_conv.ofe" > "$WORK/thmir_half.ofe"
	th_run thmir_full.ofe
	rf=$(th_r); tf=$(th_t)
	th_nface "$kind, full bar (40 x 6 faces)" 240
	th_run thmir_half.ofe
	compare "$kind, mirror : R(half bar, half current) == 2 R(full bar) [ohm]" "$(th_r)" \
		"$(awk -v r="$rf" 'BEGIN { printf "%.12e", 2 * r }')" 1e-9
	compare "$kind, mirror : T_max(half bar) == T_max(full bar) [degC]" "$(th_t)" "$tf" 1e-9
	th_nface "$kind, half bar with the mirror face switched off (40 x 3 faces)" 120
done

echo "[thermal] tempcotable : a piecewise-linear resistivity in the coupling"
# (15) 抵抗率の温度係数の表 (tempcotable) は熱連成の σ(T) も通る。対流で冷やす一様な棒 (thermal_conv.ofe) は
#       温度が一様なので、釣り合い J² ρ0 f(T) A = hP (T - T_a) を二分法で解いた値と
#       機械精度で一致する。表は 4 点で、解 (約 197 degC) は 2 区間目の途中
awk '/^tempco/{print "tempcotable = 2 0 1 100 1.2 300 1.6 600 3"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/tts_conv.ofe"
th_run tts_conv.ofe
set -- $(awk 'function f(T) { if (T <= 0) return 1; if (T <= 100) return 1 + 0.2 * T / 100
		if (T <= 300) return 1.2 + 0.4 * (T - 100) / 200; if (T <= 600) return 1.6 + 1.4 * (T - 300) / 300; return 3 }
	BEGIN { lo = 30; hi = 2000
		for (it = 0; it < 300; it++) { T = (lo + hi) / 2; if (72 * f(T) - 0.6 * (T - 30) > 0) lo = T; else hi = T }
		printf "%.12e %.12e", 0.02 * f(T), T }')
compare "tempcotable in thermal (convection-cooled bar) : R == the equilibrium [ohm]" "$(th_r)" "$1" 1e-9
compare "tempcotable in thermal (convection-cooled bar) : T == the equilibrium [degC]" "$(th_t)" "$2" 1e-9
#       表の外だけを通る棒 (どの温度でも f = 2) は α = 0 で抵抗率 2ρ0 の棒と同じで、
#       放物線が節点で厳密 : R = 2e-2、T_max = 20 + J² (2ρ0) L² / (8k) = 160.625
awk '/^tempco/{print "tempcotable = 2 -100 1.5 -50 2"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/tts_flat.ofe"
th_run tts_flat.ofe
compare "tempcotable beyond its last point in thermal : R = 2 rho0 L / A [ohm]" "$(th_r)" 2e-2 1e-10
compare "tempcotable beyond its last point in thermal : T_max (parabola) [degC]" "$(th_t)" 160.625 1e-10


# (7) 弾くべき入力 (理由まで grep する)
th_reject() {	# th_reject <ラベル> <ofe> <理由の文字列>
	mesh_reject "$1" "$2"
	if ! (cd "$WORK" && "$OFE" -n 2 "$2" 2>&1 | grep -q "$3"); then
		echo "  *** it was rejected for a different reason than: $3" >&2
		status=1
	fi
}
awk '/^analysis/{print "analysis = C"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thx_c.ofe"
th_reject "thermal with analysis C" thx_c.ofe "analysis = R alone"
awk '/^conductor = 1/{print; print "conductor = 2 1 0.005 0.005 0 1e-3 0 1e-3"; next} {print}' \
    "$SRC/thermal_bar.ofe" > "$WORK/thx_2p.ofe"
th_reject "thermal with two ports" thx_2p.ofe "exactly one port"
awk '!/^heatsink/' "$SRC/thermal_bar.ofe" > "$WORK/thx_nosink.ofe"
th_reject "thermal without a heat sink" thx_nosink.ofe "at least one heatsink"
awk '!/^kthermal/' "$SRC/thermal_bar.ofe" > "$WORK/thx_nok.ofe"
th_reject "a conductor without kthermal" thx_nok.ofe "no thermal conductivity"
awk '/^thermal/{print "frequency = 1e3"} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thx_freq.ofe"
th_reject "thermal with frequency" thx_freq.ofe "DC coupling"
awk '/^heatsink = 1/{print "heatsink = 3 20"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/thx_badsink.ofe"
th_reject "heatsink on a conductor that does not exist" thx_badsink.ofe "does not exist"
#     高次の四面体でも解ける (閉形式は [thermal higher-order tet])。ここでは 2 次の箱の恒等式だけ
cp "$SRC/box_p2.msh" "$WORK/"
awk '/^material/{print "material = 1.0 1e6"; print "kthermal = 2 10"; next}
	/^analysis/{print "heatsink = 0 20"; print "thermal = 1"; print "analysis = R"; next}
	{print}' "$SRC/box_p2.ofe" > "$WORK/thx_p2.ofe"
th_run thx_p2.ofe
th_energy "thermal on a second-order tetrahedral mesh"
#     h = 0 の対流は熱の逃げ道にならない (ヒートシンクも無ければ定常解が無い)
awk '/^convection/{print "convection = 0 30"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thx_h0.ofe"
th_reject "convection with h = 0 and no heat sink" thx_h0.ofe "at least one heatsink"
awk '/^convection/{print "convection = -200 30"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thx_hneg.ofe"
th_reject "convection with h < 0" thx_hneg.ofe "invalid convection data"
awk '/^convection/{print "convection = 200"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thx_h1.ofe"
th_reject "convection without the ambient temperature" thx_h1.ofe "invalid convection data"
awk '/^convection/{print "radiation = 0 10"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thx_e0.ofe"
th_reject "radiation with emissivity 0 and no heat sink" thx_e0.ofe "at least one heatsink"
awk '/^radiation/{print "radiation = 1.2 10"; next} {print}' "$SRC/thermal_rad.ofe" > "$WORK/thx_e12.ofe"
th_reject "radiation with emissivity > 1" thx_e12.ofe "invalid radiation data"
awk '/^radiation/{print "radiation = 0.5 -300"; next} {print}' "$SRC/thermal_rad.ofe" > "$WORK/thx_tk.ofe"
th_reject "radiation to a wall below absolute zero" thx_tk.ofe "invalid radiation data"
awk '!/^thermal/' "$SRC/thermal_transient.ofe" > "$WORK/thx_trnoth.ofe"
th_reject "transient without thermal" thx_trnoth.ofe "transient needs thermal"
awk '/^transient/{print "transient = 2 200 0.4"; next} {print}' "$SRC/thermal_transient.ofe" > "$WORK/thx_th04.ofe"
th_reject "transient with theta < 0.5 (conditionally stable)" thx_th04.ofe "invalid transient data"
awk '/^transient/{print "transient = 2 0"; next} {print}' "$SRC/thermal_transient.ofe" > "$WORK/thx_ns0.ofe"
th_reject "transient with no time step" thx_ns0.ofe "invalid transient data"
awk '!/^heatcap/' "$SRC/thermal_transient.ofe" > "$WORK/thx_noc.ofe"
th_reject "transient with a material without heatcap" thx_noc.ofe "no heat capacity"
awk '!/^transient/' "$SRC/thermal_pulse.ofe" > "$WORK/thx_wnotr.ofe"
th_reject "waveform without transient" thx_wnotr.ofe "waveform needs transient"
awk '/^waveform/{print "waveform = 0 1 2 1 1 0"; next} {print}' "$SRC/thermal_pulse.ofe" > "$WORK/thx_wdec.ofe"
th_reject "waveform with decreasing times" thx_wdec.ofe "invalid waveform data"
awk '/^waveform/{print "waveform = 0 1 2"; next} {print}' "$SRC/thermal_pulse.ofe" > "$WORK/thx_wodd.ofe"
th_reject "waveform with an odd number of values" thx_wodd.ofe "invalid waveform data"
awk '/^convection = 500/{print "convection = 500 30  0 0.01  0 1e-3  -1 -1"; next} {print}' \
	"$SRC/thermal_plate.ofe" > "$WORK/thx_box0.ofe"
#     (解く段で見つかるので理由はログに出る)
mesh_reject "a convection box that contains no exposed face" thx_box0.ofe
if ! grep -q "contains no exposed face" "$WORK/ofe.log"; then
	echo "  *** it was rejected for a different reason than the empty box" >&2
	status=1
fi
awk '/^convection = 500/{print "convection = 500 30  0.01 0  0 1e-3  0 0"; next} {print}' \
	"$SRC/thermal_plate.ofe" > "$WORK/thx_boxinv.ofe"
th_reject "a convection box with x1 > x2" thx_boxinv.ofe "invalid convection data"
awk '/^convection = 500/{print "convection = 500 30  0 0.01  0 1e-3"; next} {print}' \
	"$SRC/thermal_plate.ofe" > "$WORK/thx_box4.ofe"
th_reject "a convection box with four coordinates" thx_box4.ofe "invalid convection data"
awk '/^radiation/{print "radiation = 1.5 10  0 0.01  0 1e-3  0 0"; next} {print}' \
	"$SRC/thermal_rad.ofe" > "$WORK/thx_radbox.ofe"
th_reject "a radiation box with emissivity > 1" thx_radbox.ofe "invalid radiation data"

# ---- 電気-熱連成の非構造格子 (1 次四面体) ----
#
# data/sample/thermal_tet.ofe : thermal_bar.ofe と同じ棒 (10mm x 1mm x 1mm) を四面体で
# 切ったもの (thermal_tet20/40.msh、mkmesh.py bar)。構造格子の閉形式 (th_exact /
# th_fin / th_radfin / th_adia / th_four) がそのまま使える。熱の反復は格子に依らず、
# 四面体で分かれるのは剛性・発熱 (σ|∇φ|² V/4)・質量 (ρc V (1+δ)/20)・要素の温度
# (4 節点の平均)・露出面 (三角形) だけ。
#
# **断面の温度の一様性に注意。** 構造格子の断面 1 x 1 セルは節点が 4 隅だけなので
# 対称性から断面で温度が一様になり、1 次元の閉形式 (フィン近似) が離散解の極限に
# なる。四面体の格子 (y 1 x z 2 分割) は側面の中央に節点があるので、側面から
# 冷やすと断面の温度勾配 (Bi = hw/k) を本当に解く。h = 200 (Bi = 0.02) では
# T_max が 1 次元の式から +0.18% ずれて細分しても縮まない (実測) — 誤差ではなく
# 3 次元の物理なので、側面から冷やす閉形式は Bi の小さい h = 20 で比べる
echo "[thermal_tet] electro-thermal coupling on linear tetrahedra"
cp "$SRC/thermal_tet.ofe" "$SRC/thermal_tet20.msh" "$SRC/thermal_tet40.msh" \
	"$SRC/bar_tet.msh" "$SRC/bar_air.msh" "$WORK/"
tt_mesh() {	# tt_mesh <nx> <入力> <出力> : 格子を thermal_tet<nx>.msh に差し替える
	awk -v nx="$1" '/^mesh/{print "mesh = thermal_tet" nx ".msh"; next} {print}' "$WORK/$2" > "$WORK/$3"
}
tt_order() {	# tt_order <ラベル> <20 分割の値> <40 分割の値> <閉形式> [<基準を引く量>]
	res=$(awk -v a="$2" -v b="$3" -v e="$4" -v z="${5:-0}" 'BEGIN {
		e1 = (a - e) / (e - z); if (e1 < 0) e1 = -e1
		e2 = (b - e) / (e - z); if (e2 < 0) e2 = -e2
		q = (e2 > 0) ? e1 / e2 : 0
		printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.2)) ? "OK" : "NG", e1, e2, q }')
	echo "  $1 (20 -> 40 cells) : $res"
	case "$res" in NG*) status=1 ;; esac
}
# (1) 両端ヒートシンク (I = 7.5 A, mL = 1.5) : 閉形式・恒等式・誤差の次数。
#     許容 0.1% : h² の離散化誤差 (40 分割で R -0.020% / T_max -0.032%、実測)
th_run thermal_tet.ofe
set -- $(th_exact 7.5 4e-3)
compare "tet : R(I = 7.5 A, mL = 1.5) [ohm]" "$(th_r)" "$1" 0.001
compare "tet : T_max(I = 7.5 A) [degC]" "$(th_t)" "$2" 0.001
th_energy "tet : thermal_tet"
r40=$(th_r); rex=$1
tt_mesh 20 thermal_tet.ofe tt20.ofe
th_run tt20.ofe
tt_order "tet : second-order convergence of R" "$(th_r)" "$r40" "$rex"
#     場の出力の温度の最大 = ログの T_max
awk '/^thermal/{print; print "fieldout = 1"; next} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tt_fo.ofe"
th_run tt_fo.ofe
tvtk=$(awk '/^SCALARS T / { f = 1; next } f && /^LOOKUP/ { next }
	f && (NF == 1) { if (($1 + 0) > m) m = $1 + 0; next } f { f = 0 }
	END { printf "%.9e", m }' "$WORK/ofe_field.vtk")
compare "tet : max of the T field == T_max in the log [degC]" "$tvtk" "$(th_t)" 1e-8
# (2) α = 0 : R は熱連成なしの値に厳密に一致し、この格子では T_max も放物線の頂点に
#     全桁一致する (x の各断面に節点が並ぶので、一様な発熱の節点荷重が断面ごとに揃う)
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tt_a0.ofe"
th_run tt_a0.ofe
set -- $(th_exact 7.5 0)
compare "tet : R(alpha = 0) == the uncoupled value [ohm]" "$(th_r)" "$1" 1e-10
compare "tet : T_max(alpha = 0) == parabola, nodally exact [degC]" "$(th_t)" "$2" 1e-10
# (3) ヒートシンク + 側面の対流 (フィン、h = 20 で Bi = 0.002)。面の数は
#     y 側面 2 x 40 x 2 + z 側面 2 x 40 x 1 の四角形 = 240 枚 x 三角形 2 = 480
#     (電極の面は除く)。ヒートシンクと対流の両方があるので、対流の右辺を
#     熱の収支に入れ忘れるとエネルギーの恒等式が落ちる
awk '/^thermal/{print "convection = 20 30"} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tt_fin.ofe"
th_run tt_fin.ofe
set -- $(th_fin 7.5 20)
compare "tet : R(fin, h = 20) [ohm]" "$(th_r)" "$1" 0.001
compare "tet : T_max(fin, h = 20) [degC]" "$(th_t)" "$2" 0.001
th_energy "tet : fin, h = 20"
th_nface "tet : fin (sides only, not the electrode faces)" 480
r40=$(th_r); rex=$1
tt_mesh 20 tt_fin.ofe tt_fin20.ofe
th_run tt_fin20.ofe
tt_order "tet : second-order convergence, fin h = 20" "$(th_r)" "$r40" "$rex"
#     h = 200 は 3 次元の効果で 1 次元の式と比べられない (上の注意) ので恒等式だけ
awk '/^thermal/{print "convection = 200 30"} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tt_fin200.ofe"
th_run tt_fin200.ofe
th_energy "tet : fin, h = 200"
#     直方体で y = 0 の側面だけに対流 : 40 x 2 x 2 = 160 枚。面の中心でなく頂点で
#     判定すると z 側面の三角形が混ざる
awk '/^thermal/{print "convection = 200 30  0 0.01  0 0  0 1e-3"} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/tt_box.ofe"
th_run tt_box.ofe
th_energy "tet : convection on the y = 0 face only"
th_nface "tet : convection box on y = 0" 160
# (4) ヒートシンク + 輻射 (thermal_rad.ofe と同じ k = 1, ε = 1, I = 2.5 A)。
#     k = 1 なので断面の勾配が T_max に +0.09% 出る (R には -0.008%、実測)。
#     α = 0 では Newton の 2 次収束で 4 回 (実測)。線形化の誤りは同じ不動点へ
#     遅く収束するだけなので回数を assert する
awk '/^kthermal/{print "kthermal = 2 1"; next} /^thermal/{print "radiation = 1 10"; print "thermal = 2.5"; next} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/tt_rad.ofe"
th_run tt_rad.ofe
set -- $(th_radfin)
compare "tet : R(heat sinks + radiation) [ohm]" "$(th_r)" "$1" 0.001
compare "tet : T_max(heat sinks + radiation) [degC]" "$(th_t)" "$2" 0.002
th_energy "tet : heat sinks + radiation"
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$WORK/tt_rad.ofe" > "$WORK/tt_rad0.ofe"
th_run tt_rad0.ofe
nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
if [ -n "$nit" ] && [ "$nit" -le 5 ]; then
	echo "  tet : Newton on radiation (alpha = 0) : converged in $nit iterations (<= 5) -> OK"
else
	echo "  tet : Newton on radiation (alpha = 0) : ${nit:-?} iterations (> 5) -> NG" >&2
	status=1
fi
# (5) 過渡解析。断熱の棒は発熱 σ|∇φ|² V/4 と質量の行和 ρc V/4 が同じ節点の重みなので
#     温度が一様のまま進み、θ 法の漸化式に機械精度で一致する (四面体でも厳密)。
#     質量の (1 + δ) の重みを取り違えると行和が崩れてここで落ちる
for th in 1 0.5; do
	awk -v th="$th" '/^heatsink/{next}
		/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 10 100 " th} {print}' \
		"$SRC/thermal_tet.ofe" > "$WORK/tt_ad.ofe"
	th_run tt_ad.ofe
	set -- $(th_adia 100 "$th")
	compare "tet : adiabatic bar, theta = $th : R == the discrete recurrence [ohm]" "$(th_r)" "$1" 1e-9
	compare "tet : adiabatic bar, theta = $th : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
	th_balance "tet : adiabatic bar, theta = $th"
done
#     ヒートシンクの棒 (thermal_transient.ofe と同じ) : Fourier 級数解と空間の次数。
#     許容 0.1% : 40 分割で +0.040% (実測)
awk '/^tempco/{print "tempco = 2 0 0"; next} /^temperature/{print "temperature = 20"; next}
	/^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 2 200 0.5"; next} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/tt_tr.ofe"
th_run tt_tr.ofe
tf=$(th_four); t40=$(th_t)
compare "tet : heat-sink bar, transient : T_max - 20 vs Fourier series [K]" \
	"$(awk -v t="$t40" 'BEGIN { printf "%.12e", t - 20 }')" "$(awk -v t="$tf" 'BEGIN { printf "%.12e", t - 20 }')" 0.001
th_balance "tet : heat-sink bar, transient"
tt_mesh 20 tt_tr.ofe tt_tr20.ofe
th_run tt_tr20.ofe
tt_order "tet : second-order convergence in space, transient" "$(th_t)" "$t40" "$tf" 20
#     CN でヒートシンク + 対流 + 輻射 : 4 経路の収支 ((1-θ) 側の項は CN でしか見えない)
awk '/^thermal/{print "convection = 100 30"; print "radiation = 0.9 10"; print "heatcap = 2 3.4e6"; print "transient = 1 50 0.5"} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/tt_trall.ofe"
th_run tt_trall.ofe
th_balance "tet : transient (CN, sinks + convection + radiation)"
# (6) 熱を通さない領域 (空気、k = 0) : bar_air.msh は bar_tet.msh と同じ導体の格子に
#     空気層を載せたもの。導体と空気の界面は露出面になるので、空気の無い棒と
#     同じ答え・同じ露出面の数 (y 側面 2 x 24 x 24 + z 側面 2 x 24 x 2 の四角形 x 2
#     = 2496) にならなければならない (節点番号が違うので許容は反復の打ち切りの分)
awk '/^mesh/{print "mesh = bar_tet.msh"; next} /^heatsink = 1/{next}
	/^thermal/{print "convection = 2000 30"; print "thermal = 5"; next} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/tt_bar.ofe"
th_run tt_bar.ofe
rb=$(th_r); tb=$(th_t)
th_nface "tet : bar without air" 2496
awk '/^region/{print; print "region = 2 1"; next} /^mesh/{print "mesh = bar_air.msh"; next} {print}' \
	"$WORK/tt_bar.ofe" > "$WORK/tt_air.ofe"
th_run tt_air.ofe
compare "tet : bar + air (k = 0) == bar alone : R [ohm]" "$(th_r)" "$rb" 1e-9
compare "tet : bar + air (k = 0) == bar alone : T_max [degC]" "$(th_t)" "$tb" 1e-9
th_nface "tet : bar + air, the conductor-air interface is exposed" 2496
th_energy "tet : bar + air"

# ---- 4 次要素 (gmsh の tet35 = 型 30 / tri15 = 型 23) ----
#
# 局所節点の並びは gmsh 4.15.2 の出力を重心座標で分類した実測の表 (unstruct.c の TET4_A /
# TRI4_A)。形状関数は 3 次と同じ多重指数の Lagrange の式、積分は Duffy + 5 点。
#   (a) 箱 (box_p4.msh) の平行平板 : C が閉形式に全桁 (φ が 1 次)。2.2 / 4.1 / バイナリで一致
#   (b) analysis = P : 4 次の φ の恒等式 (g) が機械精度。**4 次特有の基底を検査できるのは
#       これだけ** (1〜3 次の場は 3 次の基底でも厳密に表せる)
#   (c) 曲がった同軸 (coax_p4.msh、h = 0.4mm) : 円の等間隔補間の幾何誤差は 4 次で h⁶
#       (偶数次の超収束。3 次・2 次は h⁴) なので、**同じ h の 3 次より有意に良い**ことを
#       assert してよい (3 次と 2 次の間ではしてはいけない。[order-3 tet] の注意)。
#       実測 C' : 3 次 3.2e-4 -> 4 次 8.5e-6 (38 倍)、体積 : 6.5e-5 -> 9.5e-7 (68 倍)。
#       h = 0.4 -> 0.2mm の比は C' 61、体積 66 (h⁶ = 64。0.2mm の格子は 1.3MB なので置かない)
#   (d) 断面 2 次元 (plate2d_p4.msh、plate2d_p3.geo を -order 4) : 10 MHz の表皮効果で
#       R -1.3e-8 / L +1.3e-9 (同じ頂点の 3 次は +6.7e-8 / +1.8e-7)。L で 10 倍を assert
#       (R は csv の 9 桁の丸めに届くので比を取らない)
echo "[order-4] 35-node tetrahedra / 15-node triangles (gmsh tet35 / tri15)"
for m in box_p4 box_p4_41 box_p4_bin_22 box_p4_bin_41 coax_p4 coax_p3 plate2d_p4 plate2d_p3; do cp "$SRC/$m.msh" "$WORK/"; done
sed 's/^mesh = .*/mesh = box_p4.msh/' "$SRC/box_p3.ofe" > "$WORK/box_p4.ofe"
(cd "$WORK" && "$OFE" -n 2 box_p4.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "C (tet35 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
grep -v '^title' "$WORK/rlc.csv" > "$WORK/p4_22.csv"
for f in box_p4_41 box_p4_bin_22 box_p4_bin_41; do
	sed "s/^mesh = .*/mesh = $f.msh/" "$WORK/box_p4.ofe" > "$WORK/p4f.ofe"
	(cd "$WORK" && "$OFE" -n 2 p4f.ofe > /dev/null && "$OFE_POST" > /dev/null)
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/p4f.csv"
	if cmp -s "$WORK/p4_22.csv" "$WORK/p4f.csv"; then echo "  35-node tetrahedra, $f == ASCII 2.2 -> OK"
	else echo "  35-node tetrahedra, $f differs from ASCII 2.2 -> NG" >&2; status=1; fi
done
if ! head -c 24 "$WORK/box_p4_bin_22.msh" | tr '\n' ' ' | grep -qa 'MeshFormat [0-9.]* 1 8'; then
	echo "  box_p4_bin_22.msh is not a binary Gmsh file -> NG" >&2; status=1
fi
run_case nodal_test_p4
if grep -q "(g) quartic" "$WORK/ofe.log" && grep -q "result : passed" "$WORK/ofe.log"; then
	echo "  P4 self test (quartic field identity, machine precision) : $(awk '/\(g\) quartic/ { print $NF; exit }' "$WORK/ofe.log") -> OK"
else
	echo "  P4 self test : the quartic-field identity was not run or failed -> NG" >&2
	status=1
fi
# (c) 同軸 : 4 次と 3 次 (同じ h = 0.4mm) の C' と体積
p4c() {	# p4c <格子> -> C' と積分した体積
	sed "s/^mesh = .*/mesh = $1.msh/" "$SRC/coax_p3.ofe" > "$WORK/p4c.ofe"
	(cd "$WORK" && "$OFE" -n 2 p4c.ofe > /dev/null && "$OFE_POST" > /dev/null)
	cc=$(value_of C)
	sed -e "s/^mesh = .*/mesh = $1.msh/" -e 's/^analysis = .*/analysis = P/' -e '/^tline/d' \
		"$SRC/coax_p3.ofe" > "$WORK/p4v.ofe"
	(cd "$WORK" && "$OFE" -n 2 p4v.ofe > /dev/null 2>&1)
	echo "$cc $(awk '/integrated volume/ { print $4; exit }' "$WORK/ofe.log")"
}
set -- $(p4c coax_p4); c4=$1; v4=$2
set -- $(p4c coax_p3); c3=$1; v3=$2
compare "C (tet35 coax, h = 0.4 mm) [F/m]" "$c4" 1.063417e-10 0.0001
res=$(awk -v a="$c3" -v b="$c4" -v va="$v3" -v vb="$v4" 'BEGIN { ce = 2 * atan2(0, -1) * 8.8541878128e-12 * 2.1 / log(3)
	ve = atan2(0, -1) * (1.5e-3 ^ 2 - 0.5e-3 ^ 2) * 0.1e-3
	e3 = (a - ce) / ce; if (e3 < 0) e3 = -e3
	e4 = (b - ce) / ce; if (e4 < 0) e4 = -e4
	w3 = (va - ve) / ve; if (w3 < 0) w3 = -w3
	w4 = (vb - ve) / ve; if (w4 < 0) w4 = -w4
	ok = (e3 > 10 * e4) && (w3 > 10 * w4)
	printf "%s (C : P3 %.2e vs P4 %.2e, volume : P3 %.2e vs P4 %.2e)", ok ? "OK" : "NG", e3, e4, w3, w4 }')
echo "  curved coax, order 4 (h^6 geometry) beats order 3 (h^4) at the same h : $res"
case "$res" in NG*) status=1 ;; esac
# (d) 断面 2 次元
set -- $(awk 'BEGIN { pi = atan2(0, -1); mu = 4e-7 * pi; s = 5.8e7; W = 1e-3; t = 0.05e-3; d = 0.2e-3
	om = 2 * pi * 1e7; de = sqrt(2 / (om * mu * s)); a = t / de
	S = (exp(2 * a) - exp(-2 * a)) / 2; C = (exp(2 * a) + exp(-2 * a)) / 2
	k = 1 / (de * s * W * (C - cos(2 * a)))
	printf "%.12e %.12e", 2 * k * (S + sin(2 * a)), (2 * k * (S - sin(2 * a)) + om * mu * d / W) / om }')
rex=$1; lex=$2
for o in 4 3; do
	sed -e "s/plate2d.msh/plate2d_p$o.msh/" -e "s/^frequency = .*/frequency = 1e7/" "$SRC/plate2d_ac.ofe" > "$WORK/p4ac.ofe"
	(cd "$WORK" && "$OFE" -n 2 p4ac.ofe > /dev/null && "$OFE_POST" > /dev/null)
	eval "r$o=\$(value_of Rf); l$o=\$(value_of Lf)"
done
compare "R(f = 10 MHz, 2-D P4) vs coth closed form [ohm/m]" "$r4" "$rex" 1e-7
compare "L(f = 10 MHz, 2-D P4) vs coth closed form [H/m]" "$l4" "$lex" 1e-7
res=$(awk -v a="$l3" -v b="$l4" -v l="$lex" 'BEGIN {
	f3 = (a - l) / l; if (f3 < 0) f3 = -f3
	f4 = (b - l) / l; if (f4 < 0) f4 = -f4
	printf "%s (L : P3 %.2e vs P4 %.2e)", (f3 > 10 * f4) ? "OK" : "NG", f3, f4 }')
echo "  2-D, order 4 beats order 3 on the same vertices (10 MHz) : $res"
case "$res" in NG*) status=1 ;; esac
sed "s/^mesh = .*/mesh = plate2d_p4.msh/" "$SRC/plate2d_dc.ofe" > "$WORK/p4dc.ofe"
(cd "$WORK" && "$OFE" -n 2 p4dc.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "Ldc (2-D P4) [H/m]" "$(value_of Ldc)" 2.932153e-07 1e-5
# (e) 弾くべき入力 : 4 次の四面体の面に 3 次の三角形 (型 23 を型 21 の 10 節点に落とす)
awk 'BEGIN { e = 0 } /^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 } e && ($2 == 23) { s = $1 " 21"; for (i = 3; i <= 3 + $3 + 10; i++) s = s " " $i; print s; next } { print }' "$WORK/box_p4.msh" > "$WORK/box_p4_tri10.msh"
sed 's/^mesh = .*/mesh = box_p4_tri10.msh/' "$WORK/box_p4.ofe" > "$WORK/p4x_tri.ofe"
th_reject "35-node tetrahedra with 10-node electrode triangles" p4x_tri.ofe "triangles are supported only as the faces of 20-node"
#   analysis E (辺要素の自己検証) は [A order 3+] の (h) で通す

# ---- 5 次要素 (gmsh の tet56 = 型 31 / tri21 = 型 25) ----
#
# 並びは gmsh 4.15.2 の出力を重心座標で分類した実測の表 (unstruct.c の TET5_A / TRI5_A。同じ分類で TET4_A /
# TRI4_A が全項一致)。積分は Duffy + 6 点 (216 点)、自己検証の右辺は 7 点。格子は gmsh_geo_order.py。
#   (a) 箱の平行平板 : C が閉形式に全桁。2.2 / 4.1 / バイナリ 4.1 で rlc.csv が一致
#   (b) analysis = P : 5 次の φ の恒等式 (h) と重心の勾配が機械精度 (**5 次特有の基底を見るのはこれだけ**)
#   (c) 曲がった同軸 (h = 0.4mm) : C' 2e-5・体積 5e-6 以内 (実測 −6.7e-6 / +1.9e-6)。奇数次の幾何誤差は 1 つ下の
#       偶数次と同じ h⁶ なので、**4 次に勝つことは assert しない** (実測 4 次 +8.5e-6 / −9.5e-7、偽陽性になる)
#   (d) 断面 2 次元の 10 MHz の表皮効果 : R が閉形式に 2e-9 以内 (実測 −2.9e-10)、同じ頂点の 4 次 (−1.3e-8) より 5 倍以上よい
#   (e) 5 次の四面体の面に 4 次の三角形、E は弾く。面の多重指数を入れ替えると写像が裏返って setup で落ちる (実測)
echo "[order-5] 56-node tetrahedra / 21-node triangles (gmsh tet56 / tri21)"
for m in box_p5 box_p5_41 box_p5_bin_41 coax_p5 plate2d_p5; do cp "$SRC/$m.msh" "$WORK/"; done
sed 's/^mesh = .*/mesh = box_p5.msh/' "$SRC/box_p3.ofe" > "$WORK/box_p5.ofe"
(cd "$WORK" && "$OFE" -n 2 box_p5.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "C (tet56 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
grep -v '^title' "$WORK/rlc.csv" > "$WORK/p5_22.csv"
for f in box_p5_41 box_p5_bin_41; do
	sed "s/^mesh = .*/mesh = $f.msh/" "$WORK/box_p5.ofe" > "$WORK/p5f.ofe"
	(cd "$WORK" && "$OFE" -n 2 p5f.ofe > /dev/null && "$OFE_POST" > /dev/null)
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/p5f.csv"
	if cmp -s "$WORK/p5_22.csv" "$WORK/p5f.csv"; then echo "  56-node tetrahedra, $f == ASCII 2.2 -> OK"
	else echo "  56-node tetrahedra, $f differs from ASCII 2.2 -> NG" >&2; status=1; fi
done
run_case nodal_test_p5
if grep -q "(h) quintic" "$WORK/ofe.log" && grep -q "centroid gradient (quintic" "$WORK/ofe.log" && grep -q "result : passed" "$WORK/ofe.log"; then
	echo "  P5 self test (quintic field identity, machine precision) : $(awk '/\(h\) quintic/ { print $NF; exit }' "$WORK/ofe.log") -> OK"
else
	echo "  P5 self test : the quintic-field identity was not run or failed -> NG" >&2
	status=1
fi
set -- $(p4c coax_p5)
res=$(awk -v c="$1" -v v="$2" 'BEGIN { ce = 2 * atan2(0, -1) * 8.8541878128e-12 * 2.1 / log(3)
	ve = atan2(0, -1) * (1.5e-3 ^ 2 - 0.5e-3 ^ 2) * 0.1e-3
	e = (c - ce) / ce; w = (v - ve) / ve; ae = (e < 0) ? -e : e; aw = (w < 0) ? -w : w
	printf "%s (C %+.2e, volume %+.2e)", ((ae <= 2e-5) && (aw <= 5e-6)) ? "OK" : "NG", e, w }')
echo "  curved coax, order 5 (h = 0.4 mm) : C' within 2e-5 and the isoparametric volume within 5e-6 : $res"
case "$res" in NG*) status=1 ;; esac
sed -e "s/plate2d.msh/plate2d_p5.msh/" -e "s/^frequency = .*/frequency = 1e7/" "$SRC/plate2d_ac.ofe" > "$WORK/p5ac.ofe"
(cd "$WORK" && "$OFE" -n 2 p5ac.ofe > /dev/null && "$OFE_POST" > /dev/null)
r5=$(value_of Rf)
compare "L(f = 10 MHz, 2-D P5) vs coth closed form [H/m]" "$(value_of Lf)" "$lex" 1e-7
res=$(awk -v a="$r4" -v b="$r5" -v r="$rex" 'BEGIN {
	f4 = (a - r) / r; if (f4 < 0) f4 = -f4
	f5 = (b - r) / r; if (f5 < 0) f5 = -f5
	printf "%s (R : P4 %.2e vs P5 %.2e)", ((f5 <= 2e-9) && (f4 > 5 * f5)) ? "OK" : "NG", f4, f5 }')
echo "  2-D, order 5 : R(10 MHz) within 2e-9 of the closed form and 5x better than order 4 : $res"
case "$res" in NG*) status=1 ;; esac
sed "s/^mesh = .*/mesh = plate2d_p5.msh/" "$SRC/plate2d_dc.ofe" > "$WORK/p5dc.ofe"
(cd "$WORK" && "$OFE" -n 2 p5dc.ofe > /dev/null && "$OFE_POST" > /dev/null)
compare "Ldc (2-D P5) [H/m]" "$(value_of Ldc)" 2.932153e-07 1e-5
awk 'BEGIN { e = 0 } /^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 } e && ($2 == 25) { s = $1 " 23"; for (i = 3; i <= 3 + $3 + 15; i++) s = s " " $i; print s; next } { print }' "$WORK/box_p5.msh" > "$WORK/box_p5_tri15.msh"
sed 's/^mesh = .*/mesh = box_p5_tri15.msh/' "$WORK/box_p5.ofe" > "$WORK/p5x_tri.ofe"
mesh_reject "56-node tetrahedra with 15-node electrode triangles" p5x_tri.ofe
#   analysis E (辺要素の自己検証) は [A order 3+] の (h) で通す

# ---- 6 次 / 7 次要素 (gmsh の tet84 = 型 71 / tet120 = 型 72 / tri28 = 型 42 / tri36 = 型 43) ----
#
# 並びは表でなく gmsh の並びの**再帰の規則** (unstruct.c の tet_rule / tri_rule : 頂点、辺、面に次数 p−3 の三角形、
# 内部に次数 p−4 の四面体)。規則は gmsh 4.15.2 の出力 (立方体の 1140 要素 / 540 面) を重心座標で分類したものと
# 2〜7 次の全項で一致し (gmsh_ho_order.py)、起動時に 3〜5 次の実測の表と全項一致することを確かめる。積分は
# Duffy + p+1 点 (8 点以上の Gauss-Legendre は Newton で根を求める。以前は表に無い点数を黙って 5 点に落としていた)。
#   (a) 小さな箱 (box_ho.geo、24 要素) の平行平板 : C が閉形式に全桁、バイナリ 2.2 (6 次) / 4.1 (7 次) で rlc.csv が一致
#   (b) analysis = P : 6 次 (i) / 7 次 (j) の φ の恒等式と重心の勾配 (**その次数特有の基底を見るのはこれだけ**)
#   (c) 曲がった同軸 (6 次、h = 0.8mm) : 幾何誤差が h⁸ なので、**倍の要素の大きさで 4 次 (h = 0.4mm) に勝つ**
#       (実測 C' −2.4e-6 / 体積 +4.7e-7 に対し 4 次 +8.5e-6 / −9.5e-7)
#   (d) 断面 2 次元の 400 MHz の表皮効果 (10 MHz では 5 次が rlc.csv の桁 (3e-10) に届いていて差が見えない) :
#       R の誤差が 4 次 8.0e-4、5 次 5.4e-5、6 次 8.2e-7、7 次 5.0e-8 (実測)。6 次は 5 次の 10 倍、7 次は 6 次の 5 倍よい
#   (e) 6 次の四面体の面に 5 次の三角形、E は弾く
echo "[order-6/7] 84 / 120-node tetrahedra, 28 / 36-node triangles (gmsh tet84 / tet120 / tri28 / tri36)"
for m in box_p6 box_p6_bin_22 box_p7 box_p7_bin_41 coax_p6 plate2d_p6 plate2d_p7; do cp "$SRC/$m.msh" "$WORK/"; done
for f in box_p6:box_p6_bin_22 box_p7:box_p7_bin_41; do
	a=${f%%:*}; b=${f#*:}
	sed "s/^mesh = .*/mesh = $a.msh/" "$SRC/box_p3.ofe" > "$WORK/$a.ofe"
	(cd "$WORK" && "$OFE" -n 2 "$a.ofe" > /dev/null && "$OFE_POST" > /dev/null)
	compare "C ($a box, 24 elements) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/p67a.csv"
	sed "s/^mesh = .*/mesh = $b.msh/" "$WORK/$a.ofe" > "$WORK/p67f.ofe"
	rm -f "$WORK/rlc.csv"
	(cd "$WORK" && "$OFE" -n 2 p67f.ofe > /dev/null && "$OFE_POST" > /dev/null)
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/p67f.csv" 2> /dev/null
	if [ -s "$WORK/p67f.csv" ] && cmp -s "$WORK/p67a.csv" "$WORK/p67f.csv"; then echo "  $b == ASCII 2.2 -> OK"
	else echo "  $b differs from ASCII 2.2 (or did not run) -> NG" >&2; status=1; fi
done
for p in 6:i:sextic 7:j:septic; do
	o=${p%%:*}; r=${p#*:}; c=${r%%:*}; w=${r#*:}
	run_case "nodal_test_p$o"
	if grep -q "($c) $w" "$WORK/ofe.log" && grep -q "centroid gradient ($w" "$WORK/ofe.log" && grep -q "result : passed" "$WORK/ofe.log"; then
		echo "  P$o self test ($w field identity, machine precision) : $(awk -v c="($c)" '$1 == c { print $NF; exit }' "$WORK/ofe.log") -> OK"
	else
		echo "  P$o self test : the $w-field identity was not run or failed -> NG" >&2
		status=1
	fi
done
set -- $(p4c coax_p6)
res=$(awk -v c="$1" -v v="$2" -v c4="$c4" -v v4="$v4" 'function ab(x) { return (x < 0) ? -x : x }
	BEGIN { ce = 2 * atan2(0, -1) * 8.8541878128e-12 * 2.1 / log(3); ve = atan2(0, -1) * (1.5e-3 ^ 2 - 0.5e-3 ^ 2) * 0.1e-3
	e = (c - ce) / ce; w = (v - ve) / ve; e4 = (c4 - ce) / ce; w4 = (v4 - ve) / ve
	ok = (c != "") && (ab(e) <= 1e-5) && (ab(w) <= 2e-6) && (ab(e) < ab(e4)) && (ab(w) < ab(w4))
	printf "%s (P6 h = 0.8mm : C %+.2e, volume %+.2e; P4 h = 0.4mm : C %+.2e, volume %+.2e)", ok ? "OK" : "NG", e, w, e4, w4 }')
echo "  curved coax, order 6 at twice the element size beats order 4 (geometry error h^8) : $res"
case "$res" in NG*) status=1 ;; esac
set -- $(awk 'BEGIN { pi = atan2(0, -1); mu = 4e-7 * pi; s = 5.8e7; W = 1e-3; t = 0.05e-3; d = 0.2e-3
	om = 2 * pi * 4e8; de = sqrt(2 / (om * mu * s)); a = t / de
	S = (exp(2 * a) - exp(-2 * a)) / 2; C = (exp(2 * a) + exp(-2 * a)) / 2
	k = 1 / (de * s * W * (C - cos(2 * a)))
	printf "%.12e", 2 * k * (S + sin(2 * a)) }')
r400=$1
for o in 5 6 7; do
	sed -e "s/plate2d.msh/plate2d_p$o.msh/" -e "s/^frequency = .*/frequency = 4e8/" "$SRC/plate2d_ac.ofe" > "$WORK/p67ac.ofe"
	rm -f "$WORK/rlc.csv"
	(cd "$WORK" && "$OFE" -n 2 p67ac.ofe > /dev/null && "$OFE_POST" > /dev/null)
	e=$(awk -v a="$(value_of Rf)" -v r="$r400" 'BEGIN { if (a == "") { print "nan"; exit } d = (a - r) / r; if (d < 0) d = -d; printf "%.3e", d }')
	eval "r4h_$o=\$e"
done
res=$(awk -v a="$r4h_5" -v b="$r4h_6" -v c="$r4h_7" 'BEGIN {
	ok = (b <= 2e-6) && (c <= 2e-7) && (a > 10 * b) && (b > 5 * c)
	printf "%s (R error : P5 %.1e, P6 %.1e, P7 %.1e)", ok ? "OK" : "NG", a, b, c }')
echo "  2-D, 400 MHz skin effect : order 6 within 2e-6 and 10x better than order 5, order 7 within 2e-7 and 5x better than order 6 : $res"
case "$res" in NG*) status=1 ;; esac
awk 'BEGIN { e = 0 } /^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 } e && ($2 == 42) { s = $1 " 25"; for (i = 3; i <= 3 + $3 + 21; i++) s = s " " $i; print s; next } { print }' "$WORK/box_p6.msh" > "$WORK/box_p6_tri21.msh"
sed 's/^mesh = .*/mesh = box_p6_tri21.msh/' "$WORK/box_p6.ofe" > "$WORK/p6x_tri.ofe"
th_reject "84-node tetrahedra with 21-node electrode triangles" p6x_tri.ofe "21-node (order 5) triangles are supported only as the faces of 56-node"
#   analysis E (辺要素の自己検証) は [A order 3+] の (h) で通す

# ---- 次数の違う四面体 / 三角形が面を共有する格子 (低い次数を最高次数 P に上げる、unstruct.c の ho_raise) ----
#
# 辺・面の節点は「全体の頂点番号と多重指数の組」を鍵に隣と共有し、新しい節点は元の要素の q 次の写像の上に置く。
# 格子は gmsh 製の一様な高次の格子の一部を mesh_downgrade.awk で低い次数に落としたもの (使われない節点は残る)。
#   (a) 箱 (z の下半分を落とす : 界面が場を横切るので、界面が裂けると C が変わる) : 1/2・1/3・2/4 の 3 通りで
#       C が閉形式に全桁、analysis = P の P 次の場の恒等式。**P の恒等式は界面の裂け目を見ない** (要素毎の和。鍵の
#       頂点番号を並べ替えない誤りでも通る、実測) ので C と (c) が要る
#   (b) 曲がった同軸 (coax_p4 の x < 0 を 2 次に落とす) : 新しい節点は 2 次の写像の上 (実測 −1.7e-4。同じ頂点の
#       一様な 2 次 +1.4e-4、1 次 −4.1e-3)。許容 5e-4。頂点だけから置く誤りは格子が裏返って弾かれる (実測)
#   (c) 断面 2 次元 (plate2d_p3 の y < 0 を 1 次に落とす) : 直線の三角形なので上げた格子は一様な 3 次と同じ離散問題に
#       なり、400 MHz の rlc.csv が一様な 3 次とバイト一致
#   (d) 三角形が四面体の最高次数より高い、ピラミッドのある格子、は弾く
#   (e) 六面体・角柱 (z の下半分を 1 次に落とす) : 六面体 1 -> 2 (27 節点)・1 -> 3・1 -> 4、角柱 1 -> 3 で C が全桁、P 次の場の
#       恒等式。電極の四角形・三角形は要素の面の節点に置き直す (**目標と同じ次数の面も**。面の下の要素が上げられると元の面の
#       節点は要素のものでない。実測 : 曲がった同軸の 9 節点の電極を残したら C' −16%)。種別の違う要素の次数の混在は [thermal
#       mixed] の同軸 (角柱 + 六面体) でも見る
echo "[mixed orders] elements of different orders sharing faces (raised to the highest order)"
cp "$SRC/mesh_downgrade.awk" "$SRC/box_p2.msh" "$SRC/box_p3.msh" "$SRC/box_p4.msh" "$SRC/coax_p4.msh" "$SRC/plate2d_p3.msh" "$WORK/"
mx_down() {	# mx_down <入力> <出力> <列> <しきい値> <四面体の型> <三角形の型> <1|2>
	awk -v col="$3" -v xc="$4" -v tf="$5" -v ff="$6" -v to="$7" -f "$WORK/mesh_downgrade.awk" "$WORK/$1" > "$WORK/$2"
}
for c in box_p2:11:9:1:1:2:c:quadratic box_p3:29:21:1:1:3:f:cubic box_p4:30:23:2:2:4:g:quartic; do
	IFS=: read -r m tf ff to q p lt fw <<EOF2
$c
EOF2
	mx_down "$m.msh" mx_$m.msh 4 1e-4 "$tf" "$ff" "$to"
	sed -e "s/^mesh = .*/mesh = mx_$m.msh/" "$SRC/box_p3.ofe" > "$WORK/mx.ofe"
	(cd "$WORK" && "$OFE" -n 2 mx.ofe > /dev/null && "$OFE_POST" > /dev/null) || { echo "  mx_$m : the run failed -> NG" >&2; status=1; }
	if grep -q "Mixed orders : [0-9]* tetrahedra of order $q and up raised to order $p" "$WORK/ofe.log"; then
		echo "  mx_$m : $(grep 'Mixed orders' "$WORK/ofe.log" | sed 's/^Mixed orders : //') -> OK"
	else echo "  mx_$m : the order-$q tetrahedra were not raised to order $p -> NG" >&2; status=1; fi
	compare "mx_$m (orders $q and $p, interface across the field) : C (parallel plate) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	sed -e "s/^analysis = .*/analysis = P/" "$WORK/mx.ofe" > "$WORK/mxp.ofe"
	(cd "$WORK" && "$OFE" -n 2 mxp.ofe > /dev/null 2>&1)
	if grep -q "($lt) $fw" "$WORK/ofe.log" && grep -q "result : passed" "$WORK/ofe.log"; then
		echo "  mx_$m (analysis P) : order-$p field identity on the raised mesh ($(awk -v c="($lt)" '$1 == c { print $NF; exit }' "$WORK/ofe.log")) -> OK"
	else echo "  mx_$m (analysis P) : the order-$p field identity failed or did not run -> NG" >&2; status=1; fi
done
mx_down coax_p4.msh mx_coax24.msh 2 0 30 23 2
sed "s/^mesh = .*/mesh = mx_coax24.msh/" "$SRC/coax_p2.ofe" > "$WORK/mxc.ofe"
rm -f "$WORK/rlc.csv"
(cd "$WORK" && "$OFE" -n 2 mxc.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  mx_coax24 : the run failed -> NG" >&2; status=1; }
compare "curved coax, half order 2 raised to order 4 : C' [F/m]" "$(value_of C)" 1.063416612e-10 5e-4
mx_down plate2d_p3.msh mx_plate13.msh 3 0 -1 21 1
for m in plate2d_p3 mx_plate13; do
	sed -e "s/plate2d.msh/$m.msh/" -e "s/^frequency = .*/frequency = 4e8/" "$SRC/plate2d_ac.ofe" > "$WORK/mxa.ofe"
	rm -f "$WORK/rlc.csv"
	(cd "$WORK" && "$OFE" -n 2 mxa.ofe > /dev/null && "$OFE_POST" > /dev/null) || { echo "  $m : the run failed -> NG" >&2; status=1; }
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/$m.mx.csv" 2> /dev/null
done
if [ -s "$WORK/plate2d_p3.mx.csv" ] && cmp -s "$WORK/plate2d_p3.mx.csv" "$WORK/mx_plate13.mx.csv"; then
	echo "  2-D, triangles of orders 1 and 3 (raised to 3) == uniform order 3 : rlc.csv identical at 400 MHz -> OK"
else echo "  2-D mixed-order plate differs from the uniform order-3 plate -> NG" >&2; status=1; fi
mx_down box_p4.msh mx_t1.msh 4 1e-4 30 -1 2
mx_down mx_t1.msh mx_tri.msh 4 1 30 -1 1
sed "s/^mesh = .*/mesh = mx_tri.msh/" "$SRC/box_p3.ofe" > "$WORK/mxr1.ofe"
th_reject "mixed-order tetrahedra (1 and 2) with 15-node electrode triangles" mxr1.ofe "a triangle is order 4 but the highest tetrahedron order is 2"
cp "$SRC/box_hexpyrtet.msh" "$SRC/box_hexpyrtet.ofe" "$WORK/"
awk '/^\$Nodes/ { print; getline; n = $1; print n + 10; next }
	/^\$EndNodes/ { for (i = 0; i < 10; i++) printf "%d %.6g %.6g %.6g\n", n + 1 + i, 5e-3 + (i == 1) * 1e-4, 5e-3 + (i == 2) * 1e-4, (i == 3) * 1e-4; print; next }
	/^\$Elements/ { print; getline; m = $1; print m + 1; next }
	/^\$EndElements/ { s = (m + 1) " 11 2 1 1"; for (i = 1; i <= 10; i++) s = s " " (n + i); print s; print; next }
	{ print }' "$WORK/box_hexpyrtet.msh" > "$WORK/mx_pyr.msh"
sed "s/^mesh = .*/mesh = mx_pyr.msh/" "$WORK/box_hexpyrtet.ofe" > "$WORK/mxr2.ofe"
th_reject "elements of different orders with pyramids" mxr2.ofe "pyramids are first order only"
for c in box_hex2:12:-1:10:hex:2:quadratic box_hex3:92:-1:36:hex:3:cubic box_hex4:93:-1:37:hex:4:quartic box_prism3:90:21:36:prism:3:cubic; do
	IFS=: read -r m tf ff qf kd p fw <<EOF2
$c
EOF2
	cp "$SRC/$m.msh" "$WORK/"
	awk -v col=4 -v xc=1e-4 -v tf="$tf" -v ff="$ff" -v qf="$qf" -v kind="$kd" -v to=1 -f "$WORK/mesh_downgrade.awk" "$WORK/$m.msh" > "$WORK/mx_$m.msh"
	sed -e "s/box_hex2.msh/mx_$m.msh/" -e "s/^analysis.*/analysis = C/" "$SRC/box_hex2.ofe" > "$WORK/mxh.ofe"
	rm -f "$WORK/rlc.csv"
	(cd "$WORK" && "$OFE" -n 2 mxh.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  mx_$m : the run failed -> NG" >&2; status=1; }
	grep -q "^Mixed orders : [0-9]* elements of order 1 and up raised to order $p" "$WORK/ofe.log" || { echo "  mx_$m : not raised to order $p -> NG" >&2; status=1; }
	compare "mx_$m ($kd, orders 1 and $p, interface across the field) : C (parallel plate) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	sed -e "s/^analysis.*/analysis = P/" "$WORK/mxh.ofe" > "$WORK/mxhp.ofe"
	(cd "$WORK" && "$OFE" -n 2 mxhp.ofe > /dev/null 2>&1)
	if grep -q '^\*\*\*' "$WORK/ofe.log" || ! grep -q "$fw field : phi" "$WORK/ofe.log"; then
		echo "  mx_$m (analysis P) : the order-$p field identity failed or did not run -> NG" >&2; status=1
	else echo "  mx_$m (analysis P) : order-$p field identity on the raised mesh ($(awk -v w="$fw" '$1 == w && $2 == "field" { print $NF; exit }' "$WORK/ofe.log")) -> OK"; fi
done

# ---- 熱連成の六面体・角柱・ピラミッド (1 次、混在可) ----
#
# 剛性・発熱・質量・要素の温度を 1 つの積分点 API (unstruct.c の elem3d_ip) で組む
# (thermal.c の gen_assemble の注釈)。六面体は構造格子と同じ 2x2x2 Gauss なので、
# **同じ格子を六面体で書いたものは構造格子と同じ離散問題**になり全桁一致する
# (閉形式も許容誤差も要らない恒等式)。角柱・ピラミッドは閉形式・誤差の次数・
# エネルギーの恒等式・露出面の数で見る。
echo "[thermal 3-D elements] hexahedra / prisms / pyramids (mixed) in the electro-thermal coupling"
for m in thermal_hex40 thermal_pri20 thermal_pri40 thermal_pyr80 thermal_pyr160 \
	thermal_neck_hex thermal_neck_pri thermal_neck_pyr thermal_neck_mix; do
	cp "$SRC/$m.msh" "$WORK/"
done
cp "$SRC/thermal_hex.ofe" "$SRC/thermal_neck_mix.ofe" "$SRC/thermal_bar.ofe" "$SRC/thermal_conv.ofe" \
	"$SRC/thermal_rad.ofe" "$SRC/thermal_transient.ofe" "$SRC/thermal_neck.ofe" "$WORK/"
t3_mesh() {	# t3_mesh <格子> <入力> <出力> : mesh の行を差し替える
	awk -v m="$1" '/^mesh/{print "mesh = " m ".msh"; next} {print}' "$WORK/$2" > "$WORK/$3"
}
t3_tm() { awk '/^  T_max = / { print $(NF-1); exit }' "$WORK/ofe.log"; }
t3_same() {	# t3_same <ラベル> <構造格子の入力> <六面体の入力> : R / T_max / T_mean が全桁一致
	th_run "$2"; r0=$(th_r); t0=$(th_t); m0=$(t3_tm)
	th_run "$3"
	compare "$1 : hexahedra == structured grid : R [ohm]" "$(th_r)" "$r0" 1e-10
	compare "$1 : hexahedra == structured grid : T_max [degC]" "$(th_t)" "$t0" 1e-10
	compare "$1 : hexahedra == structured grid : T_mean [degC]" "$(t3_tm)" "$m0" 1e-10
}
t3_iter() {	# t3_iter <ラベル> <上限> : Newton の反復回数
	nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
	if [ -n "$nit" ] && [ "$nit" -le "$2" ]; then
		echo "  $1 : converged in $nit iterations (<= $2) -> OK"
	else
		echo "  $1 : ${nit:-?} iterations (> $2) -> NG" >&2
		status=1
	fi
}
# (1) 六面体 == 構造格子 (同じ 40 x 1 x 1 の不等間隔の棒) : 両端ヒートシンク、側面の
#     対流だけ、ヒートシンク + 輻射、過渡 (CN)。対流・輻射は四角形の露出面 (2x2 Gauss)
#     を通る。不等間隔なのは要素の体積を揃えないため (T_mean を体積でなく要素の数で
#     平均する誤りは等間隔の格子では等価)
t3_grid() {	# t3_grid <入力> <出力> : 構造格子の x を thermal_hex40.msh と同じ不等間隔に
	awk '/^xmesh/{print "xmesh = 0 10 4e-3 30 0.01"; next} {print}' "$WORK/$1" > "$WORK/$2"
}
t3_grid thermal_bar.ofe t3g_bar.ofe
t3_same "bar, heat sinks" t3g_bar.ofe thermal_hex.ofe
th_energy "hex : thermal_hex"
awk '{ gsub(/0\.5e-3/, "1e-3"); print }' "$SRC/thermal_conv.ofe" > "$WORK/t3g_conv0.ofe"
t3_grid t3g_conv0.ofe t3g_conv.ofe
awk '/^heatsink/{next} /^thermal/{print "convection = 200 30"; print "thermal = 6"; next} {print}' \
	"$SRC/thermal_hex.ofe" > "$WORK/t3h_conv.ofe"
t3_same "bar, convection only" t3g_conv.ofe t3h_conv.ofe
th_nface "hex : convection, sides only (not the electrode faces)" 160
th_energy "hex : convection only"
awk '/^kthermal/{print "kthermal = 2 1"; next} /^thermal/{print "radiation = 1 10"; print "thermal = 2.5"; next} {print}' \
	"$SRC/thermal_hex.ofe" > "$WORK/t3h_rad.ofe"
t3_grid thermal_rad.ofe t3g_rad.ofe
t3_same "bar, heat sinks + radiation" t3g_rad.ofe t3h_rad.ofe
th_energy "hex : heat sinks + radiation"
awk '/^tempco/{print "tempco = 2 0 0"; next} /^temperature/{print "temperature = 20"; next}
	/^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 2 200 0.5"; next} {print}' \
	"$SRC/thermal_hex.ofe" > "$WORK/t3h_tr.ofe"
t3_grid thermal_transient.ofe t3g_tr.ofe
t3_same "bar, transient (CN)" t3g_tr.ofe t3h_tr.ofe
th_balance "hex : transient (CN)"
# (2) 首のある導体 (∇φ が要素内で変わる) を 2 材料で : 六面体は構造格子と全桁一致、
#     角柱・ピラミッド・混在 (六面体 + 首のピラミッド) はエネルギーの恒等式
#     (発熱を剛性と同じ積分点で組んだことの検査。1 点で取ると落ちる)
awk '/^material/{print; print "material = 1.0 2e6"; next}
	/^kthermal/{print; print "kthermal = 3 25"; print "tempco = 3 2e-3 0"; next}
	/^geometry = 2 1 1e-3/{print "geometry = 3 1 1e-3 2e-3 0.75e-3 1.25e-3 0 0.5e-3"; next}
	/^fieldout/{next} {print}' "$SRC/thermal_neck.ofe" > "$WORK/t3g_neck.ofe"
t3_mesh thermal_neck_hex thermal_neck_mix.ofe t3h_neck.ofe
t3_same "neck, two materials" t3g_neck.ofe t3h_neck.ofe
for k in pri pyr mix; do
	t3_mesh "thermal_neck_$k" thermal_neck_mix.ofe "t3_neck_$k.ofe"
	th_run "t3_neck_$k.ofe"
	th_energy "neck, two materials, $k"
done
#     場の出力の温度の最大 = ログの T_max (混在格子のセル)
awk '/^thermal/{print "fieldout = 1"} {print}' "$SRC/thermal_neck_mix.ofe" > "$WORK/t3_fo.ofe"
th_run t3_fo.ofe
tvtk=$(awk '/^SCALARS T / { f = 1; next } f && /^LOOKUP/ { next }
	f && (NF == 1) { if (($1 + 0) > m) m = $1 + 0; next } f { f = 0 }
	END { printf "%.9e", m }' "$WORK/ofe_field.vtk")
compare "mixed : max of the T field == T_max in the log [degC]" "$tvtk" "$(th_t)" 1e-8
# (3) 角柱・ピラミッドの棒 : 閉形式と h² の誤差比。α = 0 では φ が 1 次で厳密
#     (ピラミッドの求積の重み 1/2 を落とすと R が倍になる)。ピラミッドは h⁴ の項が
#     大きく、40 -> 80 で比 4.53、80 -> 160 で 4.15 と上から 4 に近づく (実測) ので
#     80 -> 160 で見て上限を 4.3 にする
set -- $(th_exact 7.5 4e-3)
rex=$1; tex=$2
t3_mesh thermal_pri40 thermal_hex.ofe t3p40.ofe
th_run t3p40.ofe
compare "prism : R(I = 7.5 A, mL = 1.5) [ohm]" "$(th_r)" "$rex" 0.001
compare "prism : T_max(I = 7.5 A) [degC]" "$(th_t)" "$tex" 0.001
th_energy "prism : bar"
r40=$(th_r)
t3_mesh thermal_pri20 thermal_hex.ofe t3p20.ofe
th_run t3p20.ofe
tt_order "prism : second-order convergence of R" "$(th_r)" "$r40" "$rex"
t3_mesh thermal_pyr160 thermal_hex.ofe t3y160.ofe
th_run t3y160.ofe
r160=$(th_r)
compare "pyramid : R(I = 7.5 A, mL = 1.5) [ohm]" "$r160" "$rex" 0.001
compare "pyramid : T_max(I = 7.5 A) [degC]" "$(th_t)" "$tex" 0.001
th_energy "pyramid : bar"
t3_mesh thermal_pyr80 thermal_hex.ofe t3y80.ofe
th_run t3y80.ofe
res=$(awk -v a="$(th_r)" -v b="$r160" -v e="$rex" 'BEGIN {
	e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 3.8) && (q < 4.3)) ? "OK" : "NG", e1, e2, q }')
echo "  pyramid : second-order convergence of R (80 -> 160 cells) : $res"
case "$res" in NG*) status=1 ;; esac
set -- $(th_exact 7.5 0)
for k in pri40 pyr80; do
	t3_mesh "thermal_$k" thermal_hex.ofe t3_a0.ofe
	awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$WORK/t3_a0.ofe" > "$WORK/t3_a0b.ofe"
	th_run t3_a0b.ofe
	compare "$k : R(alpha = 0) == the uncoupled value [ohm]" "$(th_r)" "$1" 1e-10
done
# (4) 側面の対流 : 露出面の数 (角柱は y 側面の四角形 2 x 40 + z 側面の三角形 2 x 2 x 40
#     = 240、ピラミッドは 4 側面の四角形 4 x 80 = 320) とエネルギーの恒等式。
#     輻射の Newton (α = 0) は四角形・三角形の面の接線の検査 (反復回数の上限)
for k in pri40:240 pyr80:320; do
	m=${k%%:*}; n=${k##*:}
	t3_mesh "thermal_$m" t3h_conv.ofe t3_c.ofe
	th_run t3_c.ofe
	th_nface "$m : convection, sides only" "$n"
	th_energy "$m : convection only"
	t3_mesh "thermal_$m" t3h_rad.ofe t3_r.ofe
	awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$WORK/t3_r.ofe" > "$WORK/t3_r0.ofe"
	th_run t3_r0.ofe
	th_energy "$m : heat sinks + radiation, alpha = 0"
	t3_iter "$m : Newton on radiation (alpha = 0)" 5
done
# (5) 断熱の棒の過渡 (CN) : 発熱と質量の行和が同じ節点の重み ∫N_i なので温度は一様の
#     まま進み、θ 法の漸化式に機械精度で一致する (質量を積分点で組んだことの検査)
for m in pri40 pyr80; do
	awk -v m="$m" '/^mesh/{print "mesh = thermal_" m ".msh"; next} /^heatsink/{next}
		/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 10 100 0.5"} {print}' \
		"$SRC/thermal_hex.ofe" > "$WORK/t3_ad.ofe"
	th_run t3_ad.ofe
	set -- $(th_adia 100 0.5)
	compare "$m : adiabatic bar, CN : R == the discrete recurrence [ohm]" "$(th_r)" "$1" 1e-9
	compare "$m : adiabatic bar, CN : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
	th_balance "$m : adiabatic bar, CN"
done
# (6) 2 次の六面体・角柱は [thermal second-order hex/prism] で解く (以前はここで弾いていた)
cp "$SRC/bar_hex.msh" "$WORK/"

# ---- 熱連成の高次の四面体 (2〜4 次、等パラメトリック、曲がった面) ----
#
# 剛性・発熱・質量・要素の温度は積分点の経路 (elem3d_ip、Duffy + 各方向 p + 2 点)、露出面は
# 親の四面体の形状関数を面に制限して積分する (thermal.c の surf_tri_ho)。σ(T) は要素内で一定
# (要素の平均温度) なので、温度の解く精度が上がっても R は h² (定数は 1 次の 1/10)。高次の
# 効き目は曲がった形 (面積・体積の幾何の誤差) に出る。
echo "[thermal higher-order tet] electro-thermal coupling on 10/20/35-node tetrahedra (curved)"
for m in thermal_tet20_p2 thermal_tet40_p2 coax_p3 coax_p4 cyl_p1 cyl_p2 cyl_p3; do cp "$SRC/$m.msh" "$WORK/"; done
# (1) 直線の棒 (2 次) : 閉形式と誤差比。σ が要素内で一定なので h² (20 -> 40 で比 4.2、
#     40 分割で R -1.9e-5、1 次の 1/10。実測)
set -- $(th_exact 7.5 4e-3)
rex=$1; tex=$2
for nx in 20 40; do
	awk -v m="thermal_tet${nx}_p2" '/^mesh/{print "mesh = " m ".msh"; next} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tho_bar.ofe"
	th_run tho_bar.ofe
	eval "r$nx=\$(th_r)"
done
compare "P2 bar : R(I = 7.5 A, mL = 1.5) [ohm]" "$r40" "$rex" 1e-4
th_energy "P2 bar"
tt_order "P2 bar : second-order convergence of R (sigma constant per element)" "$r20" "$r40" "$rex"
#     断熱の棒の過渡 (CN) : 2 次の頂点の ∫N は負 (-V/20) だが、発熱と質量の行和が同じ重みなので
#     温度は一様のまま進み、θ 法の漸化式に機械精度で一致する
awk '/^mesh/{print "mesh = thermal_tet40_p2.msh"; next} /^heatsink/{next}
	/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 10 100 0.5"} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tho_ad.ofe"
th_run tho_ad.ofe
set -- $(th_adia 100 0.5)
compare "P2 bar : adiabatic, CN : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
th_balance "P2 bar : adiabatic, CN"
# (2) 曲がった同軸、内外の円筒が電極 (半径方向の電流、q ∝ 1/r²)、両方ヒートシンク、α = 0 :
#       R = ln(b/a) / (2πσ l_z)、T(r) = T_s + (c/2)(ln r - ln a)(ln b - ln r)、c = (I/(2π l_z))²/(σk)
#     R の誤差は静電の C' と同じ幾何の誤差 (3 次 h⁴ : 3.3e-4、4 次 h⁶ : 8.5e-6)。体積平均の
#     温度上昇は 3 次 6.3e-4 / 4 次 2.5e-5 (実測)。許容は C' と同じ (3 次 5e-4、4 次 1e-4) と、
#     温度にその 2 倍
set -- $(awk 'BEGIN { pi = atan2(0, -1); a = 0.5e-3; b = 1.5e-3; lz = 0.1e-3; s = 1e6; k = 10; I = 30
	c = (I / (2 * pi * lz)) ^ 2 / (s * k); sa = log(a); sb = log(b); n = 20000; num = 0; den = 0
	for (i = 0; i < n; i++) { r = a + (b - a) * (i + 0.5) / n; x = log(r); num += r * (c / 2) * (x - sa) * (sb - x); den += r }
	printf "%.12e %.12e", log(b / a) / (2 * pi * s * lz), num / den }')
rcx=$1; tcx=$2
for k in 3:5e-4:1e-3 4:1e-4:2e-4; do
	o=${k%%:*}; rest=${k#*:}; rt=${rest%%:*}; tt=${rest#*:}
	printf 'OpenFEM 1 0\ntitle = radial conduction heating in a coax\nmesh = coax_p%s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nheatsink = 0 20\nheatsink = 1 20\ntemperature = 20\nthermal = 30\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$o" > "$WORK/tho_cx.ofe"
	th_run tho_cx.ofe
	compare "P$o curved coax, radial current : R = ln(b/a) / (2 pi sigma l) [ohm]" "$(th_r)" "$rcx" "$rt"
	compare "P$o curved coax, radial current : volume-mean temperature rise [K]" \
		"$(awk -v t="$(t3_tm)" 'BEGIN { printf "%.12e", t - 20 }')" "$tcx" "$tt"
	th_energy "P$o curved coax"
done
# (3) 円柱 (cyl_p*.msh、gmsh 製、半径 0.5mm) を曲がった側面だけから冷やす (端面は電極で
#     断熱)、α = 0 : T(r) = T_s + q (R² - r²)/(4k)、T_s - T_a = P / (h 2πRL)、
#     T_mean = T_s + q R²/(8k)。**曲がった面の積分を検査できるのはこれ** (面積の誤差が
#     T_s に直接入る)。1 次は折れ線の面で +3.1%、2 次・3 次は h⁴ で 1e-4 (実測)。許容 5e-4、
#     1 次との比 (> 50) も見る
set -- $(awk 'BEGIN { pi = atan2(0, -1); R0 = 0.5e-3; L = 0.5e-3; s = 1e6; k = 10; h = 2000; I = 5
	re = L / (s * pi * R0 * R0); P = I * I * re; ts = P / (h * 2 * pi * R0 * L); q = P / (pi * R0 * R0 * L)
	printf "%.12e %.12e", re, ts + q * R0 * R0 / (8 * k) }')
rcy=$1; tcy=$2
for o in 1 2 3; do
	printf 'OpenFEM 1 0\ntitle = cylinder cooled from its curved side\nmesh = cyl_p%s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 2000 20\ntemperature = 20\nthermal = 5\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$o" > "$WORK/tho_cy.ofe"
	th_run tho_cy.ofe
	eval "ecy$o=\$(awk -v t=\"\$(t3_tm)\" -v e=\"\$tcy\" 'BEGIN { d = (t - 20 - e) / e; if (d < 0) d = -d; printf \"%.6e\", d }')"
	th_nface "P$o cylinder : convection on the curved side only" 80
	th_energy "P$o cylinder"
	if [ "$o" != 1 ]; then
		compare "P$o cylinder : R = L / (sigma pi R^2) [ohm]" "$(th_r)" "$rcy" 5e-4
		compare "P$o cylinder, curved-side convection : volume-mean temperature rise [K]" \
			"$(awk -v t="$(t3_tm)" 'BEGIN { printf "%.12e", t - 20 }')" "$tcy" 5e-4
	fi
done
res=$(awk -v a="$ecy1" -v b="$ecy2" 'BEGIN { printf "%s (P1 %.2e vs P2 %.2e)", (a > 50 * b) ? "OK" : "NG", a, b }')
echo "  cylinder : the curved P2 surface beats the P1 facets (geometry h^4 vs h^2) : $res"
case "$res" in NG*) status=1 ;; esac
# (高次の四面体の相互放射は [radiation exchange] の (8))

# ---- 熱連成の 2 次の六面体・角柱 (20 / 27 / 15 / 18 節点、等パラメトリック) ----
#
# 剛性・発熱・質量・要素の温度は積分点の経路 (elem3d_ip、六面体 3x3x3 Gauss、角柱 面内 Duffy 9 点
# x ζ 3 点 = 剛性の組み立てと同じ則)、露出面は親の要素の形状関数を面に制限して積分する
# (elem3d_face_ip : 六面体・角柱の四角形 3 x 3 Gauss、角柱の三角形 Duffy 9 点。thermal.c の
# surf_face_ho)。
echo "[thermal second-order hex/prism] electro-thermal coupling on 20/27-node hexahedra and 15/18-node prisms"
for m in box_hex2 box_hex2s box_prism2 box_prism2s thermal_slab_h27 thermal_slab_p18 \
	thermal_ann_h27 thermal_ann_p18 thermal_ann_h8 thermal_ann_h27b; do cp "$SRC/$m.msh" "$WORK/"; done
h2_write() {	# h2_write <格子> <電流> <追加の行 (\n 区切り)> : 熱連成 (α = 0) の入力を書く
	printf 'OpenFEM 1 0\ntitle = second-order hex/prism thermal\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\n%b\ntemperature = 20\nthermal = %s\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' \
		"$1" "$3" "$2" > "$WORK/th2.ofe"
}
h2_ofe() { h2_write "$@"; th_run th2.ofe; }
# (1) 箱 (1mm x 1mm x 0.2mm、厚さ方向の電流、両方の電極がヒートシンク)、α = 0 :
#     T(z) = T_s + q z (t - z)/(2k) は 2 次なので 4 種とも (serendipity を含めて) 厳密に表せる
#     R = t/(σA) = 2e-4、T_max = 20 + q t²/(8k) = 25、T_mean = 20 + q t²/(12k) (q = 1e10)
for m in box_hex2 box_hex2s box_prism2 box_prism2s; do
	h2_ofe "$m" 100 'heatsink = 0 20\nheatsink = 1 20'
	compare "$m : R = t / (sigma A) [ohm]" "$(th_r)" 2e-4 1e-9
	compare "$m : T_max = T_s + q t^2 / (8k) [degC]" "$(th_t)" 25 1e-9
	compare "$m : T_mean = T_s + q t^2 / (12k) [degC]" "$(t3_tm)" 23.3333333333 1e-9
	th_energy "$m"
done
#     σ(T) あり (α = 4e-3) : 要素内一定の σ なので厳密でない。恒等式だけ
h2_write box_hex2 100 'heatsink = 0 20\nheatsink = 1 20'
awk '/^tempco/{print "tempco = 2 4e-3 20"; next} {print}' "$WORK/th2.ofe" > "$WORK/th2a.ofe"
th_run th2a.ofe
th_energy "box_hex2, alpha = 4e-3"
#     断熱の過渡 (CN、σ(T) あり) : 一様な温度のまま進むので thermal_bar と同じ漸化式
#     (断面 1mm² で J が同じ。漸化式は長さに依らない) に機械精度で一致する
for m in box_hex2 box_prism2s; do
	awk -v m="$m" '/^mesh/{print "mesh = " m ".msh"; next} /^heatsink/{next}
		/^thermal/{print "heatcap = 2 3.4e6"; print "transient = 10 100 0.5"} {print}' "$SRC/thermal_tet.ofe" > "$WORK/th2ad.ofe"
	th_run th2ad.ofe
	set -- $(th_adia 100 0.5)
	compare "$m : adiabatic, CN : T == the discrete recurrence [degC]" "$(th_t)" "$2" 1e-9
	th_balance "$m : adiabatic, CN"
done
# (2) 板 (x 向きの電流、z = 0 の面だけを直方体で選んで対流、h = 5000、z = t は断熱)、α = 0 :
#     T(z) = T_a + q t/h + q (t z - z²/2)/k が 2 次なので機械精度。六面体は四角形の面、
#     角柱は三角形の面を通る。T_max = T_a + q t/h + q t²/(2k) = 125、
#     T_mean = T_a + q t/h + q t²/(3k) = 123.333... (q = 2.5e9)。
#     **両面を冷やすと検査にならない**: 面の形状関数を ζ について裏返す (z = 0 の面の熱を
#     中央面の節点で取る) 誤りは、半分ずつが断熱の板 (厚さ t/2) になるだけで T_max・T_mean が
#     全桁一致する (実測)。片面だけ冷やして対称性を壊す
for m in thermal_slab_h27:8 thermal_slab_p18:16; do
	g=${m%%:*}; nf=${m#*:}
	h2_ofe "$g" 10 'convection = 5000 20 0 2e-3 0 1e-3 0 0'
	compare "$g : R = L / (sigma A) [ohm]" "$(th_r)" 1e-2 1e-9
	compare "$g : bottom face cooled : T_max (parabola) [degC]" "$(th_t)" 125 1e-9
	compare "$g : bottom face cooled : T_mean [degC]" "$(t3_tm)" 123.333333333 1e-9
	th_nface "$g : convection on the bottom face only" "$nf"
	th_energy "$g"
	#   輻射 (α = 0、0.5 A) : 面の Newton の接線・右辺・熱の評価点がそろっていることの恒等式
	h2_ofe "$g" 0.5 'radiation = 0.8 20 0 2e-3 0 1e-3 0 0'
	th_energy "$g : radiation"
done
# (3) 円環 (内径 0.5mm、外径 1mm、軸方向の電流、内外の円筒面 = 曲がった面を対流 h = 2000)、α = 0 :
#     T(r) = T_a - q r²/(4k) + C1 ln r + C2 (両側の Robin 条件で C1, C2)。**曲がった面の積分を
#     検査できるのはこれ** (面積の誤差が温度に直接入る)。周方向 12 分割で 2 次 +2.3e-4、24 分割で
#     +1.3e-5 (h⁴、比 18)。1 次の折れ線 (同じ頂点) は +5.7%。許容は 2 次 5e-4、比 > 50
set -- $(awk 'BEGIN { pi = atan2(0, -1); a = 0.5e-3; b = 1e-3; lz = 0.2e-3; s = 1e6; k = 10; h = 2000; I = 10
	A = pi * (b * b - a * a); q = (I / A) ^ 2 / s
	Pb = -q * b * b / (4 * k); Pa = -q * a * a / (4 * k)
	a11 = -k / b - h * log(b); a12 = -h; r1 = h * Pb - q * b / 2
	a21 = k / a - h * log(a); a22 = -h; r2 = h * Pa + q * a / 2
	d = a11 * a22 - a12 * a21; c1 = (r1 * a22 - a12 * r2) / d; c2 = (a11 * r2 - a21 * r1) / d
	n = 20000; num = 0; den = 0
	for (i = 0; i < n; i++) { r = a + (b - a) * (i + 0.5) / n; num += r * (-q * r * r / (4 * k) + c1 * log(r) + c2); den += r }
	printf "%.12e %.12e", lz / (s * A), num / den }')
ran=$1; tan=$2
for m in thermal_ann_h8 thermal_ann_h27 thermal_ann_p18 thermal_ann_h27b; do
	h2_ofe "$m" 10 'convection = 2000 20'
	e=$(awk -v t="$(t3_tm)" -v x="$tan" 'BEGIN { d = (t - 20 - x) / x; if (d < 0) d = -d; printf "%.6e", d }')
	eval "e_$m=\$e"
	th_energy "$m"
	case "$m" in
	*_h8) ;;
	*) compare "$m : R = L / (sigma pi (b^2 - a^2)) [ohm]" "$(th_r)" "$ran" 5e-4
	   compare "$m : curved faces cooled : volume-mean temperature rise [K]" \
		"$(awk -v t="$(t3_tm)" 'BEGIN { printf "%.12e", t - 20 }')" "$tan" 5e-4 ;;
	esac
	case "$m" in *_h27b) th_nface "$m : convection on the curved faces only" 96 ;; *) th_nface "$m : convection on the curved faces only" 48 ;; esac
done
res=$(awk -v a="$e_thermal_ann_h8" -v b="$e_thermal_ann_h27" -v c="$e_thermal_ann_h27b" 'BEGIN {
	ok = (a > 50 * b) && (b > 10 * c)
	printf "%s (P1 %.2e, P2 %.2e, P2 refined %.2e)", ok ? "OK" : "NG", a, b, c }')
echo "  annulus : curved P2 faces beat the P1 facets, and converge at h^4 : $res"
case "$res" in NG*) status=1 ;; esac
# (高次の六面体・角柱の相互放射は [radiation exchange] の (8) の (c)〜(e))

# ---- 3 次の六面体 (64 節点、gmsh の hexa64) と角柱 (40 節点、prism40) ----
#
# 局所節点の並びは gmsh 4.15.2 の出力を分類した実測の表 (unstruct.c の HEX3_I / PRISM3_A)。形状関数は
# 次数によらない Lagrange の式 (1 次元の等間隔の積 / 三角形の多重指数)、積分は各方向 4 点 (角柱は
# 面内 Duffy 16 点 x ζ 4 点)。格子はすべて gmsh 製 (gmsh_box_ho.py / gmsh_coax_ho.py / gmsh_slab_ho.py)。
#   (a) 箱の平行平板 : C が閉形式に全桁、analysis = P の線形場・3 次の場の恒等式 (3 次特有の基底を
#       検査できるのは 3 次の場だけ)、2.2 / 4.1 / バイナリで rlc.csv が完全一致
#   (b) 曲がった同軸 : C' が h⁴ 以上で収束 (1/4 周 2 -> 4 分割で 六面体 24 倍・角柱 26 倍、実測)
#   (c) 熱連成 : 片面を冷やす板が閉形式に機械精度 (四角形の面 / 三角形の面)、曲がった面を冷やす円環が
#       h⁴ (1 次の折れ線の 1/1000)
echo "[order-3 hex/prism] 64-node hexahedra / 40-node prisms (gmsh hexa64 / prism40)"
for m in box_hex3 box_prism3; do
	for f in "" _41 _bin_22 _bin_41; do cp "$SRC/$m$f.msh" "$WORK/"; done
done
for m in coax_hex3a coax_hex3b coax_prism3a coax_prism3b slab_hex3 slab_prism3 ann_hex3a ann_hex3b ann_prism3a ann_prism3b ann_hex1b; do
	cp "$SRC/$m.msh" "$WORK/"
done
h3_c() {	# h3_c <格子> <解析> : box_hex2.ofe の格子と解析を差し替えて解く
	sed -e "s/box_hex2.msh/$1.msh/" -e "s/^analysis.*/analysis = $2/" "$SRC/box_hex2.ofe" > "$WORK/h3.ofe"
	(cd "$WORK" && "$OFE" -n 2 h3.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  $1 : the run failed -> NG" >&2; status=1; }
}
h3_selftest() {	# h3_selftest <ラベル> <cubic|skip> : analysis = P の自己検証の結果
	if grep -q '^\*\*\*' "$WORK/ofe.log"; then echo "  $1 : self test failed -> NG" >&2; status=1; return; fi
	if ! grep -q "linear field : .*rel. error" "$WORK/ofe.log"; then echo "  $1 : no linear-field line -> NG" >&2; status=1; return; fi
	if [ "$2" = cubic ]; then
		e=$(awk '/cubic field : phi/ { print $NF }' "$WORK/ofe.log")
		if [ -z "$e" ]; then echo "  $1 : the cubic-field identity did not run -> NG" >&2; status=1; return; fi
		echo "  $1 : linear and cubic field identities (cubic rel. error $e) -> OK"
	else
		if ! grep -q "cubic field : skipped" "$WORK/ofe.log"; then echo "  $1 : the cubic identity ran on a curved mesh -> NG" >&2; status=1; return; fi
		echo "  $1 : linear field identity on a curved mesh (cubic skipped) -> OK"
	fi
}
for m in box_hex3 box_prism3; do
	h3_c "$m" C
	compare "$m : C (parallel plate) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/h3_22.csv"
	for f in _41 _bin_22 _bin_41; do
		h3_c "$m$f" C
		grep -v '^title' "$WORK/rlc.csv" > "$WORK/h3_f.csv"
		if [ -s "$WORK/h3_22.csv" ] && cmp -s "$WORK/h3_22.csv" "$WORK/h3_f.csv"; then echo "  $m : 2.2 vs ${f#_} : rlc.csv identical -> OK"
		else echo "  $m : 2.2 vs ${f#_} : differ -> NG" >&2; status=1; fi
	done
	h3_c "$m" P
	h3_selftest "$m (analysis P)" cubic
done
# (b) 曲がった同軸 (C' = 2π ε0 εr / ln 3)
cex=$(awk 'BEGIN { pi = atan2(0, -1); printf "%.12e", 2 * pi * 8.8541878128e-12 * 2.1 / log(3) }')
for k in hex prism; do
	for r in a b; do
		sed -e "s/coax_hex2.msh/coax_${k}3$r.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/h3cx.ofe"
		(cd "$WORK" && "$OFE" -n 2 h3cx.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  coax_${k}3$r : the run failed -> NG" >&2; status=1; }
		eval "c_$r=\$(value_of C)"
	done
	tol=$( [ $k = hex ] && echo 1e-4 || echo 2e-4 )
	compare "curved coax ($k, order 3, finer) : C' [F/m]" "$c_b" "$cex" "$tol"
	res=$(awk -v a="$c_a" -v b="$c_b" -v e="$cex" 'BEGIN { ea = (a - e) / e; if (ea < 0) ea = -ea; eb = (b - e) / e; if (eb < 0) eb = -eb
		printf "%s (%.2e -> %.2e, ratio %.1f)", (ea > 12 * eb) ? "OK" : "NG", ea, eb, (eb > 0) ? ea / eb : 0 }')
	echo "  curved coax ($k, order 3) : C' converges at h^4 or faster (1/4 turn 2 -> 4 cells) : $res"
	case "$res" in NG*) status=1 ;; esac
done
sed -e "s/coax_hex2.msh/coax_hex3a.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = P"; next} {print}' > "$WORK/h3cx.ofe"
(cd "$WORK" && "$OFE" -n 2 h3cx.ofe > /dev/null 2>&1) || true
h3_selftest "curved coax (hexa64, analysis P)" skip
# (c) 熱連成
for m in slab_hex3:8 slab_prism3:16; do
	g=${m%%:*}; nf=${m#*:}
	printf 'OpenFEM 1 0\ntitle = order-3 slab cooled on one face\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 5000 20 0 2e-3 0 1e-3 0 0\ntemperature = 20\nthermal = 10\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$g" > "$WORK/h3th.ofe"
	th_run h3th.ofe
	compare "$g : bottom face cooled : T_max (parabola) [degC]" "$(th_t)" 125 1e-9
	compare "$g : bottom face cooled : T_mean [degC]" "$(t3_tm)" 123.333333333 1e-9
	th_nface "$g : convection on the bottom face only" "$nf"
	th_energy "$g"
done
set -- $(awk 'BEGIN { pi = atan2(0, -1); a = 0.5e-3; b = 1.5e-3; lz = 0.1e-3; s = 1e6; k = 10; h = 2000; I = 10
	A = pi * (b * b - a * a); q = (I / A) ^ 2 / s
	Pb = -q * b * b / (4 * k); Pa = -q * a * a / (4 * k)
	a11 = -k / b - h * log(b); a12 = -h; r1 = h * Pb - q * b / 2
	a21 = k / a - h * log(a); a22 = -h; r2 = h * Pa + q * a / 2
	d = a11 * a22 - a12 * a21; c1 = (r1 * a22 - a12 * r2) / d; c2 = (a11 * r2 - a21 * r1) / d
	n = 20000; num = 0; den = 0
	for (i = 0; i < n; i++) { r = a + (b - a) * (i + 0.5) / n; num += r * (-q * r * r / (4 * k) + c1 * log(r) + c2); den += r }
	printf "%.12e %.12e", lz / (s * A), num / den }')
h3ra=$1; h3ta=$2
for m in ann_hex1b ann_hex3a ann_hex3b ann_prism3a ann_prism3b; do
	printf 'OpenFEM 1 0\ntitle = order-3 annulus cooled on its curved faces\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 2000 20\ntemperature = 20\nthermal = 10\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$m" > "$WORK/h3an.ofe"
	th_run h3an.ofe
	th_energy "$m"
	eval "e_$m=\$(awk -v t=\"\$(t3_tm)\" -v x=\"\$h3ta\" 'BEGIN { d = (t - 20 - x) / x; if (d < 0) d = -d; printf \"%.6e\", d }')"
	case "$m" in *3b) compare "$m : R = L / (sigma pi (b^2 - a^2)) [ohm]" "$(th_r)" "$h3ra" 1e-4
	                  compare "$m : curved faces cooled : volume-mean temperature rise [K]" \
	                          "$(awk -v t="$(t3_tm)" 'BEGIN { printf "%.12e", t - 20 }')" "$h3ta" 1e-4 ;; esac
done
res=$(awk -v l="$e_ann_hex1b" -v ha="$e_ann_hex3a" -v hb="$e_ann_hex3b" -v pa="$e_ann_prism3a" -v pb="$e_ann_prism3b" 'BEGIN {
	ok = (ha > 12 * hb) && (pa > 12 * pb) && (l > 100 * hb)
	printf "%s (hexa64 %.2e -> %.2e, prism40 %.2e -> %.2e, hex8 %.2e)", ok ? "OK" : "NG", ha, hb, pa, pb, l }')
echo "  annulus : order 3 converges at h^4 and beats the linear facets : $res"
case "$res" in NG*) status=1 ;; esac
# (d) 弾くべき入力 : 電極面の次数違い、辺要素
awk '/^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 }
	e && ($2 == 36) && !done { s = $1 " 10"; for (i = 3; i <= 3 + $3 + 9; i++) s = s " " $i; print s; done = 1; next } { print }' \
	"$SRC/box_hex3.msh" > "$WORK/h3x_q9.msh"
sed -e "s/box_hex2.msh/h3x_q9.msh/" "$SRC/box_hex2.ofe" > "$WORK/h3x_q9.ofe"
mesh_reject "hexa64 with a 9-node electrode quadrilateral" h3x_q9.ofe
(cd "$WORK" && "$OFE" -n 2 h3x_q9.ofe 2>&1 | grep -q "not a 16-node quadrilateral") || { echo "  *** rejected for a different reason" >&2; status=1; }
awk '/^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 }
	e && ($2 == 21) && !done { s = $1 " 9"; for (i = 3; i <= 3 + $3 + 6; i++) s = s " " $i; print s; done = 1; next } { print }' \
	"$SRC/box_prism3.msh" > "$WORK/h3x_t6.msh"
sed -e "s/box_hex2.msh/h3x_t6.msh/" "$SRC/box_hex2.ofe" > "$WORK/h3x_t6.ofe"
mesh_reject "prism40 with a 6-node electrode triangle" h3x_t6.ofe
(cd "$WORK" && "$OFE" -n 2 h3x_t6.ofe 2>&1 | grep -q "not a 10-node triangle") || { echo "  *** rejected for a different reason" >&2; status=1; }
sed -e "s/box_hex2.msh/box_hex3.msh/" -e "s/^analysis.*/analysis = E/" "$SRC/box_hex2.ofe" > "$WORK/h3x_e.ofe"
mesh_reject "analysis E on a hexa64 mesh" h3x_e.ofe
(cd "$WORK" && "$OFE" -n 2 h3x_e.ofe 2>&1 | grep -q "need a first-order mesh") || { echo "  *** rejected for a different reason" >&2; status=1; }

# ---- 4 次の六面体 (125 節点、gmsh の hexa125) と角柱 (75 節点、prism75) ----
#
# 3 次と同じ形 (表 HEX4_I / PRISM4_A、各方向 5 点、角柱は面内 Duffy 25 点 x ζ 5 点)。4 次は曲がった
# 境界の幾何の誤差が h⁶ なので、**同じ h の 3 次に勝つことを assert してよい** (四面体の 4 次と同じ)
#   (a) 平行平板の C が全桁、4 形式で完全一致、4 次の場 (x⁴, x²y², yz³, xyz² を足す) の恒等式
#   (b) 曲がった同軸 : 1/4 周 2 -> 4 分割で 六面体 64 倍・角柱 75 倍 (h⁶、実測)、同じ格子の 3 次の
#       13.7 / 10.2 倍よい
#   (c) 熱連成 : 片面を冷やす板が機械精度、曲がった面を冷やす円環が同じ格子の 3 次の 10 倍よい
echo "[order-4 hex/prism] 125-node hexahedra / 75-node prisms (gmsh hexa125 / prism75)"
for m in box_hex4 box_prism4; do
	for f in "" _41 _bin_22 _bin_41; do cp "$SRC/$m$f.msh" "$WORK/"; done
done
for m in coax_hex4a coax_hex4b coax_prism4a coax_prism4b slab_hex4 slab_prism4 ann_hex4a ann_prism4a; do cp "$SRC/$m.msh" "$WORK/"; done
for m in box_hex4 box_prism4; do
	h3_c "$m" C
	compare "$m : C (parallel plate) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/h3_22.csv"
	for f in _41 _bin_22 _bin_41; do
		h3_c "$m$f" C
		grep -v '^title' "$WORK/rlc.csv" > "$WORK/h3_f.csv"
		if [ -s "$WORK/h3_22.csv" ] && cmp -s "$WORK/h3_22.csv" "$WORK/h3_f.csv"; then echo "  $m : 2.2 vs ${f#_} : rlc.csv identical -> OK"
		else echo "  $m : 2.2 vs ${f#_} : differ -> NG" >&2; status=1; fi
	done
	h3_c "$m" P
	if grep -q '^\*\*\*' "$WORK/ofe.log" || ! grep -q "quartic field : phi" "$WORK/ofe.log"; then
		echo "  $m (analysis P) : the quartic-field identity failed or did not run -> NG" >&2; status=1
	else
		echo "  $m (analysis P) : linear and quartic field identities ($(awk '/quartic field : phi/ { print $NF }' "$WORK/ofe.log")) -> OK"
	fi
done
for k in hex prism; do
	for r in 3a 4a 4b; do
		sed -e "s/coax_hex2.msh/coax_${k}$r.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/h3cx.ofe"
		(cd "$WORK" && "$OFE" -n 2 h3cx.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  coax_${k}$r : the run failed -> NG" >&2; status=1; }
		eval "c_$r=\$(value_of C)"
	done
	compare "curved coax ($k, order 4, finer) : C' [F/m]" "$c_4b" "$cex" 1e-5
	res=$(awk -v a="$c_4a" -v b="$c_4b" -v t="$c_3a" -v e="$cex" 'BEGIN {
		ea = (a - e) / e; if (ea < 0) ea = -ea; eb = (b - e) / e; if (eb < 0) eb = -eb; et = (t - e) / e; if (et < 0) et = -et
		ok = (ea > 40 * eb) && (et > 5 * ea)
		printf "%s (order 4 : %.2e -> %.2e, ratio %.1f; order 3 on the coarse mesh %.2e)", ok ? "OK" : "NG", ea, eb, (eb > 0) ? ea / eb : 0, et }')
	echo "  curved coax ($k) : order 4 converges at h^6 and beats order 3 on the same mesh : $res"
	case "$res" in NG*) status=1 ;; esac
done
for m in slab_hex4:8 slab_prism4:16; do
	g=${m%%:*}; nf=${m#*:}
	printf 'OpenFEM 1 0\ntitle = order-4 slab cooled on one face\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 5000 20 0 2e-3 0 1e-3 0 0\ntemperature = 20\nthermal = 10\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$g" > "$WORK/h3th.ofe"
	th_run h3th.ofe
	compare "$g : bottom face cooled : T_max (parabola) [degC]" "$(th_t)" 125 1e-9
	compare "$g : bottom face cooled : T_mean [degC]" "$(t3_tm)" 123.333333333 1e-9
	th_nface "$g : convection on the bottom face only" "$nf"
	th_energy "$g"
done
for m in ann_hex4a ann_prism4a; do
	printf 'OpenFEM 1 0\ntitle = order-4 annulus cooled on its curved faces\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 2000 20\ntemperature = 20\nthermal = 10\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$m" > "$WORK/h3an.ofe"
	th_run h3an.ofe
	th_energy "$m"
	eval "e_$m=\$(awk -v t=\"\$(t3_tm)\" -v x=\"\$h3ta\" 'BEGIN { d = (t - 20 - x) / x; if (d < 0) d = -d; printf \"%.6e\", d }')"
done
res=$(awk -v h3="$e_ann_hex3a" -v h4="$e_ann_hex4a" -v p3="$e_ann_prism3a" -v p4="$e_ann_prism4a" 'BEGIN {
	ok = (h4 <= 2e-4) && (p4 <= 2e-4) && (h3 > 5 * h4) && (p3 > 5 * p4)
	printf "%s (hexa125 %.2e vs hexa64 %.2e, prism75 %.2e vs prism40 %.2e)", ok ? "OK" : "NG", h4, h3, p4, p3 }')
echo "  annulus (coarse) : order 4 beats order 3 on the curved faces (geometry h^6 vs h^4) : $res"
case "$res" in NG*) status=1 ;; esac

# ---- 5 次の六面体 (216 節点、gmsh の hexa216 = 型 94) と角柱 (126 節点、prism126 = 型 106) ----
#
# 3 次 / 4 次と同じ形 (実測の表 HEX5_I / PRISM5_A、各方向 6 点)。表は gmsh の出力を分類する gmsh_hp_order.py で作り、
# 同じ分類が 3 次・4 次の表を全項再現することを先に確かめた。境界面は 36 節点四角形 (型 38) と 21 節点三角形。
#   (a) 平行平板の C が全桁、2.2 / バイナリ 4.1 で rlc.csv が一致、5 次の場の恒等式 (許容 1e-12、実測 4e-16 / 8e-15)。
#       **5 次の場は軸ごとの大きさで割る** (薄い箱で全体の大きさで割ると Z⁵ の項が 6e-6 しか効かず、積分点を 5 点に
#       落とす誤りが 1e-11 しか動かなかった、実測)
#   (b) 曲がった同軸 : 1/4 周 2 -> 4 分割で 六面体 166 倍・角柱 182 倍 (h⁶ 以上、実測)、同じ格子の 4 次の 14 / 13 倍よい
#   (c) 熱連成 : 片面を冷やす板が閉形式に機械精度 (四角形の面 / 三角形の面)
echo "[order-5 hex/prism] 216-node hexahedra / 126-node prisms (gmsh hexa216 / prism126)"
for m in box_hex5 box_hex5_bin_41 box_prism5 box_prism5_bin_41 coax_hex5a coax_hex5b coax_prism5a coax_prism5b slab_hex5 slab_prism5; do
	cp "$SRC/$m.msh" "$WORK/"
done
for m in box_hex5 box_prism5; do
	h3_c "$m" C
	compare "$m : C (parallel plate) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/h5_22.csv"
	rm -f "$WORK/rlc.csv"
	h3_c "${m}_bin_41" C
	grep -v '^title' "$WORK/rlc.csv" > "$WORK/h5_f.csv" 2> /dev/null
	if [ -s "$WORK/h5_22.csv" ] && cmp -s "$WORK/h5_22.csv" "$WORK/h5_f.csv"; then echo "  $m : 2.2 vs bin_41 : rlc.csv identical -> OK"
	else echo "  $m : 2.2 vs bin_41 : differ -> NG" >&2; status=1; fi
	h3_c "$m" P
	if grep -q '^\*\*\*' "$WORK/ofe.log" || ! grep -q "quintic field : phi" "$WORK/ofe.log"; then
		echo "  $m (analysis P) : the quintic-field identity failed or did not run -> NG" >&2; status=1
	else
		echo "  $m (analysis P) : linear and quintic field identities ($(awk '/quintic field : phi/ { print $NF }' "$WORK/ofe.log")) -> OK"
	fi
done
for k in hex prism; do
	for r in 4a 5a 5b; do
		sed -e "s/coax_hex2.msh/coax_${k}$r.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/h5cx.ofe"
		(cd "$WORK" && "$OFE" -n 2 h5cx.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  coax_${k}$r : the run failed -> NG" >&2; status=1; }
		eval "c_$r=\$(value_of C)"
	done
	compare "curved coax ($k, order 5, finer) : C' [F/m]" "$c_5b" "$cex" 1e-6
	res=$(awk -v a="$c_5a" -v b="$c_5b" -v t="$c_4a" -v e="$cex" 'BEGIN {
		ea = (a - e) / e; if (ea < 0) ea = -ea; eb = (b - e) / e; if (eb < 0) eb = -eb; et = (t - e) / e; if (et < 0) et = -et
		ok = (ea > 64 * eb) && (et > 5 * ea)
		printf "%s (order 5 : %.2e -> %.2e, ratio %.1f; order 4 on the coarse mesh %.2e)", ok ? "OK" : "NG", ea, eb, (eb > 0) ? ea / eb : 0, et }')
	echo "  curved coax ($k) : order 5 converges at h^6 or faster and beats order 4 on the same mesh : $res"
	case "$res" in NG*) status=1 ;; esac
done
for m in slab_hex5:8 slab_prism5:16; do
	g=${m%%:*}; nf=${m#*:}
	printf 'OpenFEM 1 0\ntitle = order-5 slab cooled on one face\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 5000 20 0 2e-3 0 1e-3 0 0\ntemperature = 20\nthermal = 10\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$g" > "$WORK/h5th.ofe"
	th_run h5th.ofe
	compare "$g : bottom face cooled : T_max (parabola) [degC]" "$(th_t)" 125 1e-9
	compare "$g : bottom face cooled : T_mean [degC]" "$(t3_tm)" 123.333333333 1e-9
	th_nface "$g : convection on the bottom face only" "$nf"
	th_energy "$g"
done

# ---- 高次の種別混在 (同じ次数の四面体・六面体・角柱) ----
#
# 面のトレースは三角形なら P_p、四角形なら Q_p なので、同じ次数の種別どうしは面で適合する
# (unstruct.c の elem_finish)。格子は gmsh 製 (gmsh_coaxmix_ho.py : 同軸の内側のリングを角柱、外側を
# 六面体 / gmsh_boxmix_ho.py : 下半分を角柱、上半分を四面体)。
#   (a) 角柱 + 六面体の曲がった同軸 : 2 次 h⁴ (比 15.5)、3 次 h⁴ 以上 (比 42、実測)
#   (b) 角柱 + 四面体の箱 : C が全桁、analysis = P の線形場の恒等式 (種別ごとの要素の行列の振り分け)
#   (c) 熱連成 : 三角形の露出面が四面体のものか角柱のものかで振り分ける (混在で初めて両方が出る)。
#       恒等式と露出面の数
#   (d) 弾く : 次数の違う種別、serendipity と完全 Lagrange の混在、四角形面を三角形で覆う界面
echo "[mixed high order] tetrahedra / hexahedra / prisms of the same order in one mesh"
for m in coaxmix2a coaxmix2b coaxmix3a coaxmix3b boxmix2 boxmix3 coaxmixt2a; do cp "$SRC/$m.msh" "$WORK/"; done
for o in 2 3; do
	for r in a b; do
		sed -e "s/coax_hex2.msh/coaxmix$o$r.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/mxcx.ofe"
		(cd "$WORK" && "$OFE" -n 2 mxcx.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  coaxmix$o$r : the run failed -> NG" >&2; status=1; }
		eval "c_$r=\$(value_of C)"
	done
	tol=$( [ $o = 2 ] && echo 5e-4 || echo 2e-5 )
	compare "prisms + hexahedra, curved coax (order $o, finer) : C' [F/m]" "$c_b" "$cex" "$tol"
	res=$(awk -v a="$c_a" -v b="$c_b" -v e="$cex" 'BEGIN { ea = (a - e) / e; if (ea < 0) ea = -ea; eb = (b - e) / e; if (eb < 0) eb = -eb
		printf "%s (%.2e -> %.2e, ratio %.1f)", (ea > 12 * eb) ? "OK" : "NG", ea, eb, (eb > 0) ? ea / eb : 0 }')
	echo "  prisms + hexahedra, curved coax (order $o) : C' converges at h^4 or faster : $res"
	case "$res" in NG*) status=1 ;; esac
done
for o in 2 3; do
	h3_c "boxmix$o" C
	compare "prisms + tetrahedra (order $o) : C (parallel plate) [F]" "$(value_of C)" 1.77083756e-13 1e-8
	h3_c "boxmix$o" P
	if grep -q '^\*\*\*' "$WORK/ofe.log" || ! grep -q "mixed, order $o) self test" "$WORK/ofe.log" || ! grep -q "linear field : .*rel. error" "$WORK/ofe.log"; then
		echo "  prisms + tetrahedra (order $o, analysis P) : the self test failed -> NG" >&2; status=1
	else
		echo "  prisms + tetrahedra (order $o, analysis P) : linear field identity ($(awk '/linear field : / { print $NF }' "$WORK/ofe.log")) -> OK"
	fi
done
for m in boxmix2:36 boxmix3:36 coaxmix2a:48 coaxmix3a:48; do
	g=${m%%:*}; nf=${m#*:}
	printf 'OpenFEM 1 0\ntitle = mixed high-order thermal\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 4e-3 20\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nheatsink = 0 20\nconvection = 3000 20\nradiation = 0.5 20\ntemperature = 20\nthermal = 5\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$g" > "$WORK/mxth.ofe"
	th_run mxth.ofe
	th_energy "$g : convection + radiation"
	th_nface "$g : exposed faces (triangles of tetrahedra / prisms, quadrilaterals)" "$nf"
done
# (d)
awk '/^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 }
	e && ($2 == 13) { s = $1 " 6"; for (i = 3; i <= 3 + $3 + 6; i++) s = s " " $i; print s; next } { print }' "$SRC/coaxmix2a.msh" > "$WORK/mxx_ord.msh"
#     次数の違う種別は最高次数に上げる ([mixed orders])。角柱を 1 次に落とした同軸 : 角柱は 1 次の形のまま (内側の円が
#     8 角形、面積の等しい半径 0.949 a で C' −4.6% の見積もり) なので許容 5% (実測 −3.6%)
sed -e "s/coax_hex2.msh/mxx_ord.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/mxx_ord.ofe"
rm -f "$WORK/rlc.csv"
(cd "$WORK" && "$OFE" -n 2 mxx_ord.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  mxx_ord : the run failed -> NG" >&2; status=1; }
grep -q "^Mixed orders : 16 elements of order 1 and up raised to order 2" "$WORK/ofe.log" || { echo "  hex27 + prism6 : the prisms were not raised -> NG" >&2; status=1; }
compare "hex27 + prism6 raised to order 2 (prisms keep their first-order shape) : C' [F/m]" "$(value_of C)" "$cex" 5e-2
awk '/^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 }
	e && ($2 == 12) { s = $1 " 17"; for (i = 3; i <= 3 + $3 + 20; i++) s = s " " $i; print s; next } { print }' "$SRC/coaxmix2a.msh" > "$WORK/mxx_fam.msh"
#     serendipity と完全の混在も 27 節点に上げる : 27 節点の六面体から面心・体心を落とした格子を上げると元の格子に戻る
sed -e "s/coax_hex2.msh/mxx_fam.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/mxx_fam.ofe"
sed -e "s/coax_hex2.msh/coaxmix2a.msh/" "$SRC/coax_hex2.ofe" | awk '/^analysis/{print "analysis = C"; next} {print}' > "$WORK/mxx_fam0.ofe"
mx_same "hex20 + prism18 (serendipity and complete) raised to hex27" mxx_fam0.ofe mxx_fam.ofe
sed -e "s/coax_hex2.msh/coaxmixt2a.msh/" "$SRC/coax_hex2.ofe" > "$WORK/mxx_nc.ofe"
mesh_reject "prism quadrilaterals covered by tetrahedron triangles (order 2)" mxx_nc.ofe
(cd "$WORK" && "$OFE" -n 2 mxx_nc.ofe 2>&1 | grep -q "covered by triangular faces") || { echo "  *** rejected for a different reason" >&2; status=1; }

# ---- 交流の発熱 (thermal + analysis = A) ----
#
# 3 次元渦電流 (A-φ) の時間平均の損失 Q_e = ½ eᴴ T_e e (e = jωa + Gφ) で温度を解く。
# Σ Q_e = ½ Re(V I*) は連成系の行から厳密に出る (eddy3d.c の注釈) ので、
# エネルギーの恒等式がそのまま損失の積分の検査になる。電流は波高値。
echo "[thermal_ac] AC heating : eddy-current losses (A-phi) in the electro-thermal coupling"
cp "$SRC/thermal_ac.ofe" "$SRC/bar_tet.msh" "$WORK/"
# (1) 100 Hz (表皮深さ 50mm ≫ 断面) : ω -> 0 で A-φ は DC の導電解析に帰着するので、
#     DC の熱連成 (thermal_tet.ofe、I = 7.5 A) と波高値 7.5√2 A の交流が一致する。
#     α = 4e-3 なので σ(T) を要素ごとに A-φ へ渡していなければ R が 9% ずれる。
#     許容 1e-8 : 低周波の A-φ と DC の差 (実測 9e-10)
th_run thermal_tet.ofe
rdc=$(th_r); tdc=$(th_t)
awk '/^analysis/{print "frequency = 100"; print "awall = 20"; print "analysis = A"; next}
	/^thermal/{print "thermal = 10.606601717798213"; next} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/tac_lf.ofe"
th_run tac_lf.ofe
compare "AC at 100 Hz (I = 7.5 sqrt 2 A peak) == DC thermal (7.5 A) : R [ohm]" "$(th_r)" "$rdc" 1e-8
compare "AC at 100 Hz == DC thermal : T_max [degC]" "$(th_t)" "$tdc" 1e-8
th_energy "AC at 100 Hz (joule = I^2 R / 2 with the peak current)"
# (2) 100 kHz の銅の棒 (thermal_ac.ofe、表皮深さ 0.21mm)、α = 0。
#     R・X は熱連成なしの A (bar_eddy.ofe と同じ格子・材料) と一致し、中央の断面の
#     表面と中心の温度差が表皮効果の発熱分布の Fourier 級数解 (ac_heat.awk) に合う。
#     損失を要素で一様に配ったり (分布を潰す)、e の Gφ の符号を誤ると差が崩れる。
#     許容 0.3% : 24 分割で +0.17% (12 / 24 / 48 分割で +0.35 / +0.17 / +0.05%、実測)
awk '/^tempco/{print "tempco = 2 0 0"; next} /^thermal/{print "thermal = 100"; print "fieldout = 1"; next} {print}' \
	"$SRC/thermal_ac.ofe" > "$WORK/tac_a0.ofe"
th_run tac_a0.ofe
ra0=$(th_r)
xa0=$(awk '/^  X = / { print $3; exit }' "$WORK/ofe.log")
tma0=$(awk '/T_mean/ { print $8; exit }' "$WORK/ofe.log")
th_energy "AC, 100 kHz, alpha = 0"
dfe=$(awk -v xc=1e-3 -v z1=0 -v z2=0.5e-3 '
	/^POINTS/ { np = $2; pt = 1; n = 0; next }
	pt { for (i = 1; i <= NF; i++) { c[n] = $i; n++ } if (n >= 3 * np) pt = 0; next }
	/^SCALARS T / { st = 1; m = 0; next }
	st && /^LOOKUP/ { next }
	st && (NF != 1) { st = 0; next }
	st { x = c[3 * m]; z = c[3 * m + 2]; if (((x - xc) ^ 2 < 1e-18) && ((z - z1) ^ 2 < 1e-18)) { s1 += $1; n1++ } if (((x - xc) ^ 2 < 1e-18) && ((z - z2) ^ 2 < 1e-18)) { s2 += $1; n2++ } m++; if (m >= np) st = 0 }
	END { if ((n1 > 0) && (n2 > 0)) printf "%.12e", s1 / n1 - s2 / n2 }' "$WORK/ofe_field.vtk")
set -- $(awk -v f=1e5 -v sig=5.8e7 -v k=400 -v I=100 -v R="$ra0" -f "$SRC/ac_heat.awk")
compare "AC, 100 kHz : T(L/2, surface) - T(L/2, centre) vs Fourier series [K]" "$dfe" "$1" 0.003
compare "AC, 100 kHz : volume-mean temperature rise == q L^2 / (12 k) [K]" \
	"$(awk -v t="$tma0" 'BEGIN { printf "%.12e", t - 20 }')" "$2" 0.003
awk '!/^tempco|^kthermal|^heatsink|^temperature|^thermal/' "$SRC/thermal_ac.ofe" > "$WORK/tac_unc.ofe"
th_run tac_unc.ofe
zline=$(awk '/Z\(1,1\)/ { print; exit }' "$WORK/ofe.log")
compare "AC, alpha = 0 : R == the uncoupled A [ohm]" "$ra0" "$(echo "$zline" | awk '{ print $3 }')" 1e-6
compare "AC, alpha = 0 : X == the uncoupled A [ohm]" "$xa0" "$(echo "$zline" | awk '{ print $4 }' | tr -d '+')" 1e-6
# (3) α > 0 (thermal_ac.ofe) : 恒等式と収束。断熱の棒の過渡 (CN) で発熱 = 蓄熱
th_run thermal_ac.ofe
th_energy "thermal_ac (alpha = 3.9e-3, 300 A peak)"
compare "thermal_ac : R rises above the alpha = 0 value by the heating [ohm]" \
	"$(awk -v a="$(th_r)" -v b="$ra0" 'BEGIN { printf "%.12e", ((a > b) ? 1 : 0) }')" 1 0
awk '/^heatsink/{next} /^thermal/{print; print "heatcap = 2 3.45e6"; print "transient = 0.01 20 0.5"; next} {print}' \
	"$SRC/thermal_ac.ofe" > "$WORK/tac_tr.ofe"
th_run tac_tr.ofe
th_balance "thermal_ac, adiabatic transient (CN)"
# (2b) 10 節点四面体 (bar_tet2.msh) の交流の発熱 : 2 次の Nedelec (eddy3d_p2.c)。損失 ½σ|e|² を組み立てと同じ
#      積分点で節点の形状関数に掛けるので Σf = ½I²R が厳密。
#      1 Hz (表皮効果の補正 3e-10) の交流が同じ格子の DC の熱連成に 1e-8 で一致し、100 kHz・α = 0 では R が
#      熱連成なしの 2 次の A に一致、中央の断面の温度差が Fourier 級数解に 0.1% 以内 (実測 +0.031%、1 次 +0.17%)
cp "$SRC/bar_tet2.msh" "$WORK/"
sed 's/^mesh = .*/mesh = bar_tet2.msh/' "$SRC/thermal_ac.ofe" > "$WORK/tac2.ofe"
awk '/^frequency/{next} /^awall/{next} /^analysis/{print "analysis = R"; next} /^thermal/{print "thermal = 7.5"; next} {print}' "$WORK/tac2.ofe" > "$WORK/tac2_dc.ofe"
th_run tac2_dc.ofe
rdc2=$(th_r); tdc2=$(th_t)
awk '/^frequency/{print "frequency = 1"; next} /^thermal/{print "thermal = 10.606601717798213"; next} {print}' "$WORK/tac2.ofe" > "$WORK/tac2_lf.ofe"
th_run tac2_lf.ofe
compare "10-node tets : AC at 1 Hz (7.5 sqrt 2 A peak) == DC thermal (7.5 A) : R [ohm]" "$(th_r)" "$rdc2" 1e-8
compare "10-node tets : AC at 1 Hz == DC thermal : T_max [degC]" "$(th_t)" "$tdc2" 1e-8
th_energy "10-node tets, AC at 1 Hz"
awk '/^tempco/{print "tempco = 2 0 0"; next} /^thermal/{print "thermal = 100"; print "fieldout = 1"; next} {print}' "$WORK/tac2.ofe" > "$WORK/tac2_a0.ofe"
th_run tac2_a0.ofe
ra02=$(th_r)
tma02=$(awk '/T_mean/ { print $8; exit }' "$WORK/ofe.log")
th_energy "10-node tets, AC, 100 kHz, alpha = 0"
dfe2=$(awk -v xc=1e-3 -v z1=0 -v z2=0.5e-3 '
	/^POINTS/ { np = $2; pt = 1; n = 0; next }
	pt { for (i = 1; i <= NF; i++) { c[n] = $i; n++ } if (n >= 3 * np) pt = 0; next }
	/^SCALARS T / { st = 1; m = 0; next }
	st && /^LOOKUP/ { next }
	st && (NF != 1) { st = 0; next }
	st { x = c[3 * m]; z = c[3 * m + 2]; if (((x - xc) ^ 2 < 1e-18) && ((z - z1) ^ 2 < 1e-18)) { s1 += $1; n1++ } if (((x - xc) ^ 2 < 1e-18) && ((z - z2) ^ 2 < 1e-18)) { s2 += $1; n2++ } m++; if (m >= np) st = 0 }
	END { if ((n1 > 0) && (n2 > 0)) printf "%.12e", s1 / n1 - s2 / n2 }' "$WORK/ofe_field.vtk")
set -- $(awk -v f=1e5 -v sig=5.8e7 -v k=400 -v I=100 -v R="$ra02" -f "$SRC/ac_heat.awk")
compare "10-node tets, AC, 100 kHz : T(L/2, surface) - T(L/2, centre) vs Fourier series [K]" "$dfe2" "$1" 1e-3
compare "10-node tets, AC, 100 kHz : volume-mean temperature rise == q L^2 / (12 k) [K]" \
	"$(awk -v t="$tma02" 'BEGIN { printf "%.12e", t - 20 }')" "$2" 1e-6
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_tet2.msh"; next} {print}' > "$WORK/tac2_unc.ofe"
th_run tac2_unc.ofe
compare "10-node tets, AC, alpha = 0 : R == the uncoupled second-order A [ohm]" "$ra02" \
	"$(awk '/Z\(1,1\)/ { print $3; exit }' "$WORK/ofe.log")" 1e-6
awk '/^heatsink/{next} /^thermal/{print; print "heatcap = 2 3.45e6"; print "transient = 0.01 20 0.5"; next} {print}' \
	"$WORK/tac2.ofe" > "$WORK/tac2_tr.ofe"
th_run tac2_tr.ofe
th_balance "10-node tets, AC (alpha = 3.9e-3), adiabatic transient (CN)"
# (2c) 3 次の四面体 (bar_ho3.msh) の交流の発熱 : 次数 r の Nedelec (eddy3d_pr.c)。2 次と同じ形 (損失を組み立てと同じ積分点で
#      節点の形状関数に掛けるので Σf = ½I²R が厳密)。1 Hz の交流 = DC の熱連成 (1e-8)、α = 0 で R = 熱連成なしの 3 次の A、
#      体積平均の温度上昇 = q L²/(12 k) (格子は gmsh の非構造なので断面の点の温度差は見ない)
cp "$SRC/bar_ho3.msh" "$WORK/"
sed 's/^mesh = .*/mesh = bar_ho3.msh/' "$SRC/thermal_ac.ofe" > "$WORK/tac3.ofe"
awk '/^frequency/{next} /^awall/{next} /^analysis/{print "analysis = R"; next} /^thermal/{print "thermal = 7.5"; next} {print}' "$WORK/tac3.ofe" > "$WORK/tac3_dc.ofe"
th_run tac3_dc.ofe
rdc3=$(th_r); tdc3=$(th_t)
awk '/^frequency/{print "frequency = 1"; next} /^thermal/{print "thermal = 10.606601717798213"; next} {print}' "$WORK/tac3.ofe" > "$WORK/tac3_lf.ofe"
th_run tac3_lf.ofe
compare "order-3 tets : AC at 1 Hz (7.5 sqrt 2 A peak) == DC thermal (7.5 A) : R [ohm]" "$(th_r)" "$rdc3" 1e-8
compare "order-3 tets : AC at 1 Hz == DC thermal : T_max [degC]" "$(th_t)" "$tdc3" 1e-8
th_energy "order-3 tets, AC at 1 Hz"
awk '/^tempco/{print "tempco = 2 0 0"; next} /^thermal/{print "thermal = 100"; next} {print}' "$WORK/tac3.ofe" > "$WORK/tac3_a0.ofe"
th_run tac3_a0.ofe
ra03=$(th_r)
tma03=$(awk '/T_mean/ { print $8; exit }' "$WORK/ofe.log")
th_energy "order-3 tets, AC, 100 kHz, alpha = 0"
set -- $(awk -v f=1e5 -v sig=5.8e7 -v k=400 -v I=100 -v R="$ra03" -f "$SRC/ac_heat.awk")
compare "order-3 tets, AC, 100 kHz : volume-mean temperature rise == q L^2 / (12 k) [K]" \
	"$(awk -v t="$tma03" 'BEGIN { printf "%.12e", t - 20 }')" "$2" 1e-6
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_ho3.msh"; next} {print}' > "$WORK/tac3_unc.ofe"
th_run tac3_unc.ofe
compare "order-3 tets, AC, alpha = 0 : R == the uncoupled order-3 A [ohm]" "$ra03" \
	"$(awk '/Z\(1,1\)/ { print $3; exit }' "$WORK/ofe.log")" 1e-6
# (4) 弾くべき入力
awk '/^frequency/{print "frequencysweep = 1e4 1e5"; next} {print}' "$SRC/thermal_ac.ofe" > "$WORK/tacx_sw.ofe"
th_reject "AC heating with frequencysweep" tacx_sw.ofe "cannot be combined with frequencysweep"
awk '/^kthermal/{print; print "bh = 2 100 0.5"; print "bh = 2 1000 1.0"; next} {print}' "$SRC/thermal_ac.ofe" > "$WORK/tacx_bh.ofe"
th_reject "AC heating with a nonlinear (bh) material" tacx_bh.ofe "linear materials only"
awk '/^analysis/{print "analysis = A"; print "frequency = 1e3"; next} {print}' "$SRC/thermal_bar.ofe" > "$WORK/tacx_grid.ofe"
th_reject "AC heating on the structured grid" tacx_grid.ofe "requires an unstructured mesh"
# (5) 六面体・角柱 (辺要素の種別に依らない損失 ½ eᴴ T_e e)。awall の四角形面が要る
#     (無いと非受動 R になる。edge.c の注釈)。σ = 1e6 の棒 (表皮深さ 50mm) なら 100 Hz の
#     A は DC の熱連成と 7e-10 で一致 (許容 1e-8)。混在 (六面体 + 角柱) は界面で損失を要素で
#     一様に配る近似が DC の積分点の重みと 2e-8 違うので恒等式だけ
cp "$SRC/bar_hex.msh" "$SRC/bar_prism.msh" "$SRC/bar_hexprism.msh" "$SRC/bar_hex12.msh" "$SRC/ac_heat.awk" "$WORK/"
for m in bar_hex bar_prism bar_hexprism; do
	awk -v m="$m" '/^mesh/{print "mesh = " m ".msh"; next} /^region/{print; print "region = 2 2"; next}
		/^thermal/{print "thermal = 3"; next} {print}' "$SRC/thermal_tet.ofe" > "$WORK/tach_dc.ofe"
	th_run tach_dc.ofe
	rdc=$(th_r); tdc=$(th_t)
	awk '/^analysis/{print "frequency = 100"; print "awall = 20"; print "analysis = A"; next}
		/^thermal/{print "thermal = 4.242640687119285"; next} {print}' "$WORK/tach_dc.ofe" > "$WORK/tach_lf.ofe"
	th_run tach_lf.ofe
	if [ "$m" != bar_hexprism ]; then
		compare "$m : AC at 100 Hz (I = 3 sqrt 2 A peak) == DC thermal (3 A) : R [ohm]" "$(th_r)" "$rdc" 1e-8
		compare "$m : AC at 100 Hz == DC thermal : T_max [degC]" "$(th_t)" "$tdc" 1e-8
	fi
	th_energy "$m : AC at 100 Hz"
done
#     100 kHz の銅の棒 (bar_hex.msh、thermal_ac.ofe と同じ形)、α = 0 : 表面と中心の温度差の
#     Fourier 級数解との誤差が h² で縮む (12 -> 24 分割で比 4.04、24 で -1.1%、48 で -0.26%、
#     実測。六面体は損失を 8 節点に一様に配るので四面体より定数が大きい)。体積平均の
#     温度上昇は 0.3% で一致
tach_dt() {	# tach_dt <格子> -> "表面と中心の温度差 Fourier の値 平均の上昇 閉形式"
	awk -v m="$1" '/^mesh/{print "mesh = " m ".msh"; next} /^region/{print; print "region = 2 2"; next}
		/^tempco/{print "tempco = 2 0 0"; next} /^thermal/{print "thermal = 100"; print "fieldout = 1"; next} {print}' \
		"$SRC/thermal_ac.ofe" > "$WORK/tach_a0.ofe"
	th_run tach_a0.ofe
	fd=$(awk -v xc=1e-3 -v z1=0 -v z2=0.5e-3 '
		/^POINTS/ { np = $2; pt = 1; n = 0; next }
		pt { for (i = 1; i <= NF; i++) { c[n] = $i; n++ } if (n >= 3 * np) pt = 0; next }
		/^SCALARS T / { st = 1; m = 0; next }
		st && /^LOOKUP/ { next }
		st && (NF != 1) { st = 0; next }
		st { x = c[3 * m]; z = c[3 * m + 2]; m++; if (m >= np) st = 0 }
		st && ((x - xc) ^ 2 < 1e-18) && ((z - z1) ^ 2 < 1e-18) { s1 += $1; n1++ }
		st && ((x - xc) ^ 2 < 1e-18) && ((z - z2) ^ 2 < 1e-18) { s2 += $1; n2++ }
		END { if ((n1 > 0) && (n2 > 0)) printf "%.12e", s1 / n1 - s2 / n2 }' "$WORK/ofe_field.vtk")
	tm=$(awk '/T_mean/ { print $8; exit }' "$WORK/ofe.log")
	set -- $(awk -v f=1e5 -v sig=5.8e7 -v k=400 -v I=100 -v R="$(th_r)" -f "$SRC/ac_heat.awk")
	echo "$fd $1 $(awk -v t="$tm" 'BEGIN { printf "%.12e", t - 20 }') $2"
}
set -- $(tach_dt bar_hex12); d12=$1; f12=$2
set -- $(tach_dt bar_hex); d24=$1; f24=$2
compare "hexahedra, 100 kHz : volume-mean temperature rise == q L^2 / (12 k) [K]" "$3" "$4" 0.003
res=$(awk -v a="$d12" -v fa="$f12" -v b="$d24" -v fb="$f24" 'BEGIN {
	e1 = (a - fa) / fa; if (e1 < 0) e1 = -e1
	e2 = (b - fb) / fb; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.2f)", ((e2 < 0.02) && (q > 3.6) && (q < 4.6)) ? "OK" : "NG", e1, e2, q }')
echo "  hexahedra, 100 kHz : surface - centre temperature vs Fourier series, h^2 (12 -> 24 cells) : $res"
case "$res" in NG*) status=1 ;; esac
#     ピラミッドは A が弾く (多項式の Nedelec 基底が無い)
awk '/^mesh/{print "mesh = box_pyr.msh"; next} /^region/{print "region = 1 2"; next} {print}' "$SRC/thermal_ac.ofe" > "$WORK/tacx_pyr.ofe"
cp "$SRC/box_pyr.msh" "$WORK/"
mesh_reject "AC heating on a pyramid mesh" tacx_pyr.ofe

# ---- 交流の発熱 (thermal + analysis = F、断面 2 次元・単位長あたり) ----
#
# 1 次三角形の熱 (剛性・質量・露出辺の対流 / 輻射、ヒートシンクは導体の三角形の節点) と、
# 要素ごとの σ(T) で単一ポートの F を電流で駆動して解く eddy2d_heat_solve (solve.c の注釈)。
# 発熱は組み立てと同じ質量行列の上の ½∫σ|V' - jωAz|² なので Σf = ½|I|²R(f) が厳密。
echo "[thermal F 2-D] AC heating on a 2-D cross-section (per unit length)"
cp "$SRC/thermal_f2d.ofe" "$SRC/plate2d.msh" "$SRC/plate2d_p2.msh" "$SRC/plate_line_ac.ofe" "$WORK/"
f2_tm() { awk '/^  T_max = / { print $(NF-1); exit }' "$WORK/ofe.log"; }
# (1) 100 Hz、h = 50 : 等温の板の閉形式 (thermal_f2d.ofe の注釈)。許容は板の中の温度差
#     (Bi = ht/k = 6e-6 と端の効果、実測 2e-6) から ΔT に 1e-5、R に 1e-6
set -- $(awk 'BEGIN { I = 4.5; s0 = 5.8e7; a = 3.9e-3; W = 1e-3; t = 0.05e-3; h = 50; L = 2 * (W + t)
	c = 0.5 * I * I / s0 / (W * t * h * L); dT = c / (1 - a * c)
	printf "%.12e %.12e", dT, 2 * (1 + a * dT) / (s0 * W * t) }')
dtex=$1; rex=$2
th_run thermal_f2d.ofe
compare "2-D AC heating, isothermal plates : T_mean - T_a [K]" "$(awk -v t="$(f2_tm)" 'BEGIN { printf "%.12e", t - 20 }')" "$dtex" 1e-5
compare "2-D AC heating, isothermal plates : R(f) = 2 rho(T) / (W t) [ohm/m]" "$(th_r)" "$rex" 1e-6
th_nface "2-D AC heating : exposed edges 2 x (2 x 40 + 2 x 12)" 208
# (2) エネルギーの恒等式は h = 500 (I を √10 倍で同じ ΔT) で見る。h = 50 では k/h の比が
#     大きすぎて熱の CG の定数分の残差が 6e-10 残る (実測。許容 1e-9 の上限に近い)
awk '/^convection/{print "convection = 500 20"; next} /^thermal/{print "thermal = 14.230249470757707 200 1e-11"; next} {print}' \
	"$SRC/thermal_f2d.ofe" > "$WORK/f2_h500.ofe"
th_run f2_h500.ofe
th_energy "2-D AC heating, 100 Hz"
#     10 MHz (表皮効果・近接効果で発熱が間隙側に偏る)
awk '/^frequency/{print "frequency = 1e7"; next} {print}' "$WORK/f2_h500.ofe" > "$WORK/f2_10m.ofe"
th_run f2_10m.ofe
th_energy "2-D AC heating, 10 MHz"
r10=$(th_r)
# (3) α = 0 : 熱連成の R(f) は熱連成なしの F と一致 (σ(T) の配線の検査。rlc.csv の 8 桁)
for fq in 100 1e7; do
	awk -v f="$fq" '/^thermal|^convection/{next} /^frequency/{print "frequency = " f; next} {print}' \
		"$SRC/thermal_f2d.ofe" > "$WORK/f2_plain.ofe"
	(cd "$WORK" && "$OFE" -n 2 f2_plain.ofe > /dev/null && "$OFE_POST" > /dev/null)
	rp=$(value_of Rf)
	awk -v f="$fq" '/^frequency/{print "frequency = " f; next} /^conductortempco/{print "conductortempco = " $3 " 0 20"; next} {print}' \
		"$SRC/thermal_f2d.ofe" > "$WORK/f2_a0.ofe"
	th_run f2_a0.ofe
	compare "2-D AC heating, alpha = 0 == plain F : R(f = $fq) [ohm/m]" "$(th_r)" "$rp" 1e-8
	[ "$fq" = 1e7 ] && r10a0=$(th_r)
done
res=$(awk -v a="$r10" -v b="$r10a0" 'BEGIN { printf "%s (R %.6e with heating vs %.6e at alpha = 0)", (a > b) ? "OK" : "NG", a, b }')
echo "  2-D AC heating, 10 MHz : R(f) rises with the temperature : $res"
case "$res" in NG*) status=1 ;; esac
# (4) 輻射の辺 (2 点 Gauss の Newton) と、α = 0 での反復回数 (接線の検査)
awk '/^convection/{print; print "radiation = 0.8 20"; next} {print}' "$WORK/f2_h500.ofe" > "$WORK/f2_rad.ofe"
th_run f2_rad.ofe
th_energy "2-D AC heating, convection + radiation"
awk '/^conductortempco/{print "conductortempco = " $3 " 0 20"; next} {print}' "$WORK/f2_rad.ofe" > "$WORK/f2_rad0.ofe"
th_run f2_rad0.ofe
nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
if [ -n "$nit" ] && [ "$nit" -le 5 ]; then
	echo "  2-D : Newton on edge radiation (alpha = 0) : converged in $nit iterations (<= 5) -> OK"
else
	echo "  2-D : Newton on edge radiation (alpha = 0) : ${nit:-?} iterations (> 5) -> NG" >&2
	status=1
fi
# (5) ヒートシンク (導体 0 の三角形の節点を 20 degC に固定。導体 1 は間隙が k = 0 なので
#     対流でしか冷えない) と過渡 (CN) の収支
awk '/^convection/{print; print "heatsink = 0 20"; next} {print}' "$WORK/f2_h500.ofe" > "$WORK/f2_sink.ofe"
th_run f2_sink.ofe
th_energy "2-D AC heating, heat sink on conductor 0"
awk '/^thermal/{print; print "heatcap = 3 3.4e6"; print "transient = 0.05 20 0.5"; next} {print}' "$WORK/f2_h500.ofe" > "$WORK/f2_tr.ofe"
th_run f2_tr.ofe
th_balance "2-D AC heating, transient (CN)"
# (7) 熱伝導の剛性 : 間隙の誘電体に k_d = 0.5 を与え、導体 0 をヒートシンク、対流なし。
#     板 1 の熱はすべて間隙を 1 次元に通る (側面は断熱) ので
#       ΔT = c / (1 - α c),  c = ½ I² ρ0 d / (W t k_d W)
#     に板 1 の中の放物線 q t² / (2 k_cu) (q = 板 1 の発熱密度) を足したものが T_max。
#     面内で 90 度回した格子 (plate2d_rot.msh) でも同じ答え (剛性の面内 2 成分の両方を見る)。
#     許容は ΔT に 1e-5 (σ(T) を板の平均温度で取る近似の 1e-7 と、板の中の 1 次要素の誤差)
cp "$SRC/plate2d_rot.msh" "$WORK/"
set -- $(awk 'BEGIN { I = 4.5; s0 = 5.8e7; a = 3.9e-3; W = 1e-3; t = 0.05e-3; d = 0.2e-3; kd = 0.5; kc = 400
	c = 0.5 * I * I / s0 * d / (W * t * kd * W); dT = c / (1 - a * c)
	q = dT * kd * W / d / (W * t)
	printf "%.12e", dT + q * t * t / (2 * kc) }')
for m in plate2d plate2d_rot; do
	awk -v m="$m" '/^mesh/{print "mesh = " m ".msh"; next} /^kthermal/{print; print "kthermal = 2 0.5"; next}
		/^convection/{print "heatsink = 0 20"; next} {print}' "$SRC/thermal_f2d.ofe" > "$WORK/f2_kd.ofe"
	th_run f2_kd.ofe
	compare "2-D AC heating, conduction across the gap ($m) : T_max - T_sink [K]" \
		"$(awk -v t="$(th_t)" 'BEGIN { printf "%.12e", t - 20 }')" "$1" 1e-5
	th_energy "2-D AC heating, conduction across the gap ($m)"
done
#     断熱 (冷却なし)・α = 0 の過渡 (CN) : 面積加重の平均温度は質量の行和 1ᵀM T/ρc そのもの
#     なので、T_mean(t) = T_0 + ½ I² R t / (ρc 2 W t_p) が θ に依らず厳密 (質量の (1 + δ) の
#     重みの検査。収支の恒等式は質量の誤りを見ない)
awk '/^convection/{next} /^conductortempco/{print "conductortempco = " $3 " 0 20"; next}
	/^thermal/{print; print "heatcap = 3 3.4e6"; print "transient = 0.2 10 0.5"; next} {print}' \
	"$SRC/thermal_f2d.ofe" > "$WORK/f2_ad.ofe"
th_run f2_ad.ofe
compare "2-D AC heating, adiabatic transient (CN) : T_mean(t) - T_0 == I^2 R t / (2 rho c A) [K]" \
	"$(awk -v t="$(f2_tm)" 'BEGIN { printf "%.12e", t - 20 }')" \
	"$(awk -v r="$(th_r)" 'BEGIN { printf "%.12e", 0.5 * 4.5 * 4.5 * r * 0.2 / (3.4e6 * 2 * 1e-3 * 0.05e-3) }')" 1e-9
th_balance "2-D AC heating, adiabatic transient (CN)"
# (6) 弾くべき入力
awk '/^mesh/{print "mesh = plate2d_p2.msh"; next} {print}' "$SRC/thermal_f2d.ofe" > "$WORK/f2x_p2.ofe"
mesh_reject "AC heating (F) on 6-node triangles" f2x_p2.ofe
if ! grep -q "needs a 2-D cross-section mesh of 3-node triangles" "$WORK/ofe.log"; then
	echo "  *** it was rejected for a different reason than the element order" >&2; status=1
fi
awk '/^convection/{print; print "radiation = 0.8 20"; print "radiationexchange = 1"; next} {print}' "$SRC/thermal_f2d.ofe" > "$WORK/f2x_rx.ofe"
mesh_reject "radiationexchange on a 2-D cross-section" f2x_rx.ofe
if ! grep -q "not implemented on a 2-D cross-section" "$WORK/ofe.log"; then
	echo "  *** it was rejected for a different reason than the 2-D radiation exchange" >&2; status=1
fi
awk '/^conductor = 0/{print "material = 1.0 0"; print "kthermal = 2 400"} /^analysis/{print "convection = 50 20"; print "thermal = 1"} {print}' "$SRC/plate_line_ac.ofe" > "$WORK/f2x_grid.ofe"
mesh_reject "AC heating (F) on the structured grid" f2x_grid.ofe
if ! grep -q "needs a 2-D cross-section mesh of 3-node triangles" "$WORK/ofe.log"; then
	echo "  *** it was rejected for a different reason than the structured grid" >&2; status=1
fi

# ---- 対流の h(T) (convectiontemptable) ----
#
# h(T) = h f(T) (折れ線) で対流が非線形になり、輻射と同じく Picard の各回で Newton
# 線形化する。側面だけから冷やす一様な棒は温度が一様のままなので、釣り合いの
# 代数方程式 (f は同じ折れ線) を二分法で解いた値と機械精度で比べられる。
echo "[thermal h(T)] convection with a temperature-dependent h (convectiontemptable)"
cp "$SRC/thermal_natconv.ofe" "$WORK/"
th_convt() {	# th_convt <I> <α> <表 "T1 f1 T2 f2 ..."> -> "R T" (thermal_conv.ofe の棒)
	awk -v I="$1" -v a="$2" -v tab="$3" 'BEGIN {
		nt = split(tab, v, " ") / 2; for (k = 1; k <= nt; k++) { tt[k] = v[2 * k - 1]; ff[k] = v[2 * k] }
		r0 = 1e-6; T0 = 0; Ta = 30; h = 200; L = 0.01; A = 5e-7; P = 3e-3; J = I / A
		lo = Ta; hi = 5000
		for (it = 0; it < 300; it++) {
			T = (lo + hi) / 2
			if (T <= tt[1]) f = ff[1]; else if (T >= tt[nt]) f = ff[nt]; else { for (k = 1; k < nt; k++) if (T <= tt[k + 1]) { f = ff[k] + (T - tt[k]) / (tt[k + 1] - tt[k]) * (ff[k + 1] - ff[k]); break } }
			F = J * J * r0 * (1 + a * (T - T0)) * A - P * h * f * (T - Ta)
			if (F > 0) lo = T; else hi = T
		}
		printf "%.12e %.12e", r0 * (1 + a * (T - T0)) * L / A, T }'
}
nctab=$(awk '/^convectiontemptable/ { for (i = 3; i <= NF; i++) printf "%s ", $i; exit }' "$SRC/thermal_natconv.ofe")
# (1) 自然対流の表 (thermal_natconv.ofe)。許容 1e-9 : Picard 反復の打ち切りの分
th_run thermal_natconv.ofe
set -- $(th_convt 6 4e-3 "$nctab")
compare "h(T), natural convection, uniform bar : R [ohm]" "$(th_r)" "$1" 1e-9
compare "h(T), natural convection, uniform bar : T [degC]" "$(th_t)" "$2" 1e-9
th_energy "h(T), natural convection"
th_nface "h(T) (sides only)" 160
# (2) f ≡ 1 の表は表なし (線形の対流) と一致する (経路が違うので許容は反復の打ち切りの分)
th_run thermal_conv.ofe
rlin=$(th_r); tlin=$(th_t)
awk '/^convectiontemptable/{print "convectiontemptable = 0 1 1000 1"; next} {print}' \
	"$SRC/thermal_natconv.ofe" > "$WORK/thct1.ofe"
th_run thct1.ofe
compare "h(T) with f = 1 == constant h : R [ohm]" "$(th_r)" "$rlin" 1e-9
compare "h(T) with f = 1 == constant h : T [degC]" "$(th_t)" "$tlin" 1e-9
# (3) α = 0 なら非線形は対流だけ。Newton 線形化なら 6 回 (実測)。接線の f'(T)(T - T_a) を
#     落としても (割線にしても) 同じ不動点に遅く収束するだけなので、回数を assert する
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$SRC/thermal_natconv.ofe" > "$WORK/thct_a0.ofe"
th_run thct_a0.ofe
set -- $(th_convt 6 0 "$nctab")
compare "h(T), alpha = 0 : T [degC]" "$(th_t)" "$2" 1e-9
nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
if [ -n "$nit" ] && [ "$nit" -le 7 ]; then
	echo "  Newton on h(T) (alpha = 0) : converged in $nit iterations (<= 7) -> OK"
else
	echo "  Newton on h(T) (alpha = 0) : ${nit:-?} iterations (> 7) -> NG" >&2
	status=1
fi
# (4) 恒等式 : ヒートシンク + h(T) の対流 + 輻射の首のある導体 (構造格子)、
#     四面体のフィン、過渡 (CN)、面ごとの係数 (下面だけ冷やす板)
awk '/^thermal/{print "convection = 100 30"; print "convectiontemptable = 30 0.5 300 2"; print "radiation = 0.9 10"} {print}' \
	"$SRC/thermal_neck.ofe" > "$WORK/thct_neck.ofe"
th_run thct_neck.ofe
th_energy "h(T) : thermal_neck, sinks + h(T) convection + radiation"
th_nface "h(T) : thermal_neck" 552
awk '/^thermal/{print "convection = 20 30"; print "convectiontemptable = 20 0.5 120 3"} {print}' \
	"$SRC/thermal_tet.ofe" > "$WORK/thct_tet.ofe"
th_run thct_tet.ofe
th_energy "h(T) : tetrahedral fin"
awk '/^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 2 50 0.5"; next} {print}' \
	"$SRC/thermal_natconv.ofe" > "$WORK/thct_tr.ofe"
th_run thct_tr.ofe
th_balance "h(T) : transient (CN)"
awk '/^convection = 500/{print; print "convectiontemptable = 20 0.8 100 1.5"; next} {print}' \
	"$SRC/thermal_plate.ofe" > "$WORK/thct_plate.ofe"
th_run thct_plate.ofe
th_energy "h(T) : per-face convection (bottom face only)"
# (6) 放射率の温度依存 (radiationtemptable、ε(T) = ε f(T))。輻射だけで冷やす一様な棒の
#     釣り合い J²ρ0(1 + α(T - T0))A = P ε f(T) σ (T_K⁴ - T_r,K⁴) と 1e-9。
#     接線は q' = εσ (4 f T³ + f'(T⁴ - T_r⁴))。α = 0 で 9 回 (実測、折れ線の角のぶん
#     表なしの 7 回より多い)。f' の項を落とすと遅くなるので回数を assert する
awk '/^convection/{print "radiation = 0.8 10"; print "radiationtemptable = 0 0.5 300 1 600 1.2"; next}
	/^thermal/{print "thermal = 2"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thrt.ofe"
th_run thrt.ofe
set -- $(awk -v I=2 -v tab="0 0.5 300 1 600 1.2" -v a=4e-3 'BEGIN {
	nt = split(tab, v, " ") / 2; for (k = 1; k <= nt; k++) { tt[k] = v[2 * k - 1]; ff[k] = v[2 * k] }
	r0 = 1e-6; T0 = 0; Tr = 10; e = 0.8; L = 0.01; A = 5e-7; P = 3e-3; J = I / A; sb = 5.670374419e-8; K = 273.15
	lo = Tr; hi = 5000
	for (it = 0; it < 300; it++) {
		T = (lo + hi) / 2
		if (T <= tt[1]) f = ff[1]; else if (T >= tt[nt]) f = ff[nt]; else { for (k = 1; k < nt; k++) if (T <= tt[k + 1]) { f = ff[k] + (T - tt[k]) / (tt[k + 1] - tt[k]) * (ff[k + 1] - ff[k]); break } }
		F = J * J * r0 * (1 + a * (T - T0)) * A - P * e * f * sb * ((T + K) ^ 4 - (Tr + K) ^ 4)
		if (F > 0) lo = T; else hi = T
	}
	printf "%.12e %.12e", r0 * (1 + a * (T - T0)) * L / A, T }')
compare "e(T), radiation only, uniform bar : R [ohm]" "$(th_r)" "$1" 1e-9
compare "e(T), radiation only, uniform bar : T [degC]" "$(th_t)" "$2" 1e-9
th_energy "e(T), radiation only"
awk '/^radiationtemptable/{print "radiationtemptable = 0 1 1000 1"; next} {print}' "$WORK/thrt.ofe" > "$WORK/thrt1.ofe"
awk '!/^radiationtemptable/' "$WORK/thrt.ofe" > "$WORK/thrt0.ofe"
th_run thrt0.ofe
r0=$(th_r)
th_run thrt1.ofe
compare "e(T) with f = 1 == constant emissivity : R [ohm]" "$(th_r)" "$r0" 1e-9
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$WORK/thrt.ofe" > "$WORK/thrt_a0.ofe"
th_run thrt_a0.ofe
nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
if [ -n "$nit" ] && [ "$nit" -le 10 ]; then
	echo "  Newton on e(T) (alpha = 0) : converged in $nit iterations (<= 10) -> OK"
else
	echo "  Newton on e(T) (alpha = 0) : ${nit:-?} iterations (> 10) -> NG" >&2
	status=1
fi
awk '/^radiationtemptable/{print "radiationtemptable = 0 1.3 100 1"; next} {print}' "$WORK/thrt.ofe" > "$WORK/thrtx_gt1.ofe"
th_reject "radiationtemptable with emissivity x factor > 1" thrtx_gt1.ofe "exceeds 1"
awk '/^radiation =/{print "convection = 200 30"; next} {print}' "$WORK/thrt.ofe" > "$WORK/thrtx_norad.ofe"
# (5) 弾くべき入力と警告
awk '/^convectiontemptable/{print "convectiontemptable = 100 1 50 2"; next} {print}' \
	"$SRC/thermal_natconv.ofe" > "$WORK/thctx_dec.ofe"
th_reject "convectiontemptable with decreasing temperatures" thctx_dec.ofe "invalid convectiontemptable data"
awk '/^convectiontemptable/{print "convectiontemptable = 30 1 100"; next} {print}' \
	"$SRC/thermal_natconv.ofe" > "$WORK/thctx_odd.ofe"
th_reject "convectiontemptable with an odd number of values" thctx_odd.ofe "invalid convectiontemptable data"
awk '/^convectiontemptable/{print "convectiontemptable = 30 -1 100 1"; next} {print}' \
	"$SRC/thermal_natconv.ofe" > "$WORK/thctx_neg.ofe"
th_reject "convectiontemptable with a negative factor" thctx_neg.ofe "invalid convectiontemptable data"

# ---- 位置に依存する h / ε (convectionprofile / radiationprofile) ----
#
# h(r) = h f(r の軸成分) の折れ線 (範囲の外は端の値)。面の中心で引く (直方体で面を選ぶのと
# 同じ点。thermal.c の surf_bc)。
#   (a) 一定の倍率 f = 2 で h = 100 は h = 200 と同じ演算 -> ofe.out がバイト一致
#   (b) 急な段差 (x = L/2 の前後 1e-8 で 0.5 -> 2) は「左半分 h/2・右半分 2h の直方体 2 行」と
#       面ごとの係数が同じ -> ofe.out がバイト一致 (構造格子・六面体・四面体。輻射も)。
#       軸の取り違え・面の中心でない点で引く誤りはここで落ちる (面の境目が段差に乗らないので)
#   (c) 連続な 1 次の倍率 (f = 0.5 -> 1.5) を熱伝導の大きい棒 (k = 1e5、mL = 0.035) で :
#       温度はほぼ一様で T - T_a = I²R / (h P L <f>)、<f> = 1 (1 次の f は中点則で厳密に積分
#       される)。一様からのずれは f の勾配と (mL)² の積の 2 次で 1e-5 程度 (実測 1.0e-5) なので許容 1e-4。
#       k = 1e7 にすると系が悪条件になり、恒等式が CG の残差 (2e-8) で落ちる (実測)
#   (d) 鏡映 : f を x について裏返しても棒は左右対称なので R・T_max が一致 (1e-9)
echo "[thermal h(x)] position-dependent convection / emissivity (convectionprofile / radiationprofile)"
cp "$SRC/thermal_conv.ofe" "$SRC/thermal_hex.ofe" "$SRC/thermal_tet.ofe" "$SRC/thermal_hex40.msh" "$SRC/thermal_tet40.msh" "$WORK/"
pf_same() {	# pf_same <ラベル> <ofe a> <ofe b> : 2 つの入力の ofe.out がバイト一致すること
	rm -rf "$WORK/pfa" "$WORK/pfb"; mkdir -p "$WORK/pfa" "$WORK/pfb"
	cp "$WORK/$2" "$WORK/pfa/"; cp "$WORK/$3" "$WORK/pfb/"
	cp "$WORK"/thermal_hex40.msh "$WORK"/thermal_tet40.msh "$WORK/pfa/"; cp "$WORK"/thermal_hex40.msh "$WORK"/thermal_tet40.msh "$WORK/pfb/"
	(cd "$WORK/pfa" && "$OFE" -n 2 "$2" > /dev/null 2>&1) || true
	(cd "$WORK/pfb" && "$OFE" -n 2 "$3" > /dev/null 2>&1) || true
	if [ -s "$WORK/pfa/ofe.out" ] && cmp -s "$WORK/pfa/ofe.out" "$WORK/pfb/ofe.out"; then
		echo "  $1 : ofe.out identical -> OK"
	else
		echo "  $1 : ofe.out differs (or missing) -> NG" >&2; status=1
	fi
}
# (a)
awk '/^convection/{print "convection = 100 30"; print "convectionprofile = x 0 2 1 2"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_a2.ofe"
pf_same "constant profile f = 2 x h = 100 == h = 200" thermal_conv.ofe thpf_a2.ofe
# (b) 構造格子 (対流 / 輻射)
awk '/^convection/{print "convection = 400 30"; print "convectionprofile = x 0.00499999 0.5 0.00500001 2"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_s.ofe"
awk '/^convection/{print "convection = 200 30"; print "convection = 800 30 0.005 0.01 0 1e-3 0 0.5e-3"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_b.ofe"
pf_same "step profile == two boxes (structured grid, convection)" thpf_s.ofe thpf_b.ofe
awk '/^convection/{print "radiation = 0.8 30"; print "radiationprofile = x 0.00499999 0.5 0.00500001 1"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_rs.ofe"
awk '/^convection/{print "radiation = 0.4 30"; print "radiation = 0.8 30 0.005 0.01 0 1e-3 0 0.5e-3"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_rb.ofe"
pf_same "step profile == two boxes (structured grid, radiation)" thpf_rs.ofe thpf_rb.ofe
#     y 向きの 3 段 (棒の幅方向) : y = 0 の面 0.5、y = 1e-3 の面 2、z の面 (中心 y = 0.5e-3) は
#     平らな区間の 1.25 (傾き 0 なので丸めずに 1.25)。直方体 3 行 (500 / 200 / 800) と同じ演算
awk '/^convection/{print "convection = 400 30"; print "convectionprofile = y 0.2e-3 0.5 0.3e-3 1.25 0.7e-3 1.25 0.8e-3 2"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_y.ofe"
awk '/^convection/{print "convection = 500 30"; print "convection = 200 30 0 0.01 0 0 0 0.5e-3"; print "convection = 800 30 0 0.01 1e-3 1e-3 0 0.5e-3"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_yb.ofe"
pf_same "y profile (three levels) == three boxes" thpf_y.ofe thpf_yb.ofe
# (b) 六面体・四面体
for m in thermal_hex thermal_tet; do
	awk '/^heatsink = 1/{next} /^thermal/{print "convection = 400 30"; print "convectionprofile = x 0.00499999 0.5 0.00500001 2"} {print}' \
		"$SRC/$m.ofe" > "$WORK/thpf_${m}_s.ofe"
	awk '/^heatsink = 1/{next} /^thermal/{print "convection = 200 30"; print "convection = 800 30 0.005 0.01 -1 1 -1 1"} {print}' \
		"$SRC/$m.ofe" > "$WORK/thpf_${m}_b.ofe"
	pf_same "step profile == two boxes ($m)" "thpf_${m}_s.ofe" "thpf_${m}_b.ofe"
done
# (c) 連続な 1 次の倍率、熱伝導の大きい棒 (α = 0)
awk '/^tempco/{print "tempco = 2 0 0"; next} /^kthermal/{print "kthermal = 2 1e5"; next}
	/^convection/{print "convection = 200 30"; print "convectionprofile = x 0 0.5 0.01 1.5"; next} {print}' \
	"$SRC/thermal_conv.ofe" > "$WORK/thpf_c.ofe"
th_run thpf_c.ofe
#     R = 1e-6 x 0.01 / 5e-7 = 0.02、P = 6² x 0.02 = 0.72 W、h P L <f> = 200 x 3e-3 x 0.01 = 6e-3
compare "linear profile, high k : T_mean - T_a = I^2 R / (h P L <f>) [K]" \
	"$(awk -v t="$(t3_tm)" 'BEGIN { printf "%.12e", t - 30 }')" 120 1e-4
th_energy "linear profile, high k"
grep -q "convection h(r) = h f(x), f : (0, 0.5) (0.01, 1.5) (at the face centre)" "$WORK/ofe.log" \
	&& echo "  the profile is echoed to the log -> OK" || { echo "  the profile is not in the log -> NG" >&2; status=1; }
# (d) 鏡映
awk '/^tempco/{print "tempco = 2 4e-3 0"; next} {print}' "$WORK/thpf_c.ofe" | awk '/^kthermal/{print "kthermal = 2 10"; next} {print}' > "$WORK/thpf_m1.ofe"
awk '/^convectionprofile/{print "convectionprofile = x 0 1.5 0.01 0.5"; next} {print}' "$WORK/thpf_m1.ofe" > "$WORK/thpf_m2.ofe"
th_run thpf_m1.ofe; r1=$(th_r); t1=$(th_t)
th_run thpf_m2.ofe
compare "mirrored profile : R (mirror symmetry) [ohm]" "$(th_r)" "$r1" 1e-9
compare "mirrored profile : T_max (mirror symmetry) [degC]" "$(th_t)" "$t1" 1e-9
#     対照 : 1 次の倍率は一定の h と違う答えになること (倍率が効いていない誤りの対)
th_run thermal_conv.ofe
res=$(awk -v a="$t1" -v b="$(th_t)" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; printf "%s (%.2e)", (d > 1e-3) ? "OK" : "NG", d }')
echo "  a linear profile changes T_max against a constant h (control) : $res"
case "$res" in NG*) status=1 ;; esac
# (e) 弾くべき入力
awk '/^convection/{print; print "convectionprofile = w 0 1 1 2"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thpfx_ax.ofe"
th_reject "convectionprofile with an unknown axis" thpfx_ax.ofe "invalid convectionprofile data"
awk '/^convection/{print; print "convectionprofile = x 0.01 1 0 2"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thpfx_dec.ofe"
th_reject "convectionprofile with decreasing positions" thpfx_dec.ofe "invalid convectionprofile data"
awk '/^convection/{print "radiation = 0.8 30"; print "radiationprofile = x 0 1 0.01 1.5"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thpfx_eps.ofe"
th_reject "radiationprofile with an emissivity above 1" thpfx_eps.ofe "exceeds 1"
awk '/^convection/{print; print "convectionprofile = x 0 0 1 0"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thpfx_zero.ofe"
th_reject "convectionprofile with all factors zero" thpfx_zero.ofe "all factors are zero"

# ---- 位置の式の h / ε (convectionexpr / radiationexpr、sol/expr.c) ----
#
# 面の中心 (x, y, z) で式 f を評価して h (ε) に掛ける (折れ線もあれば両方の積)。
#   (a) f = 1 と「折れ線 2 x 式 0.5」は式なしと同じ演算 -> ofe.out がバイト一致
#   (b) step() の式は 2 進で正確な倍率 (0.5 / 1.25 / 2) になるので、直方体の行とバイト一致
#       (x の段差 : 構造格子・六面体・四面体、y の 3 段、輻射、演算子の意味)。変数の取り違えはここで落ちる
#   (c) 滑らかな式 1 + 0.5 sin(300x + 0.3) を、面の中心 (x = (i + 1/2) dx) を節点にして awk で同じ式の値を
#       入れた折れ線と比べる。節点では折れ線が値そのものなので R・T_max が丸めの範囲で一致
#   (d) 構文の誤り (入力エラー) と、面で負になる式・ε f > 1 (実行時のエラー)
echo "[thermal h(r) expression] arbitrary position-dependent convection / emissivity (convectionexpr / radiationexpr)"
awk '/^convection/{print; print "convectionexpr = 1"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_1.ofe"
pf_same "expression f = 1 == no expression" thermal_conv.ofe thex_1.ofe
awk '/^convection/{print; print "convectionprofile = x 0 2 1 2"; print "convectionexpr = 0.5"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_pe.ofe"
pf_same "profile 2 x expression 0.5 == no profile / expression" thermal_conv.ofe thex_pe.ofe
# (b)
awk '/^convection/{print "convection = 400 30"; print "convectionexpr = 0.5 + 1.5 * step(x - 0.005)"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_s.ofe"
pf_same "step(x) expression == two boxes (structured grid, convection)" thex_s.ofe thpf_b.ofe
#     演算子の意味 (^ は右結合 : 2^3^2 = 2^9 = 512、左結合なら 64)・単項の負号・min / max / abs・割り算を
#     2 進で正確に同じ倍率になる形で : 2^3^2/1024 = 0.5、max(min(3/2, 7) step, -abs(-1)) = 1.5 step
awk '/^convection/{print "convection = 400 30"; print "convectionexpr = 2^3^2 / 1024 + max(min(3 / 2, 7) * step(x - 0.005), -abs(-1))"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_op.ofe"
pf_same "operators (right-associative ^, unary minus, min / max / abs) == two boxes" thex_op.ofe thpf_b.ofe
awk '/^convection/{print "radiation = 0.8 30"; print "radiationexpr = 0.5 + 0.5*step(x-0.005)"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_rs.ofe"
pf_same "step(x) expression == two boxes (structured grid, radiation)" thex_rs.ofe thpf_rb.ofe
awk '/^convection/{print "convection = 400 30"; print "convectionexpr = 0.5 + 0.75 * step(y - 0.4e-3) + 0.75 * step(y - 0.6e-3)"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_y.ofe"
pf_same "y expression (three levels) == three boxes" thex_y.ofe thpf_yb.ofe
for m in thermal_hex thermal_tet; do
	awk '/^heatsink = 1/{next} /^thermal/{print "convection = 400 30"; print "convectionexpr = 0.5 + 1.5 * step(x - 0.005)"} {print}' \
		"$SRC/$m.ofe" > "$WORK/thex_${m}_s.ofe"
	pf_same "step(x) expression == two boxes ($m)" "thex_${m}_s.ofe" "thpf_${m}_b.ofe"
done
# (c) x を 20 セルにして面の中心を 20 点 (折れ線の上限 32 点に収める)
awk '/^xmesh/{print "xmesh = 0 20 0.01"; next} /^convection/{print; print "convectionexpr = 1 + 0.5 * sin(300 * x + 0.3)"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_c.ofe"
awk 'BEGIN { s = "convectionprofile = x"; for (i = 0; i < 20; i++) { x = (i + 0.5) * (0.01 / 20); s = s sprintf(" %.17g %.17g", x, 1 + 0.5 * sin(300 * x + 0.3)) } print s }' > "$WORK/thex_pf.txt"
awk -v pf="$WORK/thex_pf.txt" '/^xmesh/{print "xmesh = 0 20 0.01"; next} /^convection/{print; while ((getline l < pf) > 0) print l; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thex_cp.ofe"
th_run thex_cp.ofe; r1=$(th_r); t1=$(th_t)
th_run thex_c.ofe
compare "smooth expression == profile through the face centres : R [ohm]" "$(th_r)" "$r1" 1e-12
compare "smooth expression == profile through the face centres : T_max [degC]" "$(th_t)" "$t1" 1e-12
th_energy "smooth expression"
grep -q "convection h(r) = h x f(x, y, z), f = 1 + 0.5 \* sin(300 \* x + 0.3) (at the face centre)" "$WORK/ofe.log" \
	&& echo "  the expression is echoed to the log -> OK" || { echo "  the expression is not in the log -> NG" >&2; status=1; }
th_run thermal_conv.ofe
res=$(awk -v a="$t1" -v b="$(th_t)" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; printf "%s (%.2e)", (d > 1e-3) ? "OK" : "NG", d }')
echo "  the smooth expression changes T_max against a constant h (control) : $res"
case "$res" in NG*) status=1 ;; esac
# (d)
for bad in "1 +" "sin(x" "foo(x)" "x y" "min(x)" "2 * (x + 1))"; do
	awk -v e="$bad" '/^convection/{print; print "convectionexpr = " e; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thexx.ofe"
	th_reject "convectionexpr '$bad'" thexx.ofe "invalid convectionexpr"
done
awk '/^convection/{print; print "convectionexpr = 1"; print "convectionexpr = 2"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thexx2.ofe"
th_reject "convectionexpr given twice" thexx2.ofe "invalid convectionexpr"
ex_runerr() {	# ex_runerr <ラベル> <ofe> <ログの文字列> : 実行時に止まり、ログに理由が出ること
	if (cd "$WORK" && "$OFE" -n 2 "$2" > /dev/null 2>&1); then
		echo "  $1 : not rejected -> NG"; status=1
	elif grep -q "$3" "$WORK/ofe.log"; then
		echo "  $1 : rejected -> OK"
	else
		echo "  $1 : rejected for a different reason -> NG"; status=1
	fi
}
awk '/^convection/{print; print "convectionexpr = x - 0.005"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thexx_neg.ofe"
ex_runerr "convectionexpr negative on some faces" thexx_neg.ofe "is negative or not finite on 80 exposed faces"
awk '/^convection/{print "radiation = 0.8 30"; print "radiationexpr = 1 + 100 * x"; next} {print}' "$SRC/thermal_conv.ofe" > "$WORK/thexx_eps.ofe"
ex_runerr "radiationexpr with e f > 1" thexx_eps.ofe "radiationexpr = 1 + 100 \* x is negative or not finite (or e f > 1)"

# ---- 温度と位置の両方に依る式 (convectionexpr / radiationexpr に T、sol/expr.c の expr_eval_t) ----
#
# T (面の温度 [degC]) を含む式は h(T) の表と同じ Newton の経路で面の位置・積分点の温度ごとに評価し、df/dT は
# 前進型の自動微分 (2 重数) で厳密に作る。同じ関数を表と式の 2 通りに書いて比べる (閉形式が要らない) :
#   (a) h(T) = h (1 + 0.002 (T − 20)) を直線の表 (−100..2000 の 2 点) と式で : R・T_max が 1e-11、反復回数が一致
#   (b) 位置と温度の積 (0.5 + 1.5 step(x − 0.005)) (1 + 0.002 (T − 20)) を「位置の式 × 温度の表」と 1 つの式で
#   (c) 輻射の ε(T) = 0.2 + 0.004 T を表と式で (傾きが小さいと接線の誤りが反復回数に出ないので急にしてある)
# 接線の q f' を落とす変異は答えを変えず反復回数だけを変える (20 -> 13、22 -> 16、14 -> 12。実測) ので回数を assert する。
#   (d) 弾くべき入力 : 面で負になる T の式 (実行時)、相互放射との併用
echo "[thermal h(T, r) expression] convection / emissivity expressions in the face temperature T"
cp "$SRC/thermal_conv.ofe" "$SRC/thermal_rad.ofe" "$WORK/"
txr() {	# txr <ofe> : "R T_max 反復回数"
	(cd "$WORK" && "$OFE" -n 2 "$1" > /dev/null 2>&1) || { echo "x x x"; return; }
	awk '/^  R = / && !r { r = $3 } /^  T_max = / && !t { t = $3 } /converged in/ { n = $(NF - 1) } END { printf "%s %s %s", r, t, n }' "$WORK/ofe.log"
}
tx_same() {	# tx_same <ラベル> <表の ofe> <式の ofe>
	set -- "$1" $(txr "$2") $(txr "$3")
	res=$(awk -v r1="$2" -v t1="$3" -v n1="$4" -v r2="$5" -v t2="$6" -v n2="$7" 'BEGIN { if ((r1 == "x") || (r2 == "x")) { print "NG (a run failed)"; exit }
		dr = r2 / r1 - 1; dt = t2 / t1 - 1; if (dr < 0) dr = -dr; if (dt < 0) dt = -dt
		ok = ((dr < 1e-11) && (dt < 1e-11) && (n1 == n2) && (n1 > 1)) ? "OK" : "NG"
		printf "%s (R %.1e, T_max %.1e, %d / %d iterations)", ok, dr, dt, n1, n2 }')
	echo "  $1 : table == expression in T : $res"
	case "$res" in NG*) status=1 ;; esac
}
awk '/^convection/{print; print "convectiontemptable = -100 0.76  2000 4.96"; next} {print}' "$WORK/thermal_conv.ofe" > "$WORK/txa_t.ofe"
awk '/^convection/{print; print "convectionexpr = 1 + 0.002 * (T - 20)"; next} {print}' "$WORK/thermal_conv.ofe" > "$WORK/txa_e.ofe"
tx_same "(a) h(T) = h (1 + 0.002 (T - 20))" txa_t.ofe txa_e.ofe
awk '/^convection/{print; print "convectiontemptable = -100 0.76  2000 4.96"; print "convectionexpr = 0.5 + 1.5 * step(x - 0.005)"; next} {print}' "$WORK/thermal_conv.ofe" > "$WORK/txb_t.ofe"
awk '/^convection/{print; print "convectionexpr = (0.5 + 1.5 * step(x - 0.005)) * (1 + 0.002 * (T - 20))"; next} {print}' "$WORK/thermal_conv.ofe" > "$WORK/txb_e.ofe"
tx_same "(b) h(x, T) = h (0.5 + 1.5 step(x - 0.005)) (1 + 0.002 (T - 20))" txb_t.ofe txb_e.ofe
awk '/^radiation =/{print; print "radiationtemptable = 0 0.2  200 1.0"; next} {print}' "$WORK/thermal_rad.ofe" > "$WORK/txc_t.ofe"
awk '/^radiation =/{print; print "radiationexpr = 0.2 + 0.004 * T"; next} {print}' "$WORK/thermal_rad.ofe" > "$WORK/txc_e.ofe"
tx_same "(c) e(T) = 0.2 + 0.004 T" txc_t.ofe txc_e.ofe
awk '/^convection/{print; print "convectionexpr = (T - 1000) / 1000"; next} {print}' "$WORK/thermal_conv.ofe" > "$WORK/txx_neg.ofe"
ex_runerr "convectionexpr in T that is negative on the faces (runtime)" txx_neg.ofe "is negative or not finite at 640 surface points"
awk '/^radiation =/{print; print "radiationexpr = 0.2 + 0.004 * T"; print "radiationexchange = 1"; next} {print}' "$WORK/thermal_rad.ofe" > "$WORK/txx_rx.ofe"
th_reject "radiationexpr in T with radiationexchange" txx_rx.ofe "cannot be combined with a radiationexpr that depends on T"

# ---- 面どうしの相互放射 (radiationexchange、構造格子) ----
#
# 形態係数は閉形式 (周回積分 A_i F_ij = (1/2π) ∮∮ ln r dl_i・dl_j の平行な辺の対)、
# 放射度は密に一度だけ解く (thermal.c の rx_build の注釈)。検査は教科書の閉形式
# (平行な長方形・辺を共有する直交長方形) と、相互性から厳密に決まる面積の比。
echo "[thermal radiation exchange] gray-diffuse surface-to-surface radiation (view factors)"
cp "$SRC/thermal_cavity.ofe" "$SRC/thermal_rad.ofe" "$WORK/"
rx_row() { awk -v r="$1" '$0 ~ ("radiation exchange : row " r ",") { print $NF; exit }' "$WORK/ofe.log"; }
rx_sum() { awk '/sum A F_inf \/ sum A/ { print $NF; exit }' "$WORK/ofe.log"; }
set -- $(awk 'BEGIN { pi = atan2(0, -1); X = 1; Y = 1
	F = (2 / (pi * X * Y)) * (log(sqrt((1 + X * X) * (1 + Y * Y) / (1 + X * X + Y * Y))) + X * sqrt(1 + Y * Y) * atan2(X, sqrt(1 + Y * Y)) + Y * sqrt(1 + X * X) * atan2(Y, sqrt(1 + X * X)) - X * atan2(X, 1) - Y * atan2(Y, 1))
	W = 1; H = 1
	G = (1 / (pi * W)) * (W * atan2(1, W) + H * atan2(1, H) - sqrt(H * H + W * W) * atan2(1, sqrt(H * H + W * W)) + 0.25 * log(((1 + W * W) * (1 + H * H) / (1 + W * W + H * H)) * ((W * W * (1 + W * W + H * H) / ((1 + W * W) * (W * W + H * H))) ^ (W * W)) * ((H * H * (1 + H * H + W * W) / ((1 + H * H) * (H * H + W * W))) ^ (H * H))))
	printf "%.15e %.15e", F, G }')
fpar=$1; fper=$2
# (1) 開いた立方体の空洞 (thermal_cavity.ofe) : 行ごとの F∞ が閉形式に機械精度で一致
th_run thermal_cavity.ofe
compare "cavity floor -> opening == parallel squares (Hottel) [-]" "$(rx_row 2)" "$fpar" 1e-10
compare "cavity wall -> opening == perpendicular squares, common edge [-]" "$(rx_row 3)" "$fper" 1e-10
compare "cavity : sum A F_inf / sum A == opening / wall area = 0.2 [-]" "$(rx_sum)" 0.2 1e-10
th_energy "radiation exchange, open cavity"
tex=$(th_t)
awk '!/^radiationexchange/' "$SRC/thermal_cavity.ofe" > "$WORK/thrx0.ofe"
th_run thrx0.ofe
res=$(awk -v a="$tex" -v b="$(th_t)" 'BEGIN { printf "%s (T_max %.6f with exchange vs %.6f without)", (a > b) ? "OK" : "NG", a, b }')
echo "  exchange inside the cavity reduces the heat lost to the ambient : $res"
case "$res" in NG*) status=1 ;; esac
# (2) 壁で隔てた 2 つの空洞 : 空洞 1 の奥の壁 (x = 2mm) は空洞 2 の奥の壁と向き合うが、
#     間の壁で遮られるので F∞ は 1 つの空洞のときと同じ閉形式のまま (遮蔽の検査)。
#     遮蔽を見ないと F∞ が減る。面積の比 12/44 も相互性から厳密
awk '/^geometry = 0 1 4e-3/{print "geometry = 0 1 2e-3 4e-3 0.5e-3 2.5e-3 1e-3 3e-3"; print "geometry = 0 1 6e-3 8e-3 0.5e-3 2.5e-3 1e-3 3e-3"; next}
	/^radiation = /{ if (!r) { print "radiation = 0.8 10  2e-3 8e-3  0.5e-3 2.5e-3  1e-3 3e-3"; print "radiation = 0.8 10  2e-3 2e-3  0.5e-3 2.5e-3  1e-3 3e-3" } r = 1; next } {print}' \
	"$SRC/thermal_cavity.ofe" > "$WORK/thrx2.ofe"
th_run thrx2.ofe
compare "two cavities : far wall -> opening unchanged (occlusion by the wall) [-]" "$(rx_row 2)" "$fper" 1e-10
compare "two cavities : sum A F_inf / sum A == 12/44 [-]" "$(rx_sum)" "$(awk 'BEGIN { printf "%.15e", 12 / 44 }')" 1e-10
th_energy "radiation exchange, two cavities"
# (3) 凸な棒 (thermal_rad.ofe) : 見える面の対が無いので、通常の輻射と全桁一致
th_run thermal_rad.ofe
rr=$(th_r); tr0=$(th_t)
awk '/^radiation /{print; print "radiationexchange = 1"; next} {print}' "$SRC/thermal_rad.ofe" > "$WORK/thrx_cv.ofe"
th_run thrx_cv.ofe
compare "convex bar with exchange == plain radiation : R [ohm]" "$(th_r)" "$rr" 1e-12
compare "convex bar with exchange == plain radiation : T_max [degC]" "$(th_t)" "$tr0" 1e-12
# (4) 閉じた空洞 (導体の中の空隙) : Σ_j F_ij = 1 なので周囲へ逃げる熱は 0
#     (内部で熱を運ぶだけ)。許容 1e-10 (I²R に対して)
awk '/^geometry = 0 1 4e-3/{print "geometry = 0 1 4e-3 6e-3 0.5e-3 2.5e-3 0.5e-3 2.2e-3"; next}
	/^radiation = /{ if (!r) print "radiation = 0.8 10  4e-3 6e-3  0.5e-3 2.5e-3  0.5e-3 2.2e-3"; r = 1; next } {print}' \
	"$SRC/thermal_cavity.ofe" > "$WORK/thrx_void.ofe"
th_run thrx_void.ofe
res=$(awk '/joule power/ { for (i = 1; i <= NF; i++) { if ($i == "R") a = $(i+2); if ($i == "radiation") g = $(i+2) + 0 } }
	END { d = g / a; if (d < 0) d = -d; printf "%s (to radiation / I^2 R = %.1e)", (d < 1e-10) ? "OK" : "NG", d }' "$WORK/ofe.log")
echo "  closed void : no heat leaves by radiation : $res"
case "$res" in NG*) status=1 ;; esac
compare "closed void : sum A F_inf / sum A == 0 [-]" "$(awk -v s="$(rx_sum)" 'BEGIN { printf "%.15e", 1 + s }')" 1 1e-12
th_energy "radiation exchange, closed void"
# (4b) 熱平衡 : 電流 0 (waveform = 0 0) で全体が T_r なら、どんな形・ε でも q = 0 なので
#      温度は T_r のまま動かない。**反射の項 (1 - ε) の取り違えはこれでしか落ちない**
#      (形態係数の検査は幾何だけ、閉じた空洞の恒等式は放射度の行列に依らない)
awk '/^radiation = /{ sub(/ 10 /, " 20 ") } /^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 0.5 10 0.5"; print "waveform = 0 0"; next} {print}' \
	"$SRC/thermal_cavity.ofe" > "$WORK/thrx_eq.ofe"
th_run thrx_eq.ofe
#      (T_max はヒートシンクの 20 degC に張り付くので見ない。平均温度で見る)
compare "radiation exchange, equilibrium (I = 0, all at T_r) : T_mean stays 20 degC [degC]" \
	"$(awk '/T_mean/ { print $8; exit }' "$WORK/ofe.log")" 20 1e-12
# (4c) α = 0 では対角を Newton、他の面の分を Picard で 5 回 (実測)。上限を assert する
#      (接線に Q_ii でなく ε を使う誤りは Q_ii = ε(1 - O(F²)) と近いので回数も変わらず、
#      不動点も同じ。検出できない等価な変異として README に記す)
awk '/^tempco/{print "tempco = 2 0 0"; next} {print}' "$SRC/thermal_cavity.ofe" > "$WORK/thrx_a0.ofe"
th_run thrx_a0.ofe
nit=$(awk '/converged in/ { print $3; exit }' "$WORK/ofe.log")
if [ -n "$nit" ] && [ "$nit" -le 6 ]; then
	echo "  radiation exchange (alpha = 0) : converged in $nit iterations (<= 6) -> OK"
else
	echo "  radiation exchange (alpha = 0) : ${nit:-?} iterations (> 6) -> NG" >&2
	status=1
fi
# (5) 過渡 (CN) の熱の収支
awk '/^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 0.5 20 0.5"; next} {print}' \
	"$SRC/thermal_cavity.ofe" > "$WORK/thrx_tr.ofe"
th_run thrx_tr.ofe
th_balance "radiation exchange, transient (CN)"
# (7) 非構造格子 (thermal.c の rx_build_tet の注釈)。六面体の空洞は面が軸に平行なので、
#     辺の対は平行 (閉形式) か直交 (寄与 0) だけになり、構造格子と同じ離散問題になる :
#     同じ格子点の構造格子 (x を 20 分割) と R・T_max・T_mean が全桁一致 (空洞・2 つの空洞・
#     閉じた空隙)。四面体の空洞は斜めの辺の対 (数値積分) を通り、行ごとの F∞ が
#     閉形式に 1e-10 で一致する (三角形に割っても面の和は同じ形態係数)
for m in thermal_cavity_hex thermal_cavity_tet thermal_cavity2_hex thermal_void_hex; do cp "$SRC/$m.msh" "$WORK/"; done
rx_unst() {	# rx_unst <格子> <構造格子の入力> <出力> : 同じ形の非構造格子の入力にする
	awk -v m="$1" '/^xmesh|^ymesh|^zmesh|^geometry|^conductor/{next}
		/^material/{print "mesh = " m ".msh"; print; print "region = 1 2"; print "electrode = 10 0"; print "electrode = 11 1"; next}
		{print}' "$WORK/$2" > "$WORK/$3"
}
rx_x20() { awk '/^xmesh/{print "xmesh = 0 20 0.01"; next} {print}' "$WORK/$1" > "$WORK/$2"; }
rx_same() {	# rx_same <ラベル> <構造格子の入力> <六面体の格子>
	rx_x20 "$2" rxg.ofe
	th_run rxg.ofe; r0=$(th_r); t0=$(th_t); m0=$(t3_tm); f0=$(rx_sum)
	rx_unst "$3" rxg.ofe rxh.ofe
	th_run rxh.ofe
	compare "$1 : hexahedra == structured grid : R [ohm]" "$(th_r)" "$r0" 1e-10
	compare "$1 : hexahedra == structured grid : T_max [degC]" "$(th_t)" "$t0" 1e-10
	compare "$1 : hexahedra == structured grid : T_mean [degC]" "$(t3_tm)" "$m0" 1e-10
	compare "$1 : hexahedra == structured grid : sum A F_inf / sum A [-]" "$(awk -v s="$(rx_sum)" 'BEGIN { printf "%.15e", 1 + s }')" \
		"$(awk -v s="$f0" 'BEGIN { printf "%.15e", 1 + s }')" 1e-12
	th_energy "$1, hexahedra"
}
cp "$WORK/thrx2.ofe" "$WORK/thrx2s.ofe"
rx_same "open cavity" thermal_cavity.ofe thermal_cavity_hex
compare "hexahedra : cavity floor -> opening == parallel squares [-]" "$(rx_row 2)" "$fpar" 1e-10
rx_same "two cavities (occlusion by the wall)" thrx2s.ofe thermal_cavity2_hex
compare "hexahedra : two cavities, far wall -> opening (occluded) [-]" "$(rx_row 2)" "$fper" 1e-10
rx_same "closed void" thrx_void.ofe thermal_void_hex
res=$(awk '/joule power/ { for (i = 1; i <= NF; i++) { if ($i == "R") a = $(i+2); if ($i == "radiation") g = $(i+2) + 0 } }
	END { d = g / a; if (d < 0) d = -d; printf "%s (to radiation / I^2 R = %.1e)", (d < 1e-10) ? "OK" : "NG", d }' "$WORK/ofe.log")
echo "  hexahedra : closed void, no heat leaves by radiation : $res"
case "$res" in NG*) status=1 ;; esac
#     四面体 : 斜めの辺を含む三角形の面
rx_x20 thermal_cavity.ofe rxg.ofe
rx_unst thermal_cavity_tet rxg.ofe rxt.ofe
th_run rxt.ofe
compare "tetrahedra : cavity floor -> opening == parallel squares (Hottel) [-]" "$(rx_row 2)" "$fpar" 1e-10
compare "tetrahedra : cavity wall -> opening == perpendicular squares [-]" "$(rx_row 3)" "$fper" 1e-10
compare "tetrahedra : sum A F_inf / sum A == 0.2 [-]" "$(rx_sum)" 0.2 1e-10
th_energy "radiation exchange, tetrahedra"
tex=$(th_t)
awk '!/^radiationexchange/' "$WORK/rxt.ofe" > "$WORK/rxt0.ofe"
th_run rxt0.ofe
res=$(awk -v a="$tex" -v b="$(th_t)" 'BEGIN { printf "%s (T_max %.6f with exchange vs %.6f without)", (a > b) ? "OK" : "NG", a, b }')
echo "  tetrahedra : exchange inside the cavity reduces the heat lost to the ambient : $res"
case "$res" in NG*) status=1 ;; esac
awk '/^thermal/{print; print "heatcap = 2 3.4e6"; print "transient = 0.5 20 0.5"; next} {print}' \
	"$WORK/rxt.ofe" > "$WORK/rxt_tr.ofe"
th_run rxt_tr.ofe
th_balance "radiation exchange, tetrahedra, transient (CN)"
#     2 値の判定 (向き合うか・遮られるか) が丸めで割れないこと : 節点の座標を 1 ulp 程度
#     揺らした格子でも同じ答え (同一平面の対の n・d = ±1e-20 や、共有辺を通る線分。
#     どちらも CI の Linux / macOS でだけ落ちたことがある)
rx_pert() {	# rx_pert <格子> : 座標に相対 2e-16 の乱数を掛けた <格子>_p.msh
	awk 'BEGIN { srand(11) } /^\$Nodes/ { n = 1; print; getline; print; next } /^\$EndNodes/ { n = 0 }
		n { printf "%d %.17g %.17g %.17g\n", $1, $2 * (1 + (rand() - 0.5) * 4e-16), $3 * (1 + (rand() - 0.5) * 4e-16), $4 * (1 + (rand() - 0.5) * 4e-16); next }
		{ print }' "$WORK/$1.msh" > "$WORK/$1_p.msh"
}
rx_pert thermal_cavity_tet
awk '/^mesh/{print "mesh = thermal_cavity_tet_p.msh"; next} {print}' "$WORK/rxt.ofe" > "$WORK/rxtp.ofe"
th_run rxtp.ofe
compare "tetrahedra, perturbed by 1 ulp : cavity floor -> opening (Hottel) [-]" "$(rx_row 2)" "$fpar" 1e-10
compare "tetrahedra, perturbed by 1 ulp : sum A F_inf / sum A == 0.2 [-]" "$(rx_sum)" 0.2 1e-10
rx_pert thermal_cavity2_hex
rx_x20 thrx2s.ofe rxg.ofe
th_run rxg.ofe; f0=$(rx_sum)
rx_unst thermal_cavity2_hex_p rxg.ofe rxhp.ofe
th_run rxhp.ofe
compare "two cavities, hexahedra perturbed by 1 ulp == structured grid : sum A F_inf / sum A [-]" \
	"$(awk -v s="$(rx_sum)" 'BEGIN { printf "%.15e", 1 + s }')" "$(awk -v s="$f0" 'BEGIN { printf "%.15e", 1 + s }')" 1e-12
# (8) 高次の四面体の曲がった面 (thermal.c の rx_build_tet)。境界は辺を 2p 分割した折れ線 (周回積分は境界だけで
#     決まる)、面積は曲がった面の積分、凹の面は自分を見る (F_ii を面の積分点の対で直接積分)、遮蔽物も曲がった面
#     ((2p)² 枚の三角形)。
#     (a) 平らな面の 2 次の格子 (thermal_cavity_tet2.msh、gmsh_cavity.py) : 行ごとの F∞ が Hottel に 1e-10
#         (折れ線が直線の上に乗るので周回積分は 1 次と同じ値)
#     (b) 閉じた球の空隙 (thermal_sphere_void.ofe) : Σ A F∞ / Σ A が 2 次 4e-3・3 次はその 0.6 倍以下 (実測 2.7e-3 /
#         1.3e-3)、周囲へ逃げる熱 / I²R が 1e-5 以下、熱の収支の恒等式。F_ii を落とすと 1.2e-2、遮蔽を平らな面で
#         見ると対の 96% が遮られ F∞ = 0.96 (どちらも実測)
cp "$SRC/thermal_cavity_tet2.msh" "$SRC/sphere_void_p2.msh" "$SRC/sphere_void_p3.msh" "$SRC/thermal_sphere_void.ofe" "$WORK/"
rx_x20 thermal_cavity.ofe rxg.ofe
rx_unst thermal_cavity_tet2 rxg.ofe rxt2.ofe
th_run rxt2.ofe
compare "second-order tetrahedra (flat faces) : cavity floor -> opening == parallel squares (Hottel) [-]" "$(rx_row 2)" "$fpar" 1e-10
compare "second-order tetrahedra (flat faces) : cavity wall -> opening == perpendicular squares [-]" "$(rx_row 3)" "$fper" 1e-10
compare "second-order tetrahedra (flat faces) : sum A F_inf / sum A == 0.2 [-]" "$(rx_sum)" 0.2 1e-10
th_energy "radiation exchange, second-order tetrahedra"
rxs_loss() {	# 周囲へ逃げる熱 / I²R
	awk '/joule power/ { for (i = 1; i <= NF; i++) { if ($i == "R") a = $(i+2); if ($i == "radiation") g = $(i+2) + 0 } }
		END { d = g / a; if (d < 0) d = -d; printf "%.3e", d }' "$WORK/ofe.log"
}
for o in 2 3; do
	sed "s/^mesh = .*/mesh = sphere_void_p$o.msh/" "$WORK/thermal_sphere_void.ofe" > "$WORK/rxs$o.ofe"
	th_run "rxs$o.ofe"
	th_energy "radiation exchange, spherical void, order $o"
	eval "rxs_c$o=\$(rx_sum); rxs_l$o=\$(rxs_loss)"
done
res=$(awk -v a="$rxs_c2" -v b="$rxs_c3" -v la="$rxs_l2" -v lb="$rxs_l3" 'function ab(x) { return (x < 0) ? -x : x }
	BEGIN { ok = (a != "") && (ab(a) <= 4e-3) && (ab(b) <= 0.6 * ab(a)) && (la <= 1e-5) && (lb <= 1e-5)
	printf "%s (sum A F_inf / sum A : P2 %.2e, P3 %.2e; to radiation / I^2 R : %.1e, %.1e)", ok ? "OK" : "NG", a, b, la, lb }')
echo "  closed spherical void (curved faces) : the enclosure closes as the order rises, no heat leaves : $res"
case "$res" in NG*) status=1 ;; esac
#     (c) 高次の六面体・角柱の平らな面 (thermal_cavity_hex2s.msh = 20 節点 / thermal_cavity_prism2s.msh = 15 節点、
#         mkmesh.py cavity_hex -nx 5 [-prism 1] を gmsh の setOrder(2) (SecondOrderIncomplete) で 2 次にしたもの) :
#         行ごとの F∞ が Hottel に 1e-10。角柱は底面・天井が三角形の面、壁が四角形の面 (面の角の参照座標の対応
#         elem3d_face_point の誤りはここで落ちる)
cp "$SRC/thermal_cavity_hex2s.msh" "$SRC/thermal_cavity_prism2s.msh" "$WORK/"
rx_x20 thermal_cavity.ofe rxg.ofe
for m in thermal_cavity_hex2s thermal_cavity_prism2s; do
	rx_unst "$m" rxg.ofe rxh2.ofe
	th_run rxh2.ofe
	compare "$m (flat faces) : cavity floor -> opening == parallel squares (Hottel) [-]" "$(rx_row 2)" "$fpar" 1e-10
	compare "$m (flat faces) : cavity wall -> opening == perpendicular squares [-]" "$(rx_row 3)" "$fper" 1e-10
	compare "$m (flat faces) : sum A F_inf / sum A == 0.2 [-]" "$(rx_sum)" 0.2 1e-10
	th_energy "radiation exchange, $m"
done
#     (d) 円環の穴の内面 (曲がった面、両端が開いた円筒) : 側面から両端の開口への形態係数は
#         F∞ = a (1 − F12) / L、F12 = ½ [X − √(X² − 4)]、X = 2 + L²/a² (円板どうし、面積の相互性)。
#         輻射を内面だけ (|x|, |y| < 0.6mm の直方体、外側の面の頂点の平均は 0.68mm 以上) に課す。
#         境界の折れ線 (各辺 2p 分割) の誤差は弧の中心角 θ の分割で (θ / 2p)² / 8 程度 (2 次の 12 分割・3 次の 8 分割で
#         2.1e-3) なので許容は 2e-3。実測 (1 次の六面体の折れ線は −7.0e-3 / −1.5e-2) :
#           2 次 (thermal_ann_h27 / _h27b、12 / 24 分割) −8.4e-4 / −3.0e-4、3 次 (ann_hex3a / 3b、8 / 16 分割)
#           −3.0e-4 / −1.3e-4、4 次 (ann_hex4a、8 分割) 1.5e-6
#         z に積んだ面は中心の法線が中心を結ぶ線と直交し、中心だけでは「同一平面」になるので、面の積分点で
#         向き合うかを見る (rx_vis_curved。それを落とすと積んだ対の F が抜けて 3 次 +3.7e-3、4 次 +4.0e-3、実測)。
#         六面体と角柱 (側面は四角形) は同じ頂点・同じ面なので F∞ が全桁一致する
for m in thermal_ann_h27 thermal_ann_h27b thermal_ann_p18 ann_hex1b ann_hex3a ann_hex3b ann_prism3a ann_prism3b ann_hex4a ann_prism4a; do
	cp "$SRC/$m.msh" "$WORK/"
	case "$m" in thermal_*) lz=2e-4 ;; *) lz=1e-4 ;; esac
	printf 'OpenFEM 1 0\ntitle = radiation exchange inside the hole of an annulus\nmesh = %s.msh\nmaterial = 1.0 1e6\ntempco = 2 0 0\nkthermal = 2 10\nregion = 1 2\nelectrode = 10 0\nelectrode = 11 1\nconvection = 2000 20\nradiation = 0.8 20 -0.6e-3 0.6e-3 -0.6e-3 0.6e-3 0 %s\nradiationexchange = 1\ntemperature = 20\nthermal = 10\nanalysis = R\nsolver = 100000 20000 1e-12\nend\n' "$m" "$lz" > "$WORK/rxan.ofe"
	th_run rxan.ofe
	th_energy "radiation exchange, hole of $m"
	e=$(awk -v s="$(rx_sum)" -v L="$lz" 'BEGIN { if (s == "") { print "nan"; exit }
		a = 0.5e-3; X = 2 + L * L / (a * a); f12 = 0.5 * (X - sqrt(X * X - 4)); f = a * (1 - f12) / L
		printf "%.6e", (s - f) / f }')
	eval "rxa_$m=\$e; rxs_$m=\$(rx_sum)"
	case "$m" in
	ann_hex1b) ;;
	*) res=$(awk -v e="$e" 'BEGIN { a = (e < 0) ? -e : e; printf "%s (%+.2e)", (a <= 2e-3) ? "OK" : "NG", e }')
	   echo "  $m : hole side -> both openings == a (1 - F12) / L : $res"
	   case "$res" in NG*) status=1 ;; esac ;;
	esac
done
res=$(awk -v a2="$rxa_thermal_ann_h27" -v b2="$rxa_thermal_ann_h27b" -v a3="$rxa_ann_hex3a" -v b3="$rxa_ann_hex3b" \
	-v l1="$rxa_ann_hex1b" -v a4="$rxa_ann_hex4a" 'function ab(x) { return (x < 0) ? -x : x }
	BEGIN { ok = (ab(b2) <= 0.5 * ab(a2)) && (ab(b3) <= 0.5 * ab(a3)) && (ab(b3) * 20 <= ab(l1)) && (ab(a4) < ab(a3))
	printf "%s (P2 %.1e -> %.1e, P3 %.1e -> %.1e, P1 facets %.1e, P4 %.1e)", ok ? "OK" : "NG", a2, b2, a3, b3, l1, a4 }')
echo "  hole of an annulus : the curved faces converge under refinement and beat the P1 facets : $res"
case "$res" in NG*) status=1 ;; esac
for p in thermal_ann_h27:thermal_ann_p18 ann_hex3a:ann_prism3a ann_hex3b:ann_prism3b ann_hex4a:ann_prism4a; do
	eval "x=\$rxs_${p%%:*}; y=\$rxs_${p#*:}"
	compare "${p#*:} == ${p%%:*} (same vertices, quadrilateral side faces) : sum A F_inf / sum A [-]" \
		"$(awk -v s="$y" 'BEGIN { printf "%.15e", 1 + s }')" "$(awk -v s="$x" 'BEGIN { printf "%.15e", 1 + s }')" 1e-12
done
#     (e) 係数は面を選んだ点 (頂点の平均) で引く : 直方体の境目が外側の面の頂点の平均 (0.683mm) と曲がった面の
#         中央 (0.707mm) の間にあると、中央で引いた係数は無い (実測 : 異常終了)。走り切って熱の収支が閉じること
sed 's/-0.6e-3 0.6e-3 -0.6e-3 0.6e-3/-0.7e-3 0.7e-3 -0.7e-3 0.7e-3/; s/^mesh = .*/mesh = thermal_ann_h27.msh/; s/ 0 1e-4$/ 0 2e-4/' "$WORK/rxan.ofe" > "$WORK/rxan7.ofe"
th_run rxan7.ofe
th_energy "radiation exchange, box edge between the corner mean and the curved centre"
# (6) 弾くべき入力
awk '/^radiationexchange/{print; print "radiationtemptable = 0 1 500 1.2"; next} {print}' "$SRC/thermal_cavity.ofe" > "$WORK/thrxx_tab.ofe"
th_reject "radiationexchange with radiationtemptable" thrxx_tab.ofe "cannot be combined with radiationtemptable"
awk '/^radiationexchange/{print "radiationexchange = 2"; next} {print}' "$SRC/thermal_cavity.ofe" > "$WORK/thrxx_2.ofe"
th_reject "radiationexchange = 2" thrxx_2.ofe "invalid radiationexchange data"

# ---- 静磁場の過渡解析 (analysis = M + transient : 渦電流 + ヒステリシス) ----
#
# -∇・(ν∇A) + σ(∂A/∂t - c_k) = J_s(t) を後退 Euler で解く (solve.c の mag_transient の注釈)。
# 孤立した導電領域ごとに正味の電流 0 (c_k を射影質量で消去)。刻みに ΔA を掛けると
#   入力 b・ΔA = 磁気の仕事 ΔA・r_ja + 渦電流損 ΔA・M_eff ΔA/Δt
# が反復の打ち切りを除いて厳密 (ログの energy の行)。
echo "[transient M] eddy currents + hysteresis in the time domain (Az, backward Euler)"
cp "$SRC/lamination_step.ofe" "$SRC/plate_line_ja_eddy.ofe" "$SRC/plate_line_ja.ofe" "$WORK/"
mt_b() { awk -F, 'END { print $5 }' "$WORK/ofe_transient.csv"; }
mt_energy() {	# mt_energy <ラベル> <許容> : 入力 = 磁気の仕事 + 渦電流損
	res=$(awk -v tol="$2" '/energy \(per unit length\)/ { for (i = 1; i <= NF; i++) {
			if ($i == "input") a = $(i+2) + 0
			if ($i == "work") b = $(i+2) + 0
			if ($i == "loss") c = $(i+2) + 0 } got = 1 }
		END { if (!got) { printf "NG (no energy line)"; exit }
		      d = (a - b - c) / a; if (d < 0) d = -d
		      printf "%s (%.1e)", (d < tol) ? "OK" : "NG", d }' "$WORK/ofe.log")
	echo "  $1 : input = magnetic work + eddy loss : $res"
	case "$res" in NG*) status=1 ;; esac
}
mt_loss() { awk '/energy \(per unit length\)/ { for (i = 1; i <= NF; i++) if ($i == "loss") print $(i+2) + 0 }' "$WORK/ofe.log"; }
# (1) 線形の積層板のステップ応答 (lamination_step.ofe) : Fourier 級数解と、後退 Euler の
#     1 次の時間誤差。100 / 200 刻みの Richardson 補外 2 B200 - B100 は空間の h² 誤差だけ
#     (+6e-5) が残るので許容 2e-4
bex=$(awk 'BEGIN { pi = atan2(0, -1); mu = 4e-7 * pi * 1000; tau = 1e6 * mu * 0.2e-3 * 0.2e-3; t = 1e-5
	S = 0; for (k = 1; k < 4000; k += 2) S += 8 / (k * k * pi * pi) * exp(-k * k * pi * pi * t / tau)
	printf "%.12e", mu * 1000 * (1 - S) }')
th_run lamination_step.ofe
b100=$(mt_b)
mt_energy "lamination, 100 steps" 1e-9
awk '/^transient/{print "transient = 1e-5 200"; next} {print}' "$SRC/lamination_step.ofe" > "$WORK/lam200.ofe"
th_run lam200.ofe
b200=$(mt_b)
compare "lamination step : B_mean, Richardson 2 B(200) - B(100) vs Fourier series [T]" \
	"$(awk -v a="$b100" -v b="$b200" 'BEGIN { printf "%.12e", 2 * b - a }')" "$bex" 2e-4
res=$(awk -v a="$b100" -v b="$b200" -v e="$bex" 'BEGIN {
	e1 = (a - e) / e; if (e1 < 0) e1 = -e1
	e2 = (b - e) / e; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.3f)", ((q > 1.85) && (q < 2.25)) ? "OK" : "NG", e1, e2, q }')
echo "  lamination step : first-order in time (100 -> 200 steps) : $res"
case "$res" in NG*) status=1 ;; esac
#     正味の電流 0 の射影を外すと、板の電流が回路の外へ逃げて Fourier 級数から外れる
# (2) σ = 0 の J-A : 刻みを電流掃引の点 (plate_line_ja.ofe) に揃えると解の列が一致する
#     (刻みの式が静磁場そのものになる)。B の列を印字の全桁で比べる
th_run plate_line_ja.ofe
awk '/Jiles-Atherton hysteresis/ { f = 1; next } f && (NF == 6) && ($1 ~ /^[0-9]+$/) && ($2 ~ /e/) { print $1, $2, $4 }' "$WORK/ofe.log" > "$WORK/mt_sweep.txt"
awk '/^currentsweep|^jasub/{next} /^current =/{print; print "transient = 11 11"; print "waveform = 0 0 1 0.5 2 1 3 2 4 1 5 0 6 -1 7 -2 8 -1 9 0 10 1 11 2"; next}
	/^nlsolver/{print; print "jasub = 40"; next} {print}' "$SRC/plate_line_ja.ofe" > "$WORK/mt_ja0.ofe"
th_run mt_ja0.ofe
awk '/H_mean/ { f = 1; next } f && (NF == 7) && ($1 ~ /^[0-9]+$/) && ($3 ~ /e/) { print $1, $3, $6 }' "$WORK/ofe.log" > "$WORK/mt_tr.txt"
if [ -s "$WORK/mt_sweep.txt" ] && cmp -s "$WORK/mt_sweep.txt" "$WORK/mt_tr.txt"; then
	echo "  sigma = 0 transient == quasi-static current sweep (11 steps, all printed digits) -> OK"
else
	echo "  sigma = 0 transient differs from the quasi-static current sweep -> NG" >&2
	diff "$WORK/mt_sweep.txt" "$WORK/mt_tr.txt" | head -4 >&2
	status=1
fi
mt_energy "sigma = 0 hysteresis (input == loop work)" 1e-8
# (3) J-A + 渦電流 (plate_line_ja_eddy.ofe、1 区間 1 ms と 2 ms) : 恒等式、渦電流損が
#     周期に反比例 (比 2)、終わりの B が独立な ODE 積分の値 1.342754 T に 1/T で近づく
th_run plate_line_ja_eddy.ofe
mt_energy "hysteresis + eddy currents, 1 ms per segment" 1e-9
l1=$(mt_loss); bb1=$(mt_b)
awk '/^transient/{print "transient = 0.022 110"; next}
	/^waveform/{ s = "waveform ="; for (i = 3; i <= NF; i += 2) s = s " " (2 * $i) " " $(i+1); print s; next } {print}' \
	"$SRC/plate_line_ja_eddy.ofe" > "$WORK/mt_je2.ofe"
th_run mt_je2.ofe
mt_energy "hysteresis + eddy currents, 2 ms per segment" 1e-9
l2=$(mt_loss); bb2=$(mt_b)
res=$(awk -v a="$l1" -v b="$l2" 'BEGIN { q = (b > 0) ? a / b : 0
	printf "%s (%.4e -> %.4e J/m, ratio %.4f)", ((q > 1.98) && (q < 2.02)) ? "OK" : "NG", a, b, q }')
echo "  classical eddy loss halves when the period doubles : $res"
case "$res" in NG*) status=1 ;; esac
res=$(awk -v a="$bb1" -v b="$bb2" 'BEGIN { r = 1.342754
	e1 = a - r; if (e1 < 0) e1 = -e1
	e2 = b - r; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (|B - B_ODE| %.2e -> %.2e, ratio %.2f)", ((e1 < 5e-5) && (q > 1.7) && (q < 2.3)) ? "OK" : "NG", e1, e2, q }')
echo "  final B approaches the independent J-A ODE value as 1/T : $res"
case "$res" in NG*) status=1 ;; esac
# (4) 弾くべき入力
awk '/^transient/{print "transient = 1e-5 100 0.5"; next} {print}' "$SRC/lamination_step.ofe" > "$WORK/mtx_cn.ofe"
th_reject "transient M with theta = 0.5" mtx_cn.ofe "theta = 1 (backward Euler) only"
awk '/^transient/{print; print "currentsweep = 0.5 1"; next} {print}' "$SRC/lamination_step.ofe" > "$WORK/mtx_sw.ofe"
th_reject "transient M with currentsweep" mtx_sw.ofe "cannot be combined with currentsweep"
awk '/^mur/{print "bh = 2 100 0.5"; print "bh = 2 1000 1.0"; next} {print}' "$SRC/lamination_step.ofe" > "$WORK/mtx_bh.ofe"
th_reject "transient M with a bh curve" mtx_bh.ofe "linear and Jiles-Atherton materials only"

# ---- 3 次四面体 (20 節点、gmsh の tet20 / tri10) ----
#
# 局所節点の並びは gmsh 4.15.2 の出力を重心座標で分類した実測の表 (unstruct.c の
# TET3_A)。形状関数はその多重指数から次数によらない Lagrange の式で作る。
#   (a) 箱 (box_p3.msh) の平行平板 : C が閉形式に全桁一致 (φ が 1 次なので厳密)
#   (b) 2.2 / 4.1 の ASCII とバイナリ : rlc.csv が完全一致 (型 29 / 21 の節点数表)
#   (c) analysis = P : 3 次の φ の恒等式 (機械精度)。**3 次特有の基底を検査できるのは
#       これだけ** (1 次・2 次の場は 2 次の基底でも厳密に表せる)
#   (d) 曲がった同軸 : 3 次が h⁴ で収束すること (体積と C')。2 次に勝つことは
#       assert しない (円の等間隔補間の面積誤差は 2 次も 3 次も h⁴、実測で 2 次の方が
#       小さい。ここを取り違えると正しい実装が落ちる偽陽性になる)
#   (e) 弾くべき入力 (理由の文字列まで)
echo "[order-3 tet] 20-node tetrahedra (gmsh tet20 / tri10)"
for m in "$SRC"/box_p3*.msh "$SRC"/coax_p3*.msh; do cp "$m" "$WORK/"; done
run_case box_p3
compare "C (tet20 box) [F]" "$(value_of C)" 1.77083756e-13 1e-8
bin_same "order-3 tetrahedra, binary 2.2" box_p3.ofe box_p3_bin.msh box_p3_bin_22.msh
bin_same "order-3 tetrahedra, binary 4.1" box_p3.ofe box_p3_bin.msh box_p3_bin_41.msh
sed 's/^mesh = .*/mesh = box_p3_41.msh/' "$SRC/box_p3.ofe" > "$WORK/p3_41.ofe"
(cd "$WORK" && "$OFE" -n 2 p3_41.ofe > /dev/null && "$OFE_POST" > /dev/null)
grep -v '^title' "$WORK/rlc.csv" > "$WORK/p3_41.csv"
if cmp -s "$WORK/bm_a.csv" "$WORK/p3_41.csv"; then echo "  order-3 tetrahedra, ASCII 4.1 : identical -> OK"
else echo "  order-3 tetrahedra, ASCII 4.1 : differ -> NG" >&2; status=1; fi
run_case nodal_test_p3
if grep -q "(f) cubic" "$WORK/ofe.log" && grep -q "result : passed" "$WORK/ofe.log"; then
	echo "  P3 self test (cubic field identity, machine precision) : $(awk '/\(f\) cubic/ { print $NF; exit }' "$WORK/ofe.log") -> OK"
else
	echo "  P3 self test : the cubic-field identity was not run or failed -> NG" >&2
	status=1
fi
# (d) 同軸 : 閉形式と比べ (h = 0.4mm で 0.032%、0.2mm で 0.0017%)、誤差比と体積の比が h⁴
run_case coax_p3
cp3=$(value_of C)
compare "C (tet20 coax, h = 0.4 mm) [F/m]" "$cp3" 1.063417e-10 0.0005
compare "L (tet20 coax, h = 0.4 mm) [H/m]" "$(value_of L)" 2.1972246e-07 0.0005
sed 's/^mesh = .*/mesh = coax_p3f.msh/' "$SRC/coax_p3.ofe" > "$WORK/coax_p3f.ofe"
(cd "$WORK" && "$OFE" -n 2 coax_p3f.ofe > /dev/null && "$OFE_POST" > /dev/null)
cp3f=$(value_of C)
compare "C (tet20 coax, h = 0.2 mm) [F/m]" "$cp3f" 1.063417e-10 0.0001
res=$(awk -v a="$cp3" -v b="$cp3f" 'BEGIN { ce = 2 * atan2(0, -1) * 8.8541878128e-12 * 2.1 / log(3)
	e1 = (a - ce) / ce; if (e1 < 0) e1 = -e1
	e2 = (b - ce) / ce; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.1f, h^4 = 16)", ((q > 12) && (q < 28)) ? "OK" : "NG", e1, e2, q }')
echo "  tet20 coax, C' converges as h^4 (0.4 -> 0.2 mm) : $res"
case "$res" in NG*) status=1 ;; esac
for h in coax_p3 coax_p3f; do
	sed -e "s/^mesh = .*/mesh = $h.msh/" -e 's/^analysis = .*/analysis = P/' -e '/^tline/d' \
		"$SRC/coax_p3.ofe" > "$WORK/${h}_vol.ofe"
	(cd "$WORK" && "$OFE" -n 2 "${h}_vol.ofe" > /dev/null 2>&1)
	eval "v_$h=\$(awk '/integrated volume/ { print \$4; exit }' \"\$WORK/ofe.log\")"
done
res=$(awk -v a="$v_coax_p3" -v b="$v_coax_p3f" 'BEGIN { ve = atan2(0, -1) * (1.5e-3 ^ 2 - 0.5e-3 ^ 2) * 0.1e-3
	e1 = (a - ve) / ve; if (e1 < 0) e1 = -e1
	e2 = (b - ve) / ve; if (e2 < 0) e2 = -e2
	q = (e2 > 0) ? e1 / e2 : 0
	printf "%s (err %.2e -> %.2e, ratio %.1f)", ((e1 < 1e-4) && (q > 13) && (q < 19)) ? "OK" : "NG", e1, e2, q }')
echo "  tet20 coax, isoparametric volume vs pi (b^2 - a^2) lz, h^4 : $res"
case "$res" in NG*) status=1 ;; esac
# (e) 弾くべき入力
# 3 次元の格子で 10 節点の三角形が 3 次でない四面体と混ざる (四面体だけ 4 節点にする)
awk 'BEGIN { e = 0 } /^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 } e && ($2 == 29) { s = $1 " 4"; for (i = 3; i <= 3 + $3 + 4; i++) s = s " " $i; print s; next } { print }' "$WORK/box_p3.msh" > "$WORK/box_p3_tet4.msh"
sed 's/^mesh = .*/mesh = box_p3_tet4.msh/' "$SRC/box_p3.ofe" > "$WORK/p3x_t4.ofe"
th_reject "10-node triangles on 4-node tetrahedra" p3x_t4.ofe "supported only as the faces of 20-node tetrahedra"
#   analysis E (辺要素の自己検証) は [A order 3+] の (h) で通す
awk 'BEGIN { e = 0 } /^\$Elements/ { e = 1; print; next } /^\$EndElements/ { e = 0 } e && ($2 == 21) { s = $1 " 2"; for (i = 3; i <= 3 + $3 + 3; i++) s = s " " $i; print s; next } { print }' "$WORK/box_p3.msh" > "$WORK/box_p3_tri3.msh"
sed 's/^mesh = .*/mesh = box_p3_tri3.msh/' "$SRC/box_p3.ofe" > "$WORK/p3x_tri.ofe"
th_reject "20-node tetrahedra with 3-node electrode triangles" p3x_tri.ofe "triangle is not"
awk '/^material/{print "material = 1.0 1e6"; print "kthermal = 2 10"; next}
	/^analysis/{print "heatsink = 0 20"; print "thermal = 1"; print "analysis = R"; next} {print}' \
	"$SRC/box_p3.ofe" > "$WORK/p3x_th.ofe"
th_run p3x_th.ofe
th_energy "thermal on a 20-node tetrahedral mesh"

# ---- ベクトルヒステリシス (javector = 1、等方性のベクトル J-A、回転磁界) ----
#
# assemble.c の ja_chi_vec / ja_eval_vec (B 駆動の逆。dH = (μ0 (I + χ))⁻¹ dB、dM = χ dH を Heun 法で積分)。
#   (a) 場が 1 方向の板 (plate_line_ja.ofe) ではスカラーの J-A の表と B・L・反復回数が印字の全桁で一致
#       (テンソルの行列は場が 1 方向なら縦方向がスカラーの dH/dB そのもの)。H は I = 0 の点が
#       ±6e-8 A/m の丸めの雑音なので絶対値 1e-6 A/m で比べる
#   (b) 端子の仕事 Σ Ī ΔΛ = 磁気の仕事 ∫ H・ΔB (台形則。収束した各点で Galerkin の残差が 0 なので
#       反復の打ち切りを除いて厳密)。線形のセルの分は最後の状態の ½∫ν|B|² (台形則が telescoping)
#   (c) 回転磁界 (ja_rotating.ofe、2 ポート) : japrobe の Gauss 点の H の列を折れ線でつなぎ、**H 駆動の
#       順方向**のベクトル J-A を awk で独立に積分した B (ja_vector_fwd.awk) が FEM の B に近づく。
#       FEM は B を、awk は H を刻みの中で直線にするので経路の差が刻みの 2 乗で残る (1 回転 24 点で
#       1.8e-2、48 点で 4.4e-3、比 4.2、実測)。許容は 48 点で 1e-2、比 > 3
echo "[javector] vector Jiles-Atherton hysteresis (rotating field, two-port sweep)"
cp "$SRC/plate_line_ja.ofe" "$SRC/ja_rotating.ofe" "$WORK/"
jv_table() { sed -n '/Jiles-Atherton hysteresis/,/hysteresis work/p' "$WORK/ofe.log" | awk '$1 ~ /^[0-9]+$/'; }
jv_work() {	# jv_work <ラベル> : 端子の仕事 = 磁気の仕事
	res=$(awk '/hysteresis work/ { a = $6 + 0; b = $11 + 0 }
		END { if (a == 0) { printf "NG (no work line)"; exit }
		      d = (a - b) / a; if (d < 0) d = -d
		      printf "%s (terminals %.6e, magnetic %.6e, %.1e)", (d < 1e-6) ? "OK" : "NG", a, b, d }' "$WORK/ofe.log")
	echo "  $1 : terminal work = magnetic work (int H.dB) : $res"
	case "$res" in NG*) status=1 ;; esac
}
# (a) (b)
(cd "$WORK" && "$OFE" -n 2 plate_line_ja.ofe > /dev/null 2>&1) || { echo "  plate_line_ja : the run failed -> NG" >&2; status=1; }
jv_table > "$WORK/jv_s.txt"
jv_work "scalar J-A, plate line"
awk '/^analysis/{print "javector = 1"} {print}' "$SRC/plate_line_ja.ofe" > "$WORK/jv_plate.ofe"
(cd "$WORK" && "$OFE" -n 2 jv_plate.ofe > /dev/null 2>&1) || { echo "  jv_plate : the run failed -> NG" >&2; status=1; }
jv_table > "$WORK/jv_v.txt"
jv_work "vector J-A, plate line"
res=$(awk 'NR == FNR { b[FNR] = $4; l[FNR] = $5; it[FNR] = $6; h[FNR] = $3; n = FNR; next }
	{ m++; dh = $3 - h[FNR]; if (dh < 0) dh = -dh; s = h[FNR]; if (s < 0) s = -s
	  if (($4 != b[FNR]) || ($5 != l[FNR]) || ($6 != it[FNR]) || !(dh <= 1e-6 + 1e-6 * s)) bad++ }
	END { if ((n == 0) || (m != n)) printf "NG (%d vs %d rows)", n, m; else printf "%s (%d rows, %d differ)", bad ? "NG" : "OK", n, bad + 0 }' \
	"$WORK/jv_s.txt" "$WORK/jv_v.txt")
echo "  one-directional field : vector J-A == scalar J-A (B, L, iterations; H to 1e-6 A/m) : $res"
case "$res" in NG*) status=1 ;; esac
#     過渡解析 (渦電流 + J-A、plate_line_ja_eddy.ofe) でも 1 方向の場ならスカラーと一致する
#     (反復の経路が少し違うので打ち切りの範囲、実測 4e-9。許容 1e-6)
cp "$SRC/plate_line_ja_eddy.ofe" "$WORK/"
(cd "$WORK" && "$OFE" -n 2 plate_line_ja_eddy.ofe > /dev/null 2>&1) || { echo "  plate_line_ja_eddy : the run failed -> NG" >&2; status=1; }
cp "$WORK/ofe_transient.csv" "$WORK/jv_ts.csv"
awk '/^analysis/{print "javector = 1"} {print}' "$SRC/plate_line_ja_eddy.ofe" > "$WORK/jv_eddy.ofe"
(cd "$WORK" && "$OFE" -n 2 jv_eddy.ofe > /dev/null 2>&1) || { echo "  jv_eddy : the run failed -> NG" >&2; status=1; }
res=$(awk -F, 'NR == FNR { if (FNR > 1) { l[FNR] = $3; b[FNR] = $5 }; n = FNR; next }
	FNR > 1 { m++; for (k = 0; k < 2; k++) { a = (k ? b[FNR] : l[FNR]); c = (k ? $5 : $3); d = a - c; if (d < 0) d = -d; s = (a < 0) ? -a : a
	          r = (s > 0) ? d / s : d; if (!(r <= e)) e = r } }
	END { if ((n < 2) || (m != n - 1)) printf "NG (rows %d vs %d)", n - 1, m; else printf "%s (%d steps, max rel diff %.1e)", (e <= 1e-6) ? "OK" : "NG", m, e }' \
	"$WORK/jv_ts.csv" "$WORK/ofe_transient.csv")
echo "  transient (eddy currents + J-A), one-directional field : vector == scalar (Lambda, B_mean) : $res"
case "$res" in NG*) status=1 ;; esac
# (c) 回転磁界
jv_rot() {	# jv_rot <1 回転の点数> > ofe : ja_rotating.ofe の掃引を書き換える (A = 1000、2 回転)
	awk -v N="$1" 'BEGIN { pi = atan2(0, -1); A = 1000
		for (k = 1; k <= 4; k++) { s1 = s1 " " A * k / 4; s2 = s2 " 0" }
		for (k = 1; k <= 2 * N; k++) { t = 2 * pi * k / N; s1 = s1 " " sprintf("%.15g", A * cos(t)); s2 = s2 " " sprintf("%.15g", A * sin(t)) } }
		/^currentsweep =/ { print "currentsweep =" s1; next } /^currentsweep2 =/ { print "currentsweep2 =" s2; next } { print }' "$SRC/ja_rotating.ofe"
}
for N in 24 48; do
	jv_rot "$N" > "$WORK/jv_rot$N.ofe"
	(cd "$WORK" && "$OFE" -n 2 "jv_rot$N.ofe" > /dev/null 2>&1) || { echo "  jv_rot$N : the run failed -> NG" >&2; status=1; }
	if grep -q "did not converge" "$WORK/ofe.log"; then
		echo "  rotating field ($N points per turn) : the hysteresis iteration did not converge -> NG" >&2; status=1
	fi
	jv_work "rotating field ($N points per turn)"
	if [ "$N" = 24 ]; then
		#   行列のテンソルの誤り (非対角の符号など) は答えを変えずに反復を遅くするだけなので
		#   (実測 476 回 -> 769 回、仕事の恒等式は打ち切りの分 1.1e-6 で際どく落ちるだけ)、
		#   反復の総数そのものを assert する (上限 600)
		res=$(jv_table | awk '{ s += $NF } END { printf "%s (%d iterations, at most 600)", ((s > 0) && (s <= 600)) ? "OK" : "NG", s }')
		echo "  rotating field (24 points per turn) : total fixed-point iterations (tensor tangent) : $res"
		case "$res" in NG*) status=1 ;; esac
	fi
	set -- $(awk -v ms=1.6e6 -v aa=1100 -v al=1.6e-3 -v kk=400 -v cc=0.2 -f "$SRC/ja_vector_fwd.awk" "$WORK/ofe.log")
	eval "jv_e$N=\$1"
	echo "  rotating field ($N points per turn) : |B_FEM - B_H-driven| / max|B| = $1 (max |B| = $2 T)"
done
#     異方性の J-A (jaaxis、成分ごとの Ms / a / α / k / c) : 全軸を ja と同じにすると等方性の掃引と表が一致し、
#     q 軸 (y) だけ別のパラメータにした回転磁界が、同じ式を H 駆動の順方向で独立に積分した awk
#     (ja_vector_fwd_aniso.awk) に刻みの 2 乗で近づく (実測 2.1e-2 -> 5.9e-3、比 3.6)。等方性の awk とは 0.47 ずれる
#     (パラメータが効いている)。α の添字の取り違えは 48 点で 8.3e-3、k の射影・q 軸の χ_an・c の添字の誤りは
#     24 点で 6e-2 以上になる (どれも実測)
awk '/^ja =/{print; print "jaaxis = 2 x 1.6e6 1100 1.6e-3 400 0.2"; print "jaaxis = 2 y 1.6e6 1100 1.6e-3 400 0.2"; next} {print}' "$WORK/jv_rot24.ofe" > "$WORK/jva_eq.ofe"
(cd "$WORK" && "$OFE" -n 2 jv_rot24.ofe > /dev/null 2>&1 && grep -E "^  probe" ofe.log > jva_iso.txt) || true
(cd "$WORK" && "$OFE" -n 2 jva_eq.ofe > /dev/null 2>&1 && grep -E "^  probe" ofe.log > jva_eq.txt) || true
if [ -s "$WORK/jva_iso.txt" ] && cmp -s "$WORK/jva_iso.txt" "$WORK/jva_eq.txt"; then
	echo "  anisotropic J-A with every axis equal to ja == isotropic : the probe table is identical -> OK"
else
	echo "  anisotropic J-A with every axis equal to ja == isotropic : the probe tables differ -> NG" >&2; status=1
fi
for N in 24 48; do
	awk '/^ja =/{print; print "jaaxis = 2 y 1.2e6 800 1.0e-3 250 0.3"; next} {print}' "$WORK/jv_rot$N.ofe" > "$WORK/jva_rot$N.ofe"
	(cd "$WORK" && "$OFE" -n 2 "jva_rot$N.ofe" > /dev/null 2>&1) || { echo "  jva_rot$N : the run failed -> NG" >&2; status=1; }
	set -- $(awk -v msp=1.6e6 -v aap=1100 -v alp=1.6e-3 -v kkp=400 -v ccp=0.2 -v msq=1.2e6 -v aaq=800 -v alq=1.0e-3 -v kkq=250 -v ccq=0.3 \
		-f "$SRC/ja_vector_fwd_aniso.awk" "$WORK/ofe.log")
	eval "jva_e$N=\$1"
	[ "$N" = 48 ] && jva_iso=$(awk -v ms=1.6e6 -v aa=1100 -v al=1.6e-3 -v kk=400 -v cc=0.2 -f "$SRC/ja_vector_fwd.awk" "$WORK/ofe.log" | awk '{ print $1 }')
done
res=$(awk -v a="$jva_e24" -v b="$jva_e48" -v c="$jva_iso" 'BEGIN { ok = (b <= 7e-3) && (a >= 3 * b) && (c > 0.1)
	printf "%s (24 points %.2e, 48 points %.2e, ratio %.2f; isotropic awk %.2f)", ok ? "OK" : "NG", a, b, (b > 0) ? a / b : 0, c }')
echo "  anisotropic J-A, rotating field : FEM -> independent H-driven anisotropic vector J-A (awk) : $res"
case "$res" in NG*) status=1 ;; esac
awk '/^javector/{next} /^ja =/{print; print "jaaxis = 2 y 1.2e6 800 1.0e-3 250 0.3"; next} {print}' "$SRC/ja_rotating.ofe" > "$WORK/jvax_s.ofe"
th_reject "jaaxis with the scalar J-A" jvax_s.ofe "needs javector = 1"
awk '/^ja =/{print "jaaxis = 2 y 1.2e6 800 1.0e-3 250 0.3"; next} {print}' "$SRC/ja_rotating.ofe" > "$WORK/jvax_n.ofe"
th_reject "jaaxis without a ja line" jvax_n.ofe "has no ja line"
#     sample の 24 点は jv_rot 24 と同じ掃引 (書き換えで同じ入力になること)
awk '/^currentsweep/' "$SRC/ja_rotating.ofe" > "$WORK/jv_a.txt"
awk '/^currentsweep/' "$WORK/jv_rot24.ofe" > "$WORK/jv_b.txt"
cmp -s "$WORK/jv_a.txt" "$WORK/jv_b.txt" || { echo "  ja_rotating.ofe is not the 24-point sweep -> NG" >&2; status=1; }
res=$(awk -v a="$jv_e24" -v b="$jv_e48" 'BEGIN { ok = (b <= 1e-2) && (a > 3 * b)
	printf "%s (24 points %.2e, 48 points %.2e, ratio %.2f)", ok ? "OK" : "NG", a, b, (b > 0) ? a / b : 0 }')
echo "  rotating field : FEM (B-driven) -> independent H-driven vector J-A as the turn is refined : $res"
case "$res" in NG*) status=1 ;; esac
# (d) 弾くべき入力 (理由まで grep する)
awk '!/^conductor = 2/' "$SRC/ja_rotating.ofe" > "$WORK/jvx_1p.ofe"
th_reject "currentsweep2 with one port" jvx_1p.ofe "needs exactly two ports"
awk '/^currentsweep2/{print "currentsweep2 = 0 0"; next} {print}' "$SRC/ja_rotating.ofe" > "$WORK/jvx_n.ofe"
th_reject "currentsweep2 with a different number of points" jvx_n.ofe "as many points as currentsweep"
awk '/^ja =/{next} /^currentsweep =/{next} /^javector/{next} /^japrobe/{next} /^jasub/{next} {print}' "$SRC/ja_rotating.ofe" > "$WORK/jvx_noja.ofe"
th_reject "currentsweep2 without a hysteretic material" jvx_noja.ofe "currentsweep2 needs a hysteretic"

# ---- 3 次元渦電流の 2 次の要素 (analysis = A、10 節点四面体、2 次の Nedelec + 2 次の節点、eddy3d_p2.c) ----
#
# bar_eddy と同じ棒を gmsh が 2 次にした格子 (bar_tet2.msh、gmsh_bar_p2.py) で解き、1 次元の厳密解
# (bar_eddy_z.awk) と比べる。
#   (a) 100 Hz : DC の場 (B が厚さ方向に 1 次) は 2 次の空間に厳密に入るので R・L が 1e-7 以内 (実測 1e-9)
#   (b) 100 kHz (t/δ = 4.8) : 2 次は 3e-5 以内 (実測 8e-6)、同じ頂点の 1 次 (1e-2) より 100 倍以上よい
#   (c) 剛体回転した格子 (J が対角でない) で R・L が 9 桁一致
#   (d) 内部の中間節点をずらした曲がった要素 : R の DC 極限 (1 次の φ は等パラメトリックな 2 次の空間に入る)
#       が 1e-7、100 kHz が 3e-5 以内。ずらしすぎて裏返った要素は setup で弾く (飛ばすと R +8.7%、実測)
echo "[A order 2] 3-D eddy current on 10-node tetrahedra (second-order Nedelec / nodal elements)"
cp "$SRC/bar_eddy.ofe" "$SRC/bar_tet.msh" "$SRC/bar_tet2.msh" "$WORK/"
p2_rl() {	# p2_rl <ofe> : 実行して "R L" (rlc.csv の Rf / Lf)
	(cd "$WORK" && "$OFE" -n 2 "$1" > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "x x"; return; }
	awk -F, '$1 == "Rf" { a = 1; next } a == 1 && $1 == "1" { r = $2; a = 0 } $1 == "Lf" { b = 1; next } b == 1 && $1 == "1" { l = $2; b = 0 }
		END { printf "%s %s", r, l }' "$WORK/rlc.csv"
}
p2_case() {	# p2_case <ラベル> <格子> <周波数> <許容> : 厳密解との R・L の相対誤差
	grep -v '^#' "$WORK/bar_eddy.ofe" | awk -v m="$2" -v f="$3" '/^mesh/{print "mesh = " m; next} /^frequency/{print "frequency = " f; next} {print}' > "$WORK/p2c.ofe"
	set -- "$1" "$2" "$3" "$4" $(p2_rl p2c.ofe) $(awk -v f="$3" -f "$SRC/bar_eddy_z.awk")
	res=$(awk -v r="$5" -v l="$6" -v re="$7" -v le="$8" -v tol="$4" 'BEGIN { if (r == "x") { print "NG (the run failed)"; exit }
		er = r / re - 1; el = l / le - 1; ok = "NG"; if ((er < tol) && (er > -tol) && (el < tol) && (el > -tol)) ok = "OK"
		printf "%s (R %+.2e, L %+.2e)", ok, er, el }')
	echo "  $1 : R, L == 1-D exact solution within $4 : $res"
	case "$res" in NG*) status=1 ;; esac
	P2R=$5; P2L=$6
}
p2_case "order 2, 100 Hz" bar_tet2.msh 1e2 1e-7
p2_case "order 2, 10 kHz" bar_tet2.msh 1e4 1e-6
p2_case "order 2, 100 kHz" bar_tet2.msh 1e5 3e-5
r2=$P2R; l2=$P2L
p2_case "order 1 (same vertices), 100 kHz (reference only)" bar_tet.msh 1e5 3e-2
set -- $(awk -v f=1e5 -f "$SRC/bar_eddy_z.awk")
res=$(awk -v r1="$P2R" -v l1="$P2L" -v r2="$r2" -v l2="$l2" -v re="$1" -v le="$2" 'BEGIN {
	e1 = r1 / re - 1; e2 = r2 / re - 1; f1 = l1 / le - 1; f2 = l2 / le - 1
	if (e1 < 0) e1 = -e1; if (e2 < 0) e2 = -e2; if (f1 < 0) f1 = -f1; if (f2 < 0) f2 = -f2
	ok = "NG"; if ((e2 * 100 < e1) && (f2 * 100 < f1)) ok = "OK"
	printf "%s (R %.0fx, L %.0fx)", ok, e1 / e2, f1 / f2 }')
echo "  order 2 beats order 1 on the same vertices at 100 kHz by more than 100x : $res"
case "$res" in NG*) status=1 ;; esac
# (c)
awk 'BEGIN { a = 0.5; b = 0.7; ca = cos(a); sa = sin(a); cb = cos(b); sb = sin(b) }
	/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; x1 = ca * x - sa * y; y1 = sa * x + ca * y; y2 = cb * y1 - sb * z; z2 = sb * y1 + cb * z; printf "%s %.17g %.17g %.17g\n", $1, x1, y2, z2; next }
	{ print }' "$WORK/bar_tet2.msh" > "$WORK/bar_tet2_rot.msh"
p2_case "order 2, rigidly rotated mesh, 100 kHz" bar_tet2_rot.msh 1e5 3e-5
res=$(awk -v a="$P2R" -v b="$r2" -v c="$P2L" -v d="$l2" 'BEGIN { x = a / b - 1; y = c / d - 1; if (x < 0) x = -x; if (y < 0) y = -y
	printf "%s (%.1e, %.1e)", ((x < 1e-8) && (y < 1e-8)) ? "OK" : "NG", x, y }')
echo "  rigidly rotated mesh : R, L == axis-aligned mesh : $res"
case "$res" in NG*) status=1 ;; esac
# (d) 内部の中間節点 (要素の 10 節点のうち 5..10 番目) だけを決まった擬似乱数でずらす
awk '/\$Elements/{f=1;getline;next} /\$EndElements/{f=0} f && $2==11 { nt = $3; for (k = 4 + nt + 4; k <= 4 + nt + 9; k++) print $k }' "$WORK/bar_tet2.msh" | sort -u > "$WORK/p2mid.txt"
p2_bend() {	# p2_bend <振幅 [m]> <出力>
	awk -v amp="$1" -v mf="$WORK/p2mid.txt" 'BEGIN { s = 12345; while ((getline l < mf) > 0) mid[l] = 1 }
		function rnd() { s = (s * 1103515245 + 12345) % 2147483648; return s / 2147483648 }
		/^\$Nodes/ { f = 1; print; getline; print; next } /^\$EndNodes/ { f = 0 }
		f && NF == 4 { x = $2; y = $3; z = $4; e = 1e-12
			inner = (x > e) && (x < 2e-3 - e) && (y > e) && (y < 0.25e-3 - e) && (z > e) && (z < 1e-3 - e)
			if (inner && ($1 in mid)) { x += (rnd() - 0.5) * amp; y += (rnd() - 0.5) * amp; z += (rnd() - 0.5) * amp }
			printf "%s %.17g %.17g %.17g\n", $1, x, y, z; next } { print }' "$WORK/bar_tet2.msh" > "$WORK/$2"
}
p2_bend 2e-6 bar_tet2_cm.msh
p2_case "order 2, curved elements (interior mid-edge nodes moved), 100 Hz R / L" bar_tet2_cm.msh 1e2 1e-5
grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_tet2_cm.msh"; next} /^frequency/{next} /^awall/{next} /^analysis/{print "analysis = R"; next} {print}' > "$WORK/p2r.ofe"
(cd "$WORK" && "$OFE" -n 2 p2r.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  p2r : the run failed -> NG" >&2; status=1; }
compare "curved elements : analysis R == L / (sigma W t) (a linear phi is in the isoparametric space) [ohm]" \
	"$(awk -F, '$1 == "R" { a = 1; next } a && $1 == "1" { print $2; exit }' "$WORK/rlc.csv")" 1.3793103448e-04 1e-8
p2_case "order 2, curved elements, 100 kHz" bar_tet2_cm.msh 1e5 3e-5
p2_bend 2e-5 bar_tet2_fold.msh
sed 's/^mesh = .*/mesh = bar_tet2_fold.msh/' "$WORK/p2r.ofe" > "$WORK/p2x_fold.ofe"
th_reject "10-node mesh with folded elements" p2x_fold.ofe "is degenerate or inverted"
# (e) 対応しない組み合わせ
grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_tet2.msh"; next} {print}' > "$WORK/p2base.ofe"
awk '/^frequency/{next} /^analysis/{print "gauge = 1"; print "transient = 2e-4 5"} {print}' "$WORK/p2base.ofe" > "$WORK/p2x_g.ofe"
th_reject "order-2 A with gauge = 1 in the time domain" p2x_g.ofe "second-order (10-node) mesh supports gauge = 0"
# (f) 場の出力 (重心の B・J、ゲージ不変)。100 Hz では J が一様なので 1 次と同じ恒等式が使える :
#     体積・|J| = σV/L・符号つきの平均・場の損失 Σ|J|²V/(2σ) = ½Re(V I)
sed -e "s/^frequency = .*/frequency = 1e2/" -e "s/^analysis = A/fieldout = 1\nanalysis = A/" "$WORK/p2base.ofe" > "$WORK/p2_fld.ofe"
(cd "$WORK" && "$OFE" -n 2 p2_fld.ofe > /dev/null 2>&1) || { echo "  p2_fld : the run failed -> NG" >&2; status=1; }
compare "order 2, field output : volume of the 10-node cells [m^3]" "$(vtk vol J_A_re)" 5.0e-10 1e-9
compare "order 2, field output : |J| max (sigma V / L) [A/m^2]" "$(vtk vmax J_A_re)" 2.9e10 1e-4
compare "order 2, field output : |J| min [A/m^2]" "$(vtk vmin J_A_re)" 2.9e10 1e-4
compare "order 2, field output : mean J_x (signed) [A/m^2]" "$(vtk mean0 J_A_re)" -2.9e10 1e-4
pf=$(awk -v r="$(vtk int2 J_A_re)" -v m="$(vtk int2 J_A_im)" 'BEGIN{ printf "%.10e", (r + m) / (2 * 5.8e7) }')
pt=$(awk '/port 1 : I =/{ printf "%.10e", 0.5 * $6 }' "$WORK/ofe.log")
compare "order 2, field output : ohmic loss from the field == 1/2 Re(V I) [W]" "$pf" "$pt" 1e-6
# (g) 辺要素の自己検証の 2 次版 (analysis = E、box_p2.msh) : 2 次の Nedelec 空間に入る場の L2 射影の誤差 0 と
#     回転のエネルギー。(P2a) 一般の 1 次の場、(P2b) それに S2 の 2 次の場 (面の関数が要る)。等方・異方性の ν
#     (anisomur) と剛体回転した格子で
cp "$SRC/edge_test.ofe" "$SRC/edge_test_aniso.ofe" "$SRC/box_p2.msh" "$WORK/"
awk 'BEGIN { a = 0.5; b = 0.7; ca = cos(a); sa = sin(a); cb = cos(b); sb = sin(b) }
	/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; x1 = ca * x - sa * y; y1 = sa * x + ca * y; y2 = cb * y1 - sb * z; z2 = sb * y1 + cb * z; printf "%s %.17g %.17g %.17g\n", $1, x1, y2, z2; next }
	{ print }' "$WORK/box_p2.msh" > "$WORK/box_p2_rot.msh"
for c in "edge_test box_p2" "edge_test_aniso box_p2" "edge_test box_p2_rot"; do
	set -- $c
	sed "s/^mesh = .*/mesh = $2.msh/" "$WORK/$1.ofe" > "$WORK/p2e_$1_$2.ofe"
	(cd "$WORK" && "$OFE" -n 2 "p2e_$1_$2.ofe" > /dev/null 2>&1) || true
	np=$(grep -c "^  (P2[ab]) .* passed$" "$WORK/ofe.log" || true)
	if [ "$np" = 2 ] && grep -q "second-order edge self test : passed" "$WORK/ofe.log"; then
		echo "  second-order edge self test ($1 on $2) : (P2a) linear field, (P2b) + S2 quadratic : passed -> OK"
	else
		echo "  second-order edge self test ($1 on $2) : $np of 2 identities passed -> NG" >&2; status=1
		grep "(P2" "$WORK/ofe.log" | cut -c1-160 >&2
	fi
done
#     恒等式が厳密なのは直線の要素だけなので、曲がった格子 (coax_p2) は理由をつけて弾く (誤りと区別できるように)
cp "$SRC/coax_p2.msh" "$WORK/"
sed "s/^mesh = .*/mesh = coax_p2.msh/" "$WORK/edge_test.ofe" > "$WORK/p2x_ecurv.ofe"
mesh_reject "second-order edge self test on a curved mesh" p2x_ecurv.ofe
grep -q "needs straight-sided elements" "$WORK/ofe.log" || { echo "  *** it was rejected for a different reason than: needs straight-sided elements" >&2; status=1; }

# (h) B-H 曲線 (周波数領域の実効磁気抵抗率法、ν は**積分点ごと**)。bar_air_bh を粗い格子 (bar_air_c.msh、
#     mkmesh.py bar_air -nx 8 -ny 1 -nz 8 -nza 2) を gmsh が 2 次にした bar_air_p2.msh (gmsh_bar_p2.py) で解く:
#     ・層が非線形 (H が一様) : 自己無撞着な閉形式 (bar_layer.awk) に 3 つの電流で 1e-7 (実測 4e-9。
#       層の A は 1 次なので 2 次の空間に入る。同じ頂点の 1 次は L −0.25%)
#     ・導体が非線形 (10 Hz、0.5 A。B が厚さ方向に変わり曲線の 3 区間を通る) : 独立な 1 次元の参照解
#       (bar_cond_bh.awk、n = 4000 で h² の 1e-8) に R 2e-4・L 1e-3 (実測 −6.5e-5 / −3.6e-4。曲線の折れ目で ν が
#       滑らかでないので h の 1 次でしか縮まない)、同じ頂点の 1 次 (+5.2e-3 / −7.1e-3) より 10 倍以上よい。
#       **ν を要素の 1 点で代表させる誤り (R +1.8e-3) と B の虚部を落とす誤り (R +5.0e-4) はこちらでだけ落ちる**
#       (層の B は一様なので層の検査は素通りする)
cp "$SRC/bar_air_bh.ofe" "$SRC/bar_air_c.msh" "$SRC/bar_air_p2.msh" "$WORK/"
for I in 150 200 250; do
	grep -v '^#' "$WORK/bar_air_bh.ofe" | sed -e 's/^mesh = .*/mesh = bar_air_p2.msh/' -e "s/^current = .*/current = $I/" > "$WORK/p2bh.ofe"
	set -- $(p2_rl p2bh.ofe) $(awk -v f=1e2 -v mode=nl -v I="$I" -f "$SRC/bar_layer.awk")
	compare "order 2, nonlinear layer, I = $I A : R == self-consistent closed form [ohm]" "$1" "$3" 1e-7
	compare "order 2, nonlinear layer, I = $I A : L == self-consistent closed form [H]" "$2" "$4" 1e-7
done
p2bhc() {	# p2bhc <格子> : 導体を非線形にした 10 Hz、0.5 A の "R L"
	grep -v '^#' "$WORK/bar_air_bh.ofe" | awk -v m="$1" '/^mesh/{print "mesh = " m; next} /^bh = 3/{sub(/^bh = 3/, "bh = 2")}
		/^frequency/{print "frequency = 10"; next} /^current/{print "current = 0.5"; next} {print}' > "$WORK/p2bhc.ofe"
	p2_rl p2bhc.ofe
}
set -- $(p2bhc bar_air_p2.msh) $(p2bhc bar_air_c.msh) $(awk -v f=10 -v I=0.5 -v n=4000 -f "$SRC/bar_cond_bh.awk")
res=$(awk -v r2="$1" -v l2="$2" -v r1="$3" -v l1="$4" -v re="$5" -v le="$6" 'function ab(x) { return (x < 0) ? -x : x }
	BEGIN { if ((r2 == "x") || (r1 == "x")) { print "NG (a run failed)"; exit }
	e2 = r2 / re - 1; f2 = l2 / le - 1; e1 = r1 / re - 1; f1 = l1 / le - 1
	ok = (ab(e2) <= 2e-4) && (ab(f2) <= 1e-3) && (10 * ab(e2) < ab(e1)) && (10 * ab(f2) < ab(f1))
	printf "%s (order 2 : R %+.2e, L %+.2e; order 1 : R %+.2e, L %+.2e)", ok ? "OK" : "NG", e2, f2, e1, f1 }')
echo "  order 2, nonlinear conductor (10 Hz, 0.5 A) : R, L == 1-D reference (2e-4 / 1e-3) and 10x better than order 1 : $res"
case "$res" in NG*) status=1 ;; esac
awk '/^analysis/{print "bh = 2 100 0.5"} {print}' "$WORK/edge_test.ofe" | sed 's/^mesh = .*/mesh = box_p2.msh/' > "$WORK/p2x_ebh.ofe"
th_reject "order-2 edge self test (E) with a bh material" p2x_ebh.ofe "supports linear materials"
# (i) ゲージ固定 (gauge = 1) と直接解法 (direct = 1)。Whitney は awall 優先の全域木、∇(λaλb) の自由度は全部 0。
#     Z はゲージ不変なので gauge = 0 と一致する (反復の打ち切りの範囲)。**直接解法 (ピボットなし LDLᵀ) が通ることが
#     零空間が消えた検査**: 勾配の自由度か木のどちらかを固定し忘れると特異なピボットで落ちる (どちらも実測)。
#     bar_air (空気の層つき) を粗い 2 次の格子で、100 kHz
grep -v '^#' "$SRC/bar_air.ofe" | awk '/^mesh/{print "mesh = bar_air_p2.msh"; next} /^frequency/{print "frequency = 1e5"; next} {print}' > "$WORK/p2g0.ofe"
awk '/^analysis/{print "gauge = 1"} {print}' "$WORK/p2g0.ofe" > "$WORK/p2g1.ofe"
awk '/^analysis/{print "gauge = 1"; print "direct = 1"} {print}' "$WORK/p2g0.ofe" > "$WORK/p2g2.ofe"
set -- $(p2_rl p2g0.ofe) $(p2_rl p2g1.ofe)
compare "order 2, gauge = 1 (tree + edge gradients) == gauge = 0 : R [ohm]" "$3" "$1" 1e-7
compare "order 2, gauge = 1 (tree + edge gradients) == gauge = 0 : L [H]" "$4" "$2" 1e-7
r0=$1; l0=$2
set -- $(p2_rl p2g2.ofe)
compare "order 2, gauge = 1 + direct = 1 (LDL^T, no pivoting) == gauge = 0 : R [ohm]" "$1" "$r0" 1e-6
compare "order 2, gauge = 1 + direct = 1 (LDL^T, no pivoting) == gauge = 0 : L [H]" "$2" "$l0" 1e-6
awk '/^awall/{next} /^analysis/{print "gauge = 1"} {print}' "$WORK/p2g0.ofe" > "$WORK/p2x_gw.ofe"
#     電極の面を awall にしない格子 (bar_air_yw_p2.msh、mkmesh.py bar_air -nx 4 -ny 1 -nz 4 -nza 2 -yawall 1 を
#     gmsh が 2 次に) : 電極の辺の ∇(λaλb) は ψ の泡が 0 に縛られていてゲージの自由度でないので固定しない。
#     固定すると Z が黙って変わる (この格子で R −2.7%、実測)。bar_air_p2 は電極の面も awall なので該当する辺が無い
cp "$SRC/bar_air_yw_p2.msh" "$WORK/"
for g in 0 1 2; do
	awk -v g="$g" '/^mesh/{print "mesh = bar_air_yw_p2.msh"; next} /^analysis/{ if (g > 0) print "gauge = 1"; if (g > 1) print "direct = 1" } {print}' \
		"$WORK/p2g0.ofe" > "$WORK/p2y$g.ofe"
done
set -- $(p2_rl p2y0.ofe) $(p2_rl p2y1.ofe) $(p2_rl p2y2.ofe)
compare "order 2, electrodes off awall : gauge = 1 == gauge = 0 : R [ohm]" "$3" "$1" 1e-7
compare "order 2, electrodes off awall : gauge = 1 == gauge = 0 : L [H]" "$4" "$2" 1e-7
compare "order 2, electrodes off awall : gauge = 1 + direct = 1 == gauge = 0 : R [ohm]" "$5" "$1" 1e-6
compare "order 2, electrodes off awall : gauge = 1 + direct = 1 == gauge = 0 : L [H]" "$6" "$2" 1e-6
mesh_reject "order-2 A, gauge = 1 with electrodes not connected through awall" p2x_gw.ofe
grep -q "not connected through awall" "$WORK/ofe.log" || { echo "  *** it was rejected for a different reason than: not connected through awall" >&2; status=1; }
# ---- 3 次以上の Nedelec 要素 (analysis = A / E、3〜7 次の四面体、eddy3d_pr.c) ----
#
# A を次数 r の Nedelec (第 1 種、Arnold-Falk-Winther の基底 λ^α φ_σ)、φ を次数 r の節点要素で解く。格子は
# bar_eddy と同じ棒を gmsh が幾何から次数 r で切ったもの (gmsh_bar_ho.py、h = 0.25 mm。離散の格子は gmsh が
# 3 次以上に上げられない)。bar_ho2.msh は同じ頂点の 2 次 (eddy3d_p2.c で解く対照)。
#   (a) 100 Hz : DC の場は空間に入るので R・L が 1e-7 (実測 1.5e-9)
#   (b) 100 kHz (t/δ = 4.8) : 3 次は 1e-4 以内 (実測 −1.9e-5) で同じ頂点の 2 次 (+2.7e-3) より 30 倍以上よい。
#       4 次は 2e-6 以内 (実測 2e-7) で 3 次より 30 倍以上よい (誤差は次数で指数的に減る。5 次は 3e-10 で CI では回さない)
#   (c) 剛体回転した格子 (J が対角でない) で R・L が 8 桁一致
#   (d) 内部の節点 (頂点以外) をずらした曲がった要素 : 100 Hz が 1e-5、100 kHz が 1e-4、analysis R が閉形式に 1e-8
#   (e) 空気層 (導体と空気の 2 材料、bar_air の閉形式) : 3 次が 1e-4 以内 (実測 −2.1e-5) で同じ頂点の 2 次より 30 倍よい
#   (f) 磁性導体 (mur = 50、100 Hz、t/δ = 1.1) : 閉形式に 1e-5
#   (g) 場の出力 (100 Hz、重心の B・J) : 2 次と同じ恒等式 (体積・|J| = σV/L・符号つき平均・損失 = ½Re(V I))
#   (h) 辺要素の自己検証 (analysis = E) : 3〜7 次の直線の箱 (box_p3..p7) で、次数 r − 1 の一般の場 (Pra) と
#       r' × q (q は r − 1 次の同次式、P_{r−1}³ に入らない、Prb) の σ 重み L2 射影の誤差 0・回転のエネルギー。
#       3 次は異方性の ν と剛体回転した格子でも。曲がった格子 (coax_p3) は理由をつけて弾く
#   (i) 対応しない組み合わせ (ゲージ・過渡解析・B-H 曲線) は理由をつけて弾く (交流の発熱は [thermal] の (2c)。直接解法はゲージなしの A を
#       入力段が弾く)
echo "[A order 3+] 3-D eddy current on tetrahedra of order 3 and up (order-r Nedelec / nodal elements)"
cp "$SRC/bar_ho2.msh" "$SRC/bar_ho3.msh" "$SRC/bar_ho4.msh" "$SRC/bar_air_ho2.msh" "$SRC/bar_air_ho3.msh" "$WORK/"
pr_ratio() {	# pr_ratio <ラベル> <R の誤差 (低)> <L の誤差 (低)> <R の誤差 (高)> <L の誤差 (高)> <倍率>
	res=$(awk -v a="$2" -v b="$3" -v c="$4" -v d="$5" -v k="$6" 'BEGIN { if (a < 0) a = -a; if (b < 0) b = -b; if (c < 0) c = -c; if (d < 0) d = -d
		ok = "NG"; if ((c * k < a) && (d * k < b)) ok = "OK"; printf "%s (R %.0fx, L %.0fx)", ok, a / c, b / d }')
	echo "  $1 by more than $6x : $res"
	case "$res" in NG*) status=1 ;; esac
}
pr_err() {	# pr_err <R> <L> <R 厳密> <L 厳密> : "R の相対誤差 L の相対誤差"
	awk -v r="$1" -v l="$2" -v re="$3" -v le="$4" 'BEGIN { if (r == "x") { print "1 1"; exit } printf "%.6e %.6e", r / re - 1, l / le - 1 }'
}
p2_case "order 3, 100 Hz" bar_ho3.msh 1e2 1e-7
p2_case "order 2 (same vertices), 100 kHz (reference only)" bar_ho2.msh 1e5 1e-2
set -- $(awk -v f=1e5 -f "$SRC/bar_eddy_z.awk"); bre=$1; ble=$2
set -- $(pr_err "$P2R" "$P2L" "$bre" "$ble"); e2r=$1; e2l=$2
p2_case "order 3, 100 kHz" bar_ho3.msh 1e5 1e-4
r3=$P2R; l3=$P2L
set -- $(pr_err "$P2R" "$P2L" "$bre" "$ble"); e3r=$1; e3l=$2
pr_ratio "order 3 beats order 2 on the same vertices at 100 kHz" "$e2r" "$e2l" "$e3r" "$e3l" 30
p2_case "order 4, 100 kHz" bar_ho4.msh 1e5 2e-6
set -- $(pr_err "$P2R" "$P2L" "$bre" "$ble"); e4r=$1; e4l=$2
pr_ratio "order 4 beats order 3 on the same vertices at 100 kHz" "$e3r" "$e3l" "$e4r" "$e4l" 30
# (c)
awk 'BEGIN { a = 0.5; b = 0.7; ca = cos(a); sa = sin(a); cb = cos(b); sb = sin(b) }
	/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; x1 = ca * x - sa * y; y1 = sa * x + ca * y; y2 = cb * y1 - sb * z; z2 = sb * y1 + cb * z; printf "%s %.17g %.17g %.17g\n", $1, x1, y2, z2; next }
	{ print }' "$WORK/bar_ho3.msh" > "$WORK/bar_ho3_rot.msh"
p2_case "order 3, rigidly rotated mesh, 100 kHz" bar_ho3_rot.msh 1e5 1e-4
res=$(awk -v a="$P2R" -v b="$r3" -v c="$P2L" -v d="$l3" 'BEGIN { x = a / b - 1; y = c / d - 1; if (x < 0) x = -x; if (y < 0) y = -y
	printf "%s (%.1e, %.1e)", ((x < 1e-8) && (y < 1e-8)) ? "OK" : "NG", x, y }')
echo "  order 3, rigidly rotated mesh : R, L == axis-aligned mesh : $res"
case "$res" in NG*) status=1 ;; esac
# (d) 内部の頂点以外の節点 (20 節点の 5..20 番目) を決まった擬似乱数でずらす
awk '/\$Elements/{f=1;getline;next} /\$EndElements/{f=0} f && $2==29 { nt = $3; for (k = 4 + nt + 4; k <= 4 + nt + 19; k++) print $k }' "$WORK/bar_ho3.msh" | sort -u > "$WORK/p3mid.txt"
awk -v amp=2e-6 -v mf="$WORK/p3mid.txt" 'BEGIN { s = 12345; while ((getline l < mf) > 0) mid[l] = 1 }
	function rnd() { s = (s * 1103515245 + 12345) % 2147483648; return s / 2147483648 }
	/^\$Nodes/ { f = 1; print; getline; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; e = 1e-12
		inner = (x > e) && (x < 2e-3 - e) && (y > e) && (y < 0.25e-3 - e) && (z > e) && (z < 1e-3 - e)
		if (inner && ($1 in mid)) { x += (rnd() - 0.5) * amp; y += (rnd() - 0.5) * amp; z += (rnd() - 0.5) * amp; nm++ }
		printf "%s %.17g %.17g %.17g\n", $1, x, y, z; next } { print }
	END { if (nm < 100) print "*** only " nm " nodes moved" > "/dev/stderr" }' "$WORK/bar_ho3.msh" > "$WORK/bar_ho3_cm.msh"
p2_case "order 3, curved elements (interior nodes moved), 100 Hz" bar_ho3_cm.msh 1e2 1e-5
grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_ho3_cm.msh"; next} /^frequency/{next} /^awall/{next} /^analysis/{print "analysis = R"; next} {print}' > "$WORK/p3r.ofe"
(cd "$WORK" && "$OFE" -n 2 p3r.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  p3r : the run failed -> NG" >&2; status=1; }
compare "order 3, curved elements : analysis R == L / (sigma W t) [ohm]" \
	"$(awk -F, '$1 == "R" { a = 1; next } a && $1 == "1" { print $2; exit }' "$WORK/rlc.csv")" 1.3793103448e-04 1e-8
p2_case "order 3, curved elements, 100 kHz" bar_ho3_cm.msh 1e5 1e-4
# (e) 空気層 (bar_air.ofe の閉形式、100 kHz : R = 4.7930357804e-04, L = 9.7419014024e-10)
for o in 2 3; do
	grep -v '^#' "$SRC/bar_air.ofe" | awk -v m="bar_air_ho$o.msh" '/^mesh/{print "mesh = " m; next} /^frequency/{print "frequency = 1e5"; next} {print}' > "$WORK/pra$o.ofe"
	set -- $(p2_rl "pra$o.ofe")
	set -- $(pr_err "$1" "$2" 4.7930357804e-04 9.7419014024e-10)
	eval "ea${o}r=$1; ea${o}l=$2"
done
res=$(awk -v r="$ea3r" -v l="$ea3l" 'BEGIN { ok = "NG"; if ((r < 1e-4) && (r > -1e-4) && (l < 1e-4) && (l > -1e-4)) ok = "OK"; printf "%s (R %+.2e, L %+.2e)", ok, r, l }')
echo "  order 3, conductor + air (2 materials), 100 kHz : R, L == bar_air closed form within 1e-4 : $res"
case "$res" in NG*) status=1 ;; esac
pr_ratio "order 3 beats order 2 with the air layer at 100 kHz" "$ea2r" "$ea2l" "$ea3r" "$ea3l" 30
# (f)
grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_ho3.msh"; next} /^frequency/{print "frequency = 1e2"; next} /^region = 1 2/{print; print "mur = 2 50"; next} {print}' > "$WORK/prmu.ofe"
set -- $(p2_rl prmu.ofe) $(awk -v f=1e2 -v mur=50 -f "$SRC/bar_eddy_z.awk")
set -- $(pr_err "$1" "$2" "$3" "$4")
res=$(awk -v r="$1" -v l="$2" 'BEGIN { ok = "NG"; if ((r < 1e-5) && (r > -1e-5) && (l < 1e-5) && (l > -1e-5)) ok = "OK"; printf "%s (R %+.2e, L %+.2e)", ok, r, l }')
echo "  order 3, magnetic conductor (mur = 50), 100 Hz : R, L == closed form within 1e-5 : $res"
case "$res" in NG*) status=1 ;; esac
# (g)
grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_ho3.msh"; next} /^frequency/{print "frequency = 1e2"; next} /^analysis = A/{print "fieldout = 1"} {print}' > "$WORK/pr_fld.ofe"
(cd "$WORK" && "$OFE" -n 2 pr_fld.ofe > /dev/null 2>&1) || { echo "  pr_fld : the run failed -> NG" >&2; status=1; }
compare "order 3, field output : volume of the cells [m^3]" "$(vtk vol J_A_re)" 5.0e-10 1e-9
compare "order 3, field output : |J| max (sigma V / L) [A/m^2]" "$(vtk vmax J_A_re)" 2.9e10 1e-4
compare "order 3, field output : |J| min [A/m^2]" "$(vtk vmin J_A_re)" 2.9e10 1e-4
compare "order 3, field output : mean J_x (signed) [A/m^2]" "$(vtk mean0 J_A_re)" -2.9e10 1e-4
pf=$(awk -v r="$(vtk int2 J_A_re)" -v m="$(vtk int2 J_A_im)" 'BEGIN{ printf "%.10e", (r + m) / (2 * 5.8e7) }')
pt=$(awk '/port 1 : I =/{ printf "%.10e", 0.5 * $6 }' "$WORK/ofe.log")
compare "order 3, field output : ohmic loss from the field == 1/2 Re(V I) [W]" "$pf" "$pt" 1e-6
# (h)
for m in box_p3 box_p4 box_p5 box_p6 box_p7; do cp "$SRC/$m.msh" "$WORK/"; done
cp "$SRC/edge_test.ofe" "$SRC/edge_test_aniso.ofe" "$WORK/"
awk 'BEGIN { a = 0.5; b = 0.7; ca = cos(a); sa = sin(a); cb = cos(b); sb = sin(b) }
	/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; x1 = ca * x - sa * y; y1 = sa * x + ca * y; y2 = cb * y1 - sb * z; z2 = sb * y1 + cb * z; printf "%s %.17g %.17g %.17g\n", $1, x1, y2, z2; next }
	{ print }' "$WORK/box_p3.msh" > "$WORK/box_p3_rot.msh"
for c in "edge_test box_p3 3" "edge_test_aniso box_p3 3" "edge_test box_p3_rot 3" "edge_test box_p4 4" "edge_test box_p5 5" "edge_test box_p6 6" "edge_test box_p7 7"; do
	set -- $c
	sed "s/^mesh = .*/mesh = $2.msh/" "$WORK/$1.ofe" > "$WORK/pre_$1_$2.ofe"
	(cd "$WORK" && "$OFE" -n 2 "pre_$1_$2.ofe" > /dev/null 2>&1) || true
	np=$(grep -c "^  (Pr[ab]) .* passed$" "$WORK/ofe.log" || true)
	if [ "$np" = 2 ] && grep -q "order-$3 edge self test : passed" "$WORK/ofe.log"; then
		echo "  order-$3 edge self test ($1 on $2) : (Pra) degree r-1, (Prb) + r x q : passed -> OK"
	else
		echo "  order-$3 edge self test ($1 on $2) : $np of 2 identities passed -> NG" >&2; status=1
		grep "(Pr\|\*\*\*" "$WORK/ofe.log" | cut -c1-160 >&2
	fi
done
cp "$SRC/coax_p3.msh" "$WORK/"
sed "s/^mesh = .*/mesh = coax_p3.msh/" "$WORK/edge_test.ofe" > "$WORK/prx_ecurv.ofe"
mesh_reject "order-3 edge self test on a curved mesh" prx_ecurv.ofe
grep -q "needs straight-sided elements" "$WORK/ofe.log" || { echo "  *** it was rejected for a different reason than: needs straight-sided elements" >&2; status=1; }
# (i)
grep -v '^#' "$WORK/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_ho3.msh"; next} {print}' > "$WORK/prbase.ofe"
awk '/^analysis/{print "gauge = 1"} {print}' "$WORK/prbase.ofe" > "$WORK/prx_g.ofe"
th_reject "order-3 A with gauge = 1" prx_g.ofe "supports gauge = 0 and the iterative solver only"
awk '/^region = 1 2/{print; print "bh = 2 100 0.5"; print "bh = 2 1000 1.0"; print "bh = 2 10000 1.2"; print "current = 1"; next} {print}' "$WORK/prbase.ofe" > "$WORK/prx_bh.ofe"
th_reject "order-3 A with a B-H curve" prx_bh.ofe "supports linear materials"

# ---- 3 次元渦電流の過渡解析 (analysis = A + transient、後退 Euler、電圧駆動) ----
#
# χ = Δt φ を未知数にすると係数行列が実対称・半正定値になる (eddy3d.c の solve_eddy3d_transient)。
#   (a) 電圧の段差を十分長く保つと I -> V / R_dc。DC の極限の節点の式は analysis R の導電行列
#       そのもの (Whitney の恒等式) なので、R の G に反復の打ち切りの範囲で一致する
#   (b) 刻みの仕事の恒等式 : Σ V I Δt = 磁気の仕事 Σ (ΔA)ᵀ S A^{n+1} + Joule 損 Σ zᵀ M z (厳密)
#   (c) 正弦波の定常 : 最後の周期の I の sin / cos 成分が周波数領域 (analysis A) の Y = G + jB に
#       後退 Euler の 1 次で近づく (1 周期 32 -> 64 点で誤差の比 1.9〜2.0)。Richardson 補外
#       2 I_64 - I_32 は |Y| の 1e-3 以内 (実測 9.5e-4)
echo "[A transient] 3-D eddy current in the time domain (backward Euler, voltage drive)"
cp "$SRC/bar_eddy.ofe" "$SRC/bar_tet.msh" "$WORK/"
a3_energy() {	# a3_energy <ラベル>
	res=$(awk '/^  energy : terminals/ { a = $5 + 0; b = $9 + 0; c = $12 + 0 }
		END { if (a == 0) { printf "NG (no energy line)"; exit }
		      d = (a - b - c) / a; if (d < 0) d = -d
		      printf "%s (terminals %.6e, magnetic %.6e, joule %.6e, %.1e)", (d < 1e-8) ? "OK" : "NG", a, b, c, d }' "$WORK/ofe.log")
	echo "  $1 : terminal work = magnetic work + joule loss : $res"
	case "$res" in NG*) status=1 ;; esac
}
# (a) (b)
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^analysis/{print "transient = 2e-4 100"} {print}' > "$WORK/a3t_step.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t_step.ofe > /dev/null 2>&1) || { echo "  a3t_step : the run failed -> NG" >&2; status=1; }
a3_energy "voltage step"
ilast=$(awk -F, 'NR > 1 { i = $3 } END { print i }' "$WORK/ofe_transient.csv")
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^awall/{next} /^analysis/{print "analysis = R"; next} {print}' > "$WORK/a3t_r.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t_r.ofe > /dev/null 2>&1 && "$OFE_POST" > /dev/null 2>&1) || { echo "  a3t_r : the run failed -> NG" >&2; status=1; }
gdc=$(awk -F, '$1 == "G" { f = 1; next } f && $1 == "1" { print $2; exit }' "$WORK/rlc.csv")
compare "voltage step held : I(t_end) == G_dc V (analysis R on the same mesh) [A]" "$ilast" "$gdc" 1e-9
# (c) 10 kHz の正弦波、4 周期
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{print "frequency = 1e4"; next} {print}' > "$WORK/a3t_fd.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t_fd.ofe > /dev/null 2>&1) || { echo "  a3t_fd : the run failed -> NG" >&2; status=1; }
set -- $(awk '/port 1 : I =/ { print $6, $7; exit }' "$WORK/ofe.log")
yre=$1; yim=$2
for N in 32 64; do
	awk -v N="$N" 'BEGIN { pi = atan2(0, -1); f = 1e4; T = 1 / f; s = "waveform ="
		for (k = 0; k <= 4 * N; k++) { t = k * T / N; s = s sprintf(" %.12e %.15g", t, sin(2 * pi * f * t)) } print s }' > "$WORK/a3t_wf.txt"
	grep -v '^#' "$SRC/bar_eddy.ofe" | awk -v wf="$WORK/a3t_wf.txt" -v ns=$((4 * N)) '/^frequency/{next}
		/^analysis/{ while ((getline l < wf) > 0) print l; print "transient = 4e-4 " ns } {print}' > "$WORK/a3t_sin.ofe"
	(cd "$WORK" && "$OFE" -n 2 a3t_sin.ofe > /dev/null 2>&1) || { echo "  a3t_sin ($N) : the run failed -> NG" >&2; status=1; }
	[ "$N" = 64 ] && a3_energy "sine, 64 steps per period"
	set -- $(awk -F, -v N="$N" 'BEGIN { pi = atan2(0, -1); f = 1e4 } NR > 1 { k = NR - 1
		if (k > 3 * N) { a += $3 * sin(2 * pi * f * $1); b += $3 * cos(2 * pi * f * $1) } } END { printf "%.12e %.12e", 2 * a / N, 2 * b / N }' "$WORK/ofe_transient.csv")
	eval "ar$N=\$1; ai$N=\$2"
done
res=$(awk -v gr="$yre" -v gi="$yim" -v a1="$ar32" -v b1="$ai32" -v a2="$ar64" -v b2="$ai64" 'BEGIN {
	y = sqrt(gr * gr + gi * gi)
	e1 = sqrt((a1 - gr) ^ 2 + (b1 - gi) ^ 2) / y; e2 = sqrt((a2 - gr) ^ 2 + (b2 - gi) ^ 2) / y
	ra = 2 * a2 - a1; rb = 2 * b2 - b1; er = sqrt((ra - gr) ^ 2 + (rb - gi) ^ 2) / y
	q = (e2 > 0) ? e1 / e2 : 0; ok = (q >= 1.7) && (q <= 2.3) && (er <= 2e-3)
	printf "%s (32 steps %.2e, 64 steps %.2e, ratio %.2f; Richardson %.2e)", ok ? "OK" : "NG", e1, e2, q, er }')
echo "  sine at 10 kHz : steady state -> frequency-domain Y (backward Euler, first order) : $res"
case "$res" in NG*) status=1 ;; esac
# (d) 弾くべき入力
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^analysis/{print "transient = 2e-4 10 0.5"} {print}' > "$WORK/a3tx_cn.ofe"
th_reject "A + transient with theta = 0.5" a3tx_cn.ofe "theta = 1 (backward Euler) only"
# (e) 10 節点四面体 (bar_tet2.msh) : 2 次の Nedelec の過渡解析 (eddy3d_p2.c の solve_eddy3d_p2_transient)。
#     同じ 3 つの検査 : 段差を保つと I = G_dc V (G = 7250 S、1e-9)、仕事の恒等式、10 kHz の正弦波の
#     2 周期目が 2 次の周波数領域の Y に 1 次で近づく (1 周期 16 -> 32 点で比 1.7〜2.3、Richardson 1e-2 以内、
#     実測 6.3e-2 -> 3.3e-2、比 1.93、Richardson 3.5e-3。CI の時間のため 1 次より刻みを粗く周期を短くした)
cp "$SRC/bar_tet2.msh" "$WORK/"
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_tet2.msh"; next} {print}' > "$WORK/a3t2_base.ofe"
awk '/^frequency/{next} /^analysis/{print "transient = 2e-4 50"} {print}' "$WORK/a3t2_base.ofe" > "$WORK/a3t2_step.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t2_step.ofe > /dev/null 2>&1) || { echo "  a3t2_step : the run failed -> NG" >&2; status=1; }
a3_energy "order 2, voltage step"
compare "order 2, voltage step held : I(t_end) == G_dc V [A]" "$(awk -F, 'NR > 1 { i = $3 } END { print i }' "$WORK/ofe_transient.csv")" 7250 1e-9
awk '/^frequency/{print "frequency = 1e4"; next} {print}' "$WORK/a3t2_base.ofe" > "$WORK/a3t2_fd.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t2_fd.ofe > /dev/null 2>&1) || { echo "  a3t2_fd : the run failed -> NG" >&2; status=1; }
set -- $(awk '/port 1 : I =/ { print $6, $7; exit }' "$WORK/ofe.log")
yre=$1; yim=$2
for N in 16 32; do
	awk -v N="$N" 'BEGIN { pi = atan2(0, -1); f = 1e4; T = 1 / f; s = "waveform ="
		for (k = 0; k <= 2 * N; k++) { t = k * T / N; s = s sprintf(" %.12e %.15g", t, sin(2 * pi * f * t)) } print s }' > "$WORK/a3t2_wf.txt"
	awk -v wf="$WORK/a3t2_wf.txt" -v ns=$((2 * N)) '/^frequency/{next}
		/^analysis/{ while ((getline l < wf) > 0) print l; print "transient = 2e-4 " ns } {print}' "$WORK/a3t2_base.ofe" > "$WORK/a3t2_sin.ofe"
	(cd "$WORK" && "$OFE" -n 2 a3t2_sin.ofe > /dev/null 2>&1) || { echo "  a3t2_sin ($N) : the run failed -> NG" >&2; status=1; }
	[ "$N" = 32 ] && a3_energy "order 2, sine, 32 steps per period"
	set -- $(awk -F, -v N="$N" 'BEGIN { pi = atan2(0, -1); f = 1e4 } NR > 1 { k = NR - 1
		if (k > N) { a += $3 * sin(2 * pi * f * $1); b += $3 * cos(2 * pi * f * $1) } } END { printf "%.12e %.12e", 2 * a / N, 2 * b / N }' "$WORK/ofe_transient.csv")
	eval "br$N=\$1; bi$N=\$2"
done
res=$(awk -v gr="$yre" -v gi="$yim" -v a1="$br16" -v b1="$bi16" -v a2="$br32" -v b2="$bi32" 'BEGIN {
	y = sqrt(gr * gr + gi * gi)
	e1 = sqrt((a1 - gr) ^ 2 + (b1 - gi) ^ 2) / y; e2 = sqrt((a2 - gr) ^ 2 + (b2 - gi) ^ 2) / y
	ra = 2 * a2 - a1; rb = 2 * b2 - b1; er = sqrt((ra - gr) ^ 2 + (rb - gi) ^ 2) / y
	q = (e2 > 0) ? e1 / e2 : 0; ok = "NG"; if ((q >= 1.7) && (q <= 2.3) && (er <= 1e-2)) ok = "OK"
	printf "%s (16 steps %.2e, 32 steps %.2e, ratio %.2f; Richardson %.2e)", ok, e1, e2, q, er }')
echo "  order 2, sine at 10 kHz : steady state -> frequency-domain Y (backward Euler, first order) : $res"
case "$res" in NG*) status=1 ;; esac
# (f) 3 次の四面体 (bar_ho3.msh) : 次数 r の Nedelec の過渡解析 (eddy3d_pr.c の solve_eddy3d_pr_transient)。
#     (e) と同じ 3 つの検査。電極の**全節点**を χ = Δt V にして端子電流もその行から取る (頂点だけだと
#     DC の極限が崩れる)。10 kHz の正弦波は空間の誤差が時間の誤差よりずっと小さいので 2 次と同じ数字になる
#     (実測 6.3e-2 -> 3.3e-2、比 1.93、Richardson 3.5e-3)。一時的に 2 次の格子をこの経路に通すと
#     eddy3d_p2.c の I(t) と ofe_transient.csv の桁 (1.5e-10) で一致した (同じ空間)。変異: 右辺を引き忘れる・
#     M の T/Δt や φ の行の結合を落とす・χ = V (Δt を掛け忘れ) は DC の極限で、電圧を刻みの始めで取る誤りは
#     正弦波でだけ落ちる (2 次・1 次と同じ)
cp "$SRC/bar_ho3.msh" "$WORK/"
awk '/^mesh/{print "mesh = bar_ho3.msh"; next} {print}' "$WORK/a3t2_base.ofe" > "$WORK/a3t3_base.ofe"
awk '/^frequency/{next} /^analysis/{print "transient = 2e-4 50"} {print}' "$WORK/a3t3_base.ofe" > "$WORK/a3t3_step.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t3_step.ofe > /dev/null 2>&1) || { echo "  a3t3_step : the run failed -> NG" >&2; status=1; }
a3_energy "order 3, voltage step"
compare "order 3, voltage step held : I(t_end) == G_dc V [A]" "$(awk -F, 'NR > 1 { i = $3 } END { print i }' "$WORK/ofe_transient.csv")" 7250 1e-9
awk '/^frequency/{print "frequency = 1e4"; next} {print}' "$WORK/a3t3_base.ofe" > "$WORK/a3t3_fd.ofe"
(cd "$WORK" && "$OFE" -n 2 a3t3_fd.ofe > /dev/null 2>&1) || { echo "  a3t3_fd : the run failed -> NG" >&2; status=1; }
set -- $(awk '/port 1 : I =/ { print $6, $7; exit }' "$WORK/ofe.log")
yre=$1; yim=$2
for N in 16 32; do
	awk -v N="$N" 'BEGIN { pi = atan2(0, -1); f = 1e4; T = 1 / f; s = "waveform ="
		for (k = 0; k <= 2 * N; k++) { t = k * T / N; s = s sprintf(" %.12e %.15g", t, sin(2 * pi * f * t)) } print s }' > "$WORK/a3t3_wf.txt"
	awk -v wf="$WORK/a3t3_wf.txt" -v ns=$((2 * N)) '/^frequency/{next}
		/^analysis/{ while ((getline l < wf) > 0) print l; print "transient = 2e-4 " ns } {print}' "$WORK/a3t3_base.ofe" > "$WORK/a3t3_sin.ofe"
	(cd "$WORK" && "$OFE" -n 2 a3t3_sin.ofe > /dev/null 2>&1) || { echo "  a3t3_sin ($N) : the run failed -> NG" >&2; status=1; }
	[ "$N" = 32 ] && a3_energy "order 3, sine, 32 steps per period"
	set -- $(awk -F, -v N="$N" 'BEGIN { pi = atan2(0, -1); f = 1e4 } NR > 1 { k = NR - 1
		if (k > N) { a += $3 * sin(2 * pi * f * $1); b += $3 * cos(2 * pi * f * $1) } } END { printf "%.12e %.12e", 2 * a / N, 2 * b / N }' "$WORK/ofe_transient.csv")
	eval "cr$N=\$1; ci$N=\$2"
done
res=$(awk -v gr="$yre" -v gi="$yim" -v a1="$cr16" -v b1="$ci16" -v a2="$cr32" -v b2="$ci32" 'BEGIN {
	y = sqrt(gr * gr + gi * gi)
	e1 = sqrt((a1 - gr) ^ 2 + (b1 - gi) ^ 2) / y; e2 = sqrt((a2 - gr) ^ 2 + (b2 - gi) ^ 2) / y
	ra = 2 * a2 - a1; rb = 2 * b2 - b1; er = sqrt((ra - gr) ^ 2 + (rb - gi) ^ 2) / y
	q = (e2 > 0) ? e1 / e2 : 0; ok = "NG"; if ((q >= 1.7) && (q <= 2.3) && (er <= 1e-2)) ok = "OK"
	printf "%s (16 steps %.2e, 32 steps %.2e, ratio %.2f; Richardson %.2e)", ok, e1, e2, q, er }')
echo "  order 3, sine at 10 kHz : steady state -> frequency-domain Y (backward Euler, first order) : $res"
case "$res" in NG*) status=1 ;; esac
# 3 次の過渡解析は線形材料だけ (B-H / J-A は入力エラー)
awk '/^region = 1 2/{print; print "bh = 2 100 0.5"; print "bh = 2 1000 1.0"; print "bh = 2 10000 1.2"; next} {print}' "$WORK/a3t3_step.ofe" > "$WORK/a3t3x_bh.ofe"
th_reject "order-3 A transient with a B-H curve" a3t3x_bh.ofe "supports linear materials"

# ---- 3 次元渦電流 + ヒステリシスの過渡解析 (analysis = A + transient + ja、ベクトル J-A、四面体) ----
#
# 各刻みを固定点反復 (行列は dH/dB の対称部分のテンソル、残差は確定済みの状態から積分した H) で解く。
#   (a) ほぼ線形の J-A (c -> 1、k 大、|He| < 1e-4 a で Langevin が線形) は μr = 1 + c Ms/(3a) の線形の
#       過渡解析に反復の打ち切りの範囲で一致し、行列の傾きが厳密なので 1 刻み 2 回で収束する
#       (曲率・体積・残差の組み立ての検査。傾きの誤りは答えを変えず回数だけ増やす)
#   (b) 準静的 (σ を下げて 1 Hz、t/δ ~ 0.02) の正弦波の 2 周期目のヒステリシス損 Σ V_e H・ΔB が、
#       Ampere の H_y(z) = J (z − t/2) の各点で H 駆動の順方向の J-A を awk で独立に積分した損
#       (ja_bar_loop.awk、FEM と同じ刻みで H を取り同じ右端の和) に 0.5% 以内 (実測 −0.07%)。
#       状態の確定・反転の扱い・H の積分の向きの誤りはここで落ちる (線形の (a) は履歴を見ない)
#   (c) 電圧の段差を保つと I -> G_dc V (磁気が DC を壊さない)、(d) 仕事の恒等式
#   (e) 3 成分のベクトル J-A を向きを固定した経路で動かすとスカラーの J-A と H が丸めまで一致
#       (起動時の自己検証。Heun の刻みの次数の誤りは損では 1e-6 しか動かずここでだけ落ちる)
#   (f) 剛体回転した格子で I(t) と固定点の回数が変わらないこと (行列のテンソルの非対角の検査)
echo "[A + ja] 3-D eddy current with hysteresis (vector Jiles-Atherton on tetrahedra, backward Euler)"
# (a)
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^voltage/{print "voltage = 1e-4"; next}
	/^analysis/{print "transient = 2e-4 50"; print "mur = 2 50"} {print}' > "$WORK/aja_lin.ofe"
(cd "$WORK" && "$OFE" -n 2 aja_lin.ofe > /dev/null 2>&1 && cp ofe_transient.csv aja_lin.csv) || { echo "  aja_lin : the run failed -> NG" >&2; status=1; }
jms=$(awk 'BEGIN { printf "%.15g", 49 * 3e8 / 0.999999 }')
grep -v '^#' "$SRC/bar_eddy.ofe" | awk -v ms="$jms" '/^frequency/{next} /^voltage/{print "voltage = 1e-4"; next}
	/^analysis/{print "transient = 2e-4 50"; print "ja = 2 " ms " 1e8 0 1e6 0.999999"; print "nlsolver = 50 1e-9 1"} {print}' > "$WORK/aja_nl.ofe"
(cd "$WORK" && "$OFE" -n 2 aja_nl.ofe > /dev/null 2>&1) || { echo "  aja_nl : the run failed -> NG" >&2; status=1; }
res=$(paste -d, "$WORK/aja_lin.csv" "$WORK/ofe_transient.csv" | awk -F, 'NR > 1 { d = ($3 - $10) / $3; if (d < 0) d = -d; if (d > m) m = d; n++ }
	END { ok = "NG"; if ((n == 50) && (m < 1e-7)) ok = "OK"; printf "%s (%d steps, max %.1e)", ok, n, m + 0 }')
echo "  near-linear J-A (c = 1 - 1e-6, |He| < 1e-4 a) == linear mur = 1 + c Ms/(3a) : I(t) : $res"
case "$res" in NG*) status=1 ;; esac
nfp=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")
if [ "${nfp:-0}" -ge 50 ] && [ "${nfp:-0}" -le 100 ]; then
	echo "  near-linear J-A : fixed-point iterations = $nfp (exact tangent : at most 2 per step) -> OK"
else
	echo "  near-linear J-A : fixed-point iterations = ${nfp:-none} (expected 50..100) -> NG"; status=1
fi
a3_energy "near-linear J-A"
# (b) 2 周期、1 周期 32 刻み
awk 'BEGIN { pi = atan2(0, -1); s = "waveform ="
	for (k = 0; k <= 64; k++) { t = k / 32; s = s sprintf(" %.12e %.15g", t, sin(2 * pi * t)) } print s }' > "$WORK/aja_wf.txt"
grep -v '^#' "$SRC/bar_eddy.ofe" | awk -v wf="$WORK/aja_wf.txt" '/^frequency/{next} /^material/{print "material = 1.0 1e5"; next}
	/^voltage/{print "voltage = 0.04"; next}
	/^analysis/{ while ((getline l < wf) > 0) print l; print "transient = 2 64"; print "ja = 2 1.6e6 1100 1.6e-3 400 0.2"; print "nlsolver = 200 1e-7 1" } {print}' > "$WORK/aja_qs.ofe"
(cd "$WORK" && "$OFE" -n 2 aja_qs.ofe > /dev/null 2>&1) || { echo "  aja_qs : the run failed -> NG" >&2; status=1; }
a3_energy "quasi-static sine, 2 periods"
wfem=$(awk -F, 'NR == 33 { w = $5 } END { printf "%.10e", $5 - w }' "$WORK/ofe_transient.csv")
# J0 = σ V0 / L (準静的なので電流密度は一様)、面積 = L W
wref=$(awk -v ms=1.6e6 -v aa=1100 -v al=1.6e-3 -v kk=400 -v cc=0.2 -v j0=2e6 -v tt=1e-3 -v area=5e-7 -v ncyc=2 -v nstep=32 -f "$SRC/ja_bar_loop.awk")
compare "quasi-static hysteresis loss of the 2nd period == forward J-A along H_y(z) = J (z - t/2) (awk) [J]" "$wfem" "$wref" 5e-3
nfp=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")
if [ "${nfp:-0}" -ge 64 ] && [ "${nfp:-0}" -le 1100 ]; then
	echo "  quasi-static sine : fixed-point iterations = $nfp (64 steps, at most 1100) -> OK"
else
	echo "  quasi-static sine : fixed-point iterations = ${nfp:-none} (expected 64..1100) -> NG"; status=1
fi
#   (g) 異方性の J-A (jaaxis) を粗い棒 (bar_tet_c.msh = mkmesh.py bar -nx 8 -ny 1 -nz 8) で:
#       全軸を ja と同じにすると I(t) が等方性と一致、B の向き (y) の軸だけを別のパラメータにすると損が
#       その軸のパラメータの awk に 1% (粗い格子の誤差。等方性でも −0.54%)、座標を巡回させた格子と
#       軸のパラメータを同じく巡回させた入力で I(t) と反復回数が一致 (成分と軸の対応の検査)
cp "$SRC/bar_tet_c.msh" "$WORK/"
awk '/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 } f && NF == 4 { printf "%s %.17g %.17g %.17g\n", $1, $4, $2, $3; next } { print }' \
	"$WORK/bar_tet_c.msh" > "$WORK/bar_tet_c_perm.msh"
jax() {	# jax <名前> <格子> <x の行> <y の行> <z の行> : aja_qs を粗い格子と jaaxis で解いて csv を残す
	awk -v m="$2" -v px="$3" -v py="$4" -v pz="$5" '/^mesh/{print "mesh = " m; next}
		/^ja =/{print; if (px != "") print "jaaxis = 2 x " px; if (py != "") print "jaaxis = 2 y " py; if (pz != "") print "jaaxis = 2 z " pz; next} {print}' \
		"$WORK/aja_qs.ofe" > "$WORK/jax_$1.ofe"
	(cd "$WORK" && "$OFE" -n 2 "jax_$1.ofe" > /dev/null 2>&1 && cp ofe_transient.csv "jax_$1.csv" && cp ofe.log "jax_$1.log") || { echo "  jax_$1 : the run failed -> NG" >&2; status=1; }
}
jaP="1.6e6 1100 1.6e-3 400 0.2"; jaQ="1.2e6 800 1.0e-3 250 0.3"; jaP1="2.0e6 1300 2.0e-3 500 0.15"; jaP2="1.0e6 900 1.2e-3 300 0.25"
jax iso bar_tet_c.msh "" "" ""
jax eq bar_tet_c.msh "$jaP" "$jaP" "$jaP"
jax ay bar_tet_c.msh "" "$jaQ" ""
jax a3 bar_tet_c.msh "$jaP1" "$jaQ" "$jaP2"
jax a3p bar_tet_c_perm.msh "$jaP2" "$jaP1" "$jaQ"
jax_same() {	# jax_same <ラベル> <a> <b> : I(t) と反復回数が一致
	res=$(paste -d, "$WORK/jax_$2.csv" "$WORK/jax_$3.csv" | awk -F, -v na="$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/jax_$2.log")" \
		-v nb="$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/jax_$3.log")" '
		NR > 1 { d = ($3 - $10) / 0.5; if (d < 0) d = -d; if (d > m) m = d; n++ }
		END { ok = "NG"; if ((n == 64) && (m < 1e-9) && (na + 0 > 0) && (na == nb)) ok = "OK"
		      printf "%s (max %.1e of the peak, iterations %d / %d)", ok, m + 0, na, nb }')
	echo "  anisotropic J-A : $1 : $res"
	case "$res" in NG*) status=1 ;; esac
}
jax_same "every axis equal to ja == isotropic" iso eq
jax_same "cyclically permuted mesh and axis parameters == original" a3 a3p
wfem=$(awk -F, 'NR == 33 { w = $5 } END { printf "%.10e", $5 - w }' "$WORK/jax_ay.csv")
set -- $jaQ
wq=$(awk -v ms="$1" -v aa="$2" -v al="$3" -v kk="$4" -v cc="$5" -v j0=2e6 -v tt=1e-3 -v area=5e-7 -v ncyc=2 -v nstep=32 -f "$SRC/ja_bar_loop.awk")
compare "anisotropic J-A, B along y with the y-axis parameters : loss == scalar J-A of those parameters (awk) [J]" "$wfem" "$wq" 1e-2
wiso=$(awk -F, 'NR == 33 { w = $5 } END { printf "%.10e", $5 - w }' "$WORK/jax_iso.csv")
res=$(awk -v a="$wfem" -v b="$wiso" 'BEGIN { d = (a - b) / b; if (d < 0) d = -d; printf "%s (%.0f%%)", (d > 0.1) ? "OK" : "NG", 100 * d }')
echo "  anisotropic J-A : the y-axis parameters change the loss from the isotropic one : $res"
case "$res" in NG*) status=1 ;; esac
# (c) 段差を保つ (σ = 1e5、J-A の μr ~ 1000 で τ ~ 1.3e-5 s。20 刻み x 1e-4 s で DC)
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^material/{print "material = 1.0 1e5"; next} /^voltage/{print "voltage = 0.04"; next}
	/^analysis/{print "transient = 2e-3 20"; print "ja = 2 1.6e6 1100 1.6e-3 400 0.2"; print "nlsolver = 200 1e-9 1"} {print}' > "$WORK/aja_dc.ofe"
(cd "$WORK" && "$OFE" -n 2 aja_dc.ofe > aja_dc.out 2>&1 && cp ofe_transient.csv aja_dc.csv && cp ofe.log aja_dc.log) || { echo "  aja_dc : the run failed -> NG" >&2; status=1; }
# 未使用キーの警告は両方向 : 正しいケースで 0 件、javector (3 次元は常にベクトル J-A) で 1 件
nw0=$(grep -c "\*\*\* warning" "$WORK/aja_dc.out" || true)
awk '/^analysis/{print "javector = 1"} {print}' "$WORK/aja_dc.ofe" > "$WORK/aja_jv.ofe"
nw1=$(cd "$WORK" && "$OFE" -n 2 aja_jv.ofe 2>&1 | grep -c "javector has no effect with analysis A" || true)
if [ "$nw0" = 0 ] && [ "$nw1" = 1 ]; then
	echo "  input lint : no warning for A + ja, one for javector with A -> OK"
else
	echo "  input lint : warnings A + ja = $nw0 (expected 0), javector with A = $nw1 (expected 1) -> NG"; status=1
fi
a3_energy "voltage step with hysteresis"
ilast=$(awk -F, 'NR > 1 { i = $3 } END { print i }' "$WORK/ofe_transient.csv")
idc=$(awk -v g="$gdc" 'BEGIN { printf "%.12e", g * (1e5 / 5.8e7) * 0.04 }')
compare "voltage step held with hysteresis : I(t_end) == G_dc V [A]" "$ilast" "$idc" 1e-9
res=$(awk '/J-A model check/ { v = $NF + 0; n++ } END { printf "%s (%.1e)", ((n == 1) && (v < 1e-12)) ? "OK" : "NG", v }' "$WORK/aja_dc.log")
echo "  3-D vector J-A along a fixed direction == scalar J-A (same Heun steps, H to rounding) : $res"
case "$res" in NG*) status=1 ;; esac
# (f) 剛体回転した格子 : 物理は同じなので I(t) が一致し、固定点の回数も増えない。軸に平行な格子では
#     B が 1 軸に揃って行列のテンソルの非対角が 0 になり、非対角の誤り (符号・添字) は答えも回数も
#     変えない (実測で素通り)。回すと非対角が効くので、誤りは回数の増加 (か非収束) で見える
awk 'BEGIN { a = 0.5; b = 0.7; ca = cos(a); sa = sin(a); cb = cos(b); sb = sin(b) }
	/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; x1 = ca * x - sa * y; y1 = sa * x + ca * y; y2 = cb * y1 - sb * z; z2 = sb * y1 + cb * z; printf "%s %.17g %.17g %.17g\n", $1, x1, y2, z2; next }
	{ print }' "$SRC/bar_tet.msh" > "$WORK/bar_tet_rot.msh"
nfp0=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/aja_dc.log")
sed 's/^mesh = .*/mesh = bar_tet_rot.msh/' "$WORK/aja_dc.ofe" > "$WORK/aja_dcr.ofe"
(cd "$WORK" && "$OFE" -n 2 aja_dcr.ofe > /dev/null 2>&1) || { echo "  aja_dcr : the run failed -> NG" >&2; status=1; }
res=$(paste -d, "$WORK/aja_dc.csv" "$WORK/ofe_transient.csv" | awk -F, 'NR > 1 { d = ($3 - $10) / $3; if (d < 0) d = -d; if (d > m) m = d; n++ }
	END { ok = "NG"; if ((n == 20) && (m < 1e-6)) ok = "OK"; printf "%s (%d steps, max %.1e)", ok, n, m + 0 }')
echo "  rigidly rotated mesh (hysteresis, voltage step) : I(t) == axis-aligned mesh : $res"
case "$res" in NG*) status=1 ;; esac
nfp=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")
res=$(awk -v a="${nfp:-0}" -v b="${nfp0:-0}" 'BEGIN { printf "%s", ((a > 0) && (b > 0) && (a <= 1.2 * b + 5)) ? "OK" : "NG" }')
echo "  rigidly rotated mesh : fixed-point iterations = ${nfp:-none} (axis-aligned ${nfp0:-none}, at most 1.2x + 5) -> $res"
case "$res" in NG*) status=1 ;; esac
# (e) 弾くべき入力
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^analysis/{print "ja = 2 1.6e6 1100 1.6e-3 400 0.2"} {print}' > "$WORK/ajax_fd.ofe"
th_reject "ja with analysis A in the frequency domain" ajax_fd.ofe "needs transient"
cp "$SRC/bar_prismtet.ofe" "$SRC/bar_prismtet.msh" "$WORK/"
grep -v '^#' "$SRC/bar_prismtet.ofe" | awk '/^frequency/{next} /^analysis/{print "transient = 2e-4 5"; print "ja = 2 1.6e6 1100 1.6e-3 400 0.2"} {print}' > "$WORK/ajax_pt.ofe"
if (cd "$WORK" && "$OFE" -n 2 ajax_pt.ofe > /dev/null 2>&1); then
	echo "  ja on prisms with A + transient : not rejected -> NG"; status=1
elif grep -q "needs first-order tetrahedra" "$WORK/ofe.log"; then
	echo "  ja on prisms with A + transient : rejected -> OK"
else
	echo "  ja on prisms with A + transient : rejected for a different reason -> NG"; status=1
fi

# ---- 3 次元渦電流の過渡解析 + B-H 曲線 (analysis = A + transient + bh、四面体) ----
#
# J-A と同じ固定点反復で、H = ν(|B|) B、行列は厳密なヤコビアン ν I + 2 (dν/dB²) B Bᵀ。
#   (a) 直線の曲線 (1 点) は μr の線形の過渡解析に一致し、1 刻み 2 回で収束する
#   (b) 準静的な正弦波の 2 周期目の Σ H^{n+1}ΔB (履歴が無いので連続の極限では 0、右端の和の分だけ残る) が、
#       Ampere の H_y(z) = J (z − t/2) の各点で曲線から B(H) を引いた awk の同じ和 (ja_bar_loop.awk) に 0.5% 以内
#   (c) 剛体回転した格子で I(t) と反復回数が変わらない (ヤコビアンの rank-1 項の非対角の検査)
echo "[A + bh] 3-D eddy current with a B-H curve in the time domain (tetrahedra, exact Jacobian)"
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^voltage/{print "voltage = 1e-4"; next}
	/^analysis/{print "transient = 2e-4 50"; printf "bh = 2 100 %.17g\n", 100 * 4e-7 * atan2(0, -1) * 50; print "nlsolver = 50 1e-9 1"} {print}' > "$WORK/abh_lin.ofe"
(cd "$WORK" && "$OFE" -n 2 abh_lin.ofe > /dev/null 2>&1) || { echo "  abh_lin : the run failed -> NG" >&2; status=1; }
res=$(paste -d, "$WORK/aja_lin.csv" "$WORK/ofe_transient.csv" | awk -F, 'NR > 1 { d = ($3 - $10) / $3; if (d < 0) d = -d; if (d > m) m = d; n++ }
	END { ok = "NG"; if ((n == 50) && (m < 1e-7)) ok = "OK"; printf "%s (%d steps, max %.1e)", ok, n, m + 0 }')
echo "  straight-line B-H curve == linear mur = 50 : I(t) : $res"
case "$res" in NG*) status=1 ;; esac
nfp=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")
if [ "${nfp:-0}" -ge 50 ] && [ "${nfp:-0}" -le 100 ]; then
	echo "  straight-line B-H curve : fixed-point iterations = $nfp (at most 2 per step) -> OK"
else
	echo "  straight-line B-H curve : fixed-point iterations = ${nfp:-none} (expected 50..100) -> NG"; status=1
fi
awk '/^ja =/{print "bh = 2 100 0.5"; print "bh = 2 1000 1.0"; print "bh = 2 10000 1.2"; next} {print}' "$WORK/aja_qs.ofe" > "$WORK/abh_qs.ofe"
(cd "$WORK" && "$OFE" -n 2 abh_qs.ofe > /dev/null 2>&1 && cp ofe_transient.csv abh_qs.csv && cp ofe.log abh_qs.log) || { echo "  abh_qs : the run failed -> NG" >&2; status=1; }
a3_energy "B-H curve, quasi-static sine"
wfem=$(awk -F, 'NR == 33 { w = $5 } END { printf "%.10e", $5 - w }' "$WORK/abh_qs.csv")
wref=$(awk -v j0=2e6 -v tt=1e-3 -v area=5e-7 -v ncyc=2 -v nstep=32 -v bh="100 0.5 1000 1.0 10000 1.2" -f "$SRC/ja_bar_loop.awk")
compare "B-H curve, quasi-static sine : 2nd-period sum of H dB == the same sum along H_y(z) = J (z - t/2) (awk) [J]" "$wfem" "$wref" 5e-3
nfp0=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/abh_qs.log")
if [ "${nfp0:-0}" -ge 64 ] && [ "${nfp0:-0}" -le 400 ]; then
	echo "  B-H curve, quasi-static sine : fixed-point iterations = $nfp0 (64 steps, at most 400) -> OK"
else
	echo "  B-H curve, quasi-static sine : fixed-point iterations = ${nfp0:-none} (expected 64..400) -> NG"; status=1
fi
sed 's/^mesh = .*/mesh = bar_tet_rot.msh/' "$WORK/abh_qs.ofe" > "$WORK/abh_qsr.ofe"
(cd "$WORK" && "$OFE" -n 2 abh_qsr.ofe > /dev/null 2>&1) || { echo "  abh_qsr : the run failed -> NG" >&2; status=1; }
res=$(paste -d, "$WORK/abh_qs.csv" "$WORK/ofe_transient.csv" | awk -F, 'NR > 1 { d = ($3 - $10) / 0.5; if (d < 0) d = -d; if (d > m) m = d; n++ }
	END { ok = "NG"; if ((n == 64) && (m < 1e-6)) ok = "OK"; printf "%s (%d steps, max %.1e of the peak)", ok, n, m + 0 }')
echo "  B-H curve, rigidly rotated mesh : I(t) == axis-aligned mesh : $res"
case "$res" in NG*) status=1 ;; esac
nfp=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")
res=$(awk -v a="${nfp:-0}" -v b="${nfp0:-0}" 'BEGIN { printf "%s", ((a > 0) && (b > 0) && (a <= 1.2 * b + 5)) ? "OK" : "NG" }')
echo "  B-H curve, rigidly rotated mesh : fixed-point iterations = ${nfp:-none} (axis-aligned ${nfp0:-none}, at most 1.2x + 5) -> $res"
case "$res" in NG*) status=1 ;; esac
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^frequency/{next} /^analysis/{print "transient = 2e-4 5"; print "bh = 2 X 100 0.5"; print "bh = 2 Y 100 0.5"; print "bh = 2 Z 100 0.5"} {print}' > "$WORK/abhx_ax.ofe"
th_reject "per-axis bh curves with A + transient" abhx_ax.ofe "isotropic bh curves only"
# (d) 10 節点四面体 (bar_tet2c.msh = mkmesh.py bar -nx 8 -ny 1 -nz 8 を gmsh が 2 次に) : 2 次の Nedelec で
#     **積分点ごとの状態** (eddy3d_p2.c の p2nl_*)。同じ検査を同じ形で:
#     直線の B-H 曲線 / ほぼ線形の J-A == 線形の μr = 50 (1 刻み 2 回)、準静的な正弦波の B-H の和と J-A の損が
#     awk の独立な積分に 0.5% (実測 −0.08% / −0.04%)、仕事の恒等式、剛体回転した格子で I(t) と回数が不変。
#     J-A は CI の時間のため 1 周期 16 点 (awk も同じ刻み)。ヤコビアンの非対角の符号の誤りは回転した格子で
#     収束しなくなり、J-A の状態の確定漏れは損 −45%、仕事の H の取り違えは恒等式で落ちる (どれも実測)
cp "$SRC/bar_tet2c.msh" "$WORK/"
grep -v '^#' "$SRC/bar_eddy.ofe" | awk '/^mesh/{print "mesh = bar_tet2c.msh"; next} /^frequency/{next} /^voltage/{print "voltage = 1e-4"; next}
	/^analysis/{print "transient = 2e-4 50"; print "mur = 2 50"} {print}' > "$WORK/a2nl_lin.ofe"
(cd "$WORK" && "$OFE" -n 2 a2nl_lin.ofe > /dev/null 2>&1 && cp ofe_transient.csv a2nl_lin.csv) || { echo "  a2nl_lin : the run failed -> NG" >&2; status=1; }
a2nl_same() {	# a2nl_same <ラベル> <ofe> : a2nl_lin と I(t) が 1e-7、固定点は 1 刻み 2 回まで
	(cd "$WORK" && "$OFE" -n 2 "$2" > /dev/null 2>&1) || { echo "  $2 : the run failed -> NG" >&2; status=1; }
	res=$(paste -d, "$WORK/a2nl_lin.csv" "$WORK/ofe_transient.csv" | awk -F, -v nfp="$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")" '
		NR > 1 { d = ($3 - $10) / $3; if (d < 0) d = -d; if (d > m) m = d; n++ }
		END { ok = "NG"; if ((n == 50) && (m < 1e-7) && (nfp + 0 >= 50) && (nfp + 0 <= 100)) ok = "OK"
		      printf "%s (%d steps, max %.1e, %d fixed-point iterations)", ok, n, m + 0, nfp + 0 }')
	echo "  order 2, $1 == linear mur = 50 : I(t) : $res"
	case "$res" in NG*) status=1 ;; esac
}
awk '/^mur = 2 50/{ printf "bh = 2 100 %.17g\n", 100 * 4e-7 * atan2(0, -1) * 50; print "nlsolver = 50 1e-9 1"; next } {print}' "$WORK/a2nl_lin.ofe" > "$WORK/a2nl_bhl.ofe"
a2nl_same "straight-line B-H curve" a2nl_bhl.ofe
awk -v ms="$jms" '/^mur = 2 50/{ print "ja = 2 " ms " 1e8 0 1e6 0.999999"; print "nlsolver = 50 1e-9 1"; next } {print}' "$WORK/a2nl_lin.ofe" > "$WORK/a2nl_jal.ofe"
a2nl_same "near-linear J-A" a2nl_jal.ofe
sed 's/^mesh = .*/mesh = bar_tet2c.msh/' "$WORK/abh_qs.ofe" > "$WORK/a2nl_bqs.ofe"
(cd "$WORK" && "$OFE" -n 2 a2nl_bqs.ofe > /dev/null 2>&1 && cp ofe_transient.csv a2nl_bqs.csv && cp ofe.log a2nl_bqs.log) || { echo "  a2nl_bqs : the run failed -> NG" >&2; status=1; }
a3_energy "order 2, B-H curve, quasi-static sine"
wfem=$(awk -F, 'NR == 33 { w = $5 } END { printf "%.10e", $5 - w }' "$WORK/a2nl_bqs.csv")
compare "order 2, B-H curve, quasi-static sine : 2nd-period sum of H dB == the same sum along H_y(z) (awk) [J]" "$wfem" "$wref" 5e-3
awk 'BEGIN { a = 0.5; b = 0.7; ca = cos(a); sa = sin(a); cb = cos(b); sb = sin(b) }
	/^\$Nodes/ { f = 1; print; next } /^\$EndNodes/ { f = 0 }
	f && NF == 4 { x = $2; y = $3; z = $4; x1 = ca * x - sa * y; y1 = sa * x + ca * y; y2 = cb * y1 - sb * z; z2 = sb * y1 + cb * z; printf "%s %.17g %.17g %.17g\n", $1, x1, y2, z2; next }
	{ print }' "$WORK/bar_tet2c.msh" > "$WORK/bar_tet2c_rot.msh"
sed 's/^mesh = .*/mesh = bar_tet2c_rot.msh/' "$WORK/a2nl_bqs.ofe" > "$WORK/a2nl_bqsr.ofe"
(cd "$WORK" && "$OFE" -n 2 a2nl_bqsr.ofe > /dev/null 2>&1) || { echo "  a2nl_bqsr : the run failed -> NG" >&2; status=1; }
nfa=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/a2nl_bqs.log")
nfr=$(awk '/fixed-point iterations =/ { print $(NF - 2) }' "$WORK/ofe.log")
res=$(paste -d, "$WORK/a2nl_bqs.csv" "$WORK/ofe_transient.csv" | awk -F, -v a="${nfa:-0}" -v b="${nfr:-0}" '
	NR > 1 { d = ($3 - $10) / 0.5; if (d < 0) d = -d; if (d > m) m = d; n++ }
	END { ok = "NG"; if ((n == 64) && (m < 1e-6) && (a > 0) && (b > 0) && (b <= 1.2 * a + 5)) ok = "OK"
	      printf "%s (%d steps, max %.1e of the peak; iterations %d vs %d)", ok, n, m + 0, b, a }')
echo "  order 2, B-H curve, rigidly rotated mesh : I(t) and fixed-point iterations == axis-aligned mesh : $res"
case "$res" in NG*) status=1 ;; esac
awk 'BEGIN { pi = atan2(0, -1); s = "waveform ="
	for (k = 0; k <= 32; k++) { t = k / 16; s = s sprintf(" %.12e %.15g", t, sin(2 * pi * t)) } print s }' > "$WORK/a2nl_wf.txt"
awk -v wf="$WORK/a2nl_wf.txt" '/^waveform/{ while ((getline l < wf) > 0) print l; next } /^transient/{print "transient = 2 32"; next} {print}' \
	"$WORK/aja_qs.ofe" | sed 's/^mesh = .*/mesh = bar_tet2c.msh/' > "$WORK/a2nl_jqs.ofe"
(cd "$WORK" && "$OFE" -n 2 a2nl_jqs.ofe > /dev/null 2>&1) || { echo "  a2nl_jqs : the run failed -> NG" >&2; status=1; }
a3_energy "order 2, J-A, quasi-static sine (16 steps per period)"
wfem=$(awk -F, 'NR == 17 { w = $5 } END { printf "%.10e", $5 - w }' "$WORK/ofe_transient.csv")
wref=$(awk -v ms=1.6e6 -v aa=1100 -v al=1.6e-3 -v kk=400 -v cc=0.2 -v j0=2e6 -v tt=1e-3 -v area=5e-7 -v ncyc=2 -v nstep=16 -f "$SRC/ja_bar_loop.awk")
compare "order 2, J-A, quasi-static hysteresis loss of the 2nd period == forward J-A along H_y(z) (awk) [J]" "$wfem" "$wref" 5e-3

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
awk '!/^thermal/' "$SRC/thermal_bar.ofe" | awk '!/^heatsink/' > "$WORK/lint_kth.ofe"
awk '!/^thermal/' "$SRC/thermal_bar.ofe" | awk '!/^kthermal/' > "$WORK/lint_sink.ofe"
lint_expect "kthermal without thermal"     lint_kth.ofe       yes
lint_expect "heatsink without thermal"     lint_sink.ofe      yes
awk '/^analysis = /{print "convection = 10 20"} {print}' \
    "$SRC/resistor_bar.ofe" > "$WORK/lint_conv.ofe"
lint_expect "convection without thermal"   lint_conv.ofe      yes
awk '/^analysis = /{print "radiation = 0.5 20"} {print}' \
    "$SRC/resistor_bar.ofe" > "$WORK/lint_rad.ofe"
lint_expect "radiation without thermal"    lint_rad.ofe       yes
awk '/^thermal/{print "heatcap = 2 3.4e6"} {print}' \
    "$SRC/thermal_bar.ofe" > "$WORK/lint_heatcap.ofe"
lint_expect "heatcap without transient"    lint_heatcap.ofe   yes
sed "s/^epstempco = .*/mutemptable = 2 20 1 120 0.7/" "$SRC/temp_material.ofe" > "$WORK/lint_mutab.ofe"
lint_expect "mutemptable + analysis=C"     lint_mutab.ofe     yes
sed "s/^mutempco = .*/mutemptable = 2 20 1 120 0.7/" "$SRC/temp_mur.ofe" > "$WORK/lint_mutab_ok.ofe"
lint_expect "mutemptable + analysis=M (clean)" lint_mutab_ok.ofe no
sed "s/^epstempco = .*/tempcotable = 2 20 1 120 1.3/" "$SRC/temp_material.ofe" > "$WORK/lint_sigtab.ofe"
lint_expect "tempcotable + analysis=C"     lint_sigtab.ofe    yes
awk '/^analysis = /{print "conductortempcotable = 1 20 1 120 1.3"} {print}' \
    "$SRC/bar_eddy.ofe" > "$WORK/lint_ctab_a.ofe"
lint_expect "conductortempcotable + analysis=A" lint_ctab_a.ofe yes
awk '/^convection =/{next} /^thermal/{print "heatsink = 0 20"} {print}' "$SRC/thermal_natconv.ofe" > "$WORK/lint_ctt_noconv.ofe"
lint_expect "convectiontemptable without convection" lint_ctt_noconv.ofe yes
awk '/^thermal/{next} {print}' "$SRC/thermal_natconv.ofe" | awk '!/^kthermal|^convection = /' > "$WORK/lint_ctt_noth.ofe"
lint_expect "convectiontemptable without thermal" lint_ctt_noth.ofe yes
#   (thrt*.ofe は [thermal h(T)] 節で作ったもの)
lint_expect "radiationtemptable without radiation" thrtx_norad.ofe yes
lint_expect "radiationtemptable (clean)" thrt.ofe no
#   (thpf*.ofe は [thermal h(x)] 節で作ったもの)
lint_expect "convectionprofile (clean)" thpf_s.ofe no
lint_expect "radiationprofile (clean)" thpf_rs.ofe no
lint_expect "javector + currentsweep2 + japrobe (clean)" ja_rotating.ofe no
awk '/^ja =/{next} /^currentsweep/{next} /^jasub/{next} {print}' "$SRC/ja_rotating.ofe" | awk '!/^conductor = 2/' > "$WORK/lint_jv_noja.ofe"
lint_expect "javector / japrobe without ja" lint_jv_noja.ofe yes
awk '/^convection =/{next} /^thermal/{print "heatsink = 0 20"} {print}' "$WORK/thpf_s.ofe" > "$WORK/lint_pf_noconv.ofe"
lint_expect "convectionprofile without convection" lint_pf_noconv.ofe yes
awk '/^thermal/{next} {print}' "$WORK/thpf_s.ofe" | awk '!/^kthermal|^convection = /' > "$WORK/lint_pf_noth.ofe"
lint_expect "convectionprofile without thermal" lint_pf_noth.ofe yes
lint_expect "convectionexpr (clean)" thex_s.ofe no
lint_expect "radiationexpr (clean)" thex_rs.ofe no
awk '/^convection =/{next} /^thermal/{print "heatsink = 0 20"} {print}' "$WORK/thex_s.ofe" > "$WORK/lint_ex_noconv.ofe"
lint_expect "convectionexpr without convection" lint_ex_noconv.ofe yes
awk '/^thermal/{next} {print}' "$WORK/thex_s.ofe" | awk '!/^kthermal|^convection = /' > "$WORK/lint_ex_noth.ofe"
lint_expect "convectionexpr without thermal" lint_ex_noth.ofe yes
awk '/^thermal/{print "radiationexchange = 1"} {print}' "$SRC/thermal_bar.ofe" > "$WORK/lint_rx_norad.ofe"
lint_expect "radiationexchange without radiation" lint_rx_norad.ofe yes
# 正しいケースでは 1 件も出ないこと (誤検知は警告を無視させるので同じくらい悪い)
for c in parallel_plate resistor_bar coax microstrip plate_line_dc coax_loss \
         plate_line_ac plate_line_bh dispersive_plate drude_plate colecole_plate \
         temp_resistor aniso_plate box_tet coax_tet edge_test bar_eddy \
         box_p2 coax_p2 nodal_test_p2 plate2d_dc plate2d_ac plate2d_r30 \
         bertotti_core temp_material temp_mur hn_plate plate2d_p2 plate2d_p3 nodal_test_p4 thermal_f2d \
         box_tet_41 plate2d_41 coupled_microstrip thermal_bar thermal_neck \
         thermal_conv thermal_rad thermal_transient thermal_pulse thermal_plate \
         plate_line_bh_ac bar_air_bh thermal_tet thermal_ac thermal_natconv \
         box_p3 nodal_test_p3 coax_p3 thermal_cavity lamination_step plate_line_ja_eddy; do
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
