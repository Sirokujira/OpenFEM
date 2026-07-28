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
#   plate2d        : 断面 2 次元の三角形格子で M / F を解き、構造格子版と同じ
#                    1 次元厳密解と比較する (Ldc / Rs / R(f) / L(f)、許容 0.2%)
#   anisotropic mu : 面内の異方性 ν が **B に掛かる** こと (grad(Az) ではなく)。
#                    等方性では一致するので等方性ケースでは検出できない
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
vtk() { awk -v arr="$2" -v axis="${3:-0}" -v xcut="${4:-}" -v comp="${5:-}" \
	-f "$SRC/vtkcheck.awk" "$WORK/ofe_field.vtk" | awk -v k="$1" '$1 == k { print $2 }'; }

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
         box_p2 coax_p2 nodal_test_p2 plate2d_dc plate2d_ac; do
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

if [ "$status" -ne 0 ]; then
	echo "*** RLC validation FAILED" >&2
else
	echo "RLC validation passed"
fi
exit $status
