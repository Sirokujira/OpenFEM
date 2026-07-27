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
#   input lint     : 選んだ解析が読まないキーを警告すること (5 つの罠) と、
#                    正しい 17 ケースで警告が 1 件も出ないこと
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
         temp_resistor aniso_plate box_tet coax_tet edge_test bar_eddy; do
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
