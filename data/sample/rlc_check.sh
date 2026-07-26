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
#   aniso_plate    : 異方性誘電体の C = eps0 εz A/d      (厳密、許容 0.1%)
#   aniso_rot      : 非対角テンソル (z 軸 30 度回転) でも同値 (厳密、許容 0.1%)
#   plate_line_bh_aniso : 軸毎 B-H。B は x のみなので X 曲線だけが効くこと
#   plate_line_ja  : Jiles-Atherton の履歴ループを、H 掃引で独立に積分した
#                    ODE 解と比較する (11 点、許容 0.5%)
#   box_tet        : 非構造格子 (四面体) の平行平板 C     (厳密、許容 0.1%)
#   coax_tet       : 円形境界に適合した四面体格子の同軸 C', L' (許容 1%)
#   edge_test      : Whitney 辺要素 (1 次 Nedelec) の自己検証 (機械精度)
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
# 勾配の零空間・一様場の質量・回転場の回転回転・対称性を閉形式と比較する
echo "[edge_test] Nedelec edge element self test"
run_case edge_test
grad=$(awk '/gradient null space/ { print $NF }' "$WORK/ofe.log")
mass=$(awk '/uniform field mass/  { print $NF }' "$WORK/ofe.log")
curl0=$(awk '/uniform field curl/ { print $NF }' "$WORK/ofe.log")
rot=$(awk '/rotational field/     { print $NF }' "$WORK/ofe.log")
sym=$(awk '/symmetry  / { print $NF }' "$WORK/ofe.log")
for pair in "gradient-null:$grad" "mass:$mass" "uniform-curl:$curl0" "rotational:$rot" "symmetry:$sym"; do
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
