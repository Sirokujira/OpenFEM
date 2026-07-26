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

if [ "$status" -ne 0 ]; then
	echo "*** RLC validation FAILED" >&2
else
	echo "RLC validation passed"
fi
exit $status
