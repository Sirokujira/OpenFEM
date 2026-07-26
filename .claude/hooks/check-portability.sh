#!/bin/sh
# check-portability.sh — 編集直後に「Windows CI で実際に踏んだ」移植性違反を検出する
#
# PostToolUse (Edit/Write/MultiEdit) フックとして起動され、標準入力から
# {"tool_input": {"file_path": "..."} ...} を受け取る。
# 違反を見つけたら exit 2 で標準エラーに理由を出す (Claude に差し戻される)。
#
# 検出するもの (いずれも regex で誤検出なく判定できるものだけに絞る):
#   1. OpenMP の for-init 内でのループ変数宣言 -> MSVC error C3015
#   2. <complex.h> の使用 -> MSVC 互換性のため実部・虚部の配列を使う規則
#   3. libm を直接指定した target_link_libraries -> MATH_LIB を経由する規則

set -u

payload=$(cat)

# file_path を取り出す (python3 があれば JSON として、無ければ sed で)
if command -v python3 > /dev/null 2>&1; then
	file=$(printf '%s' "$payload" | python3 -c \
		'import json,sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
ti = d.get("tool_input") or {}
print(ti.get("file_path") or ti.get("notebook_path") or "")' 2>/dev/null)
else
	file=$(printf '%s' "$payload" | sed -n 's/.*"file_path"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p')
fi

[ -n "${file:-}" ] || exit 0
[ -f "$file" ] || exit 0

status=0
report() {
	echo "$1" >&2
	status=2
}

case "$file" in
*.c|*.h)
	# 1. #pragma omp parallel for が支配する for 文でループ変数を宣言すると MSVC が落ちる。
	#    pragma の後で最初に現れる for 文だけが対象 (入れ子の内側は宣言してよい)。
	hits=$(awk '
		/#pragma[ \t]+omp[ \t]+parallel[ \t]+for/ { armed = 1; next }
		armed && /for[ \t]*\(/ {
			if ($0 ~ /for[ \t]*\([ \t]*(int|long|size_t|int64_t|unsigned)[ \t]+[A-Za-z_]/) {
				printf "%s:%d:%s\n", FILENAME, FNR, $0
			}
			armed = 0
		}
	' "$file")
	if [ -n "$hits" ]; then
		report "*** MSVC の OpenMP (C モード) は for-init 内のループ変数宣言を受け付けない (error C3015)。
    ループ変数は #pragma の前で宣言すること (.claude/rules/portability.md 参照):
$hits"
	fi

	# 2. C99 complex は MSVC 互換性のため使わない
	if grep -n '#include[ \t]*<complex\.h>' "$file" > /dev/null 2>&1; then
		report "*** <complex.h> は使わない。複素数は実部・虚部の double 配列で持つこと
    (sol/solver_cocg.c を参照)。"
	fi
	;;
CMakeLists.txt|*/CMakeLists.txt)
	# 3. libm の直接指定 (Windows に m.lib は無い)
	if grep -n 'target_link_libraries[^)]*[^A-Za-z_}]m[)) \t]' "$file" > /dev/null 2>&1; then
		report "*** libm は MATH_LIB 変数を経由してリンクすること (Windows に m.lib は無い)。"
	fi
	;;
esac

exit $status
