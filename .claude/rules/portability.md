---
paths:
  - "**/*.c"
  - "**/*.h"
  - "CMakeLists.txt"
description: Linux / macOS / Windows (MSVC) の 3 OS を通すための絶対規則
---

# 移植性の絶対規則

<important if="C / ヘッダー / CMakeLists を編集するとき">
以下はいずれも **CI の Windows ジョブで実際に落ちた** ものです。ローカル (gcc) では
通ってしまうため、書いた時点で守ること。`.claude/hooks/check-portability.sh` が
1 と 2 を編集直後に検出します。
</important>

## 1. OpenMP のループ変数は pragma の外で宣言する

MSVC の OpenMP は C モードで for-init 内の宣言を受け付けません
(`error C3015: initialization in OpenMP 'for' statement has improper form`)。

```c
/* NG */                        /* OK */
#pragma omp parallel for        int i;
for (int i = 0; i < n; i++)     #pragma omp parallel for
                                for (i = 0; i < n; i++)
```

入れ子の内側 (pragma が支配しないループ) は `for (int j = ...)` のままで構いません。

また MSVC の OpenMP は 64bit のループ変数を受け付けないため、節点ループは
`const int n = (int)A->n;` として `int` で回します
(節点数が `INT_MAX` 未満であることは `setup()` で確認済み)。

## 2. C99 complex は使わない

複素数は実部・虚部の `double` 配列を別々に持ちます (`sol/solver_cocg.c` 参照)。
`<complex.h>` を include しないこと。

## 3. C99 VLA 禁止

MSVC の C2057 / C2466。`malloc` + 明示インデックスのフラット配列を使います。
`double ke[8][8]` のような**定数長**の配列は問題ありません。

## 4. libm は MATH_LIB を経由する

Windows に `m.lib` は無いので、CMakeLists の `MATH_LIB` 変数を通してリンクします
(MSVC では空文字列)。継続行の `m` も置換対象になるので注意。

## 5. float\* / double\* の取り違え禁止

配列の実型と読み出しポインタ型の不一致は Windows で 0xC0000005 クラッシュになります
(glibc は偶然耐えてしまう)。本リポジトリは**全て `double`** で統一しています。

## 6. MSVC フラグは CMakeLists の既存ブロックに従う

`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`, `/STACK:16777216`。

## 7. 数学定数は既存マクロを使う

`PI` / `EPS0` / `MU0` / `C0` / `ETA0` (`include/fem.h`)。
