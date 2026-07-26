# OpenFEM

準静的 FEM による回路パラメータ (R/L/C) 抽出ソルバー (C)。
OpenFDTD (FDTD 電磁界ソルバー) と同じ構成・入力書式を踏襲した別リポジトリ。
実装は CPU (OpenMP) のみ。

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 回帰 (microstrip): Z0 / eps_eff が変更前と一致すること
mkdir -p /tmp/smoke && cp data/sample/microstrip.ofe /tmp/smoke/ && cd /tmp/smoke
$OLDPWD/bin/ofe -n 2 microstrip.ofe && grep "normal end" ofe.log

# 検証 (解析解との比較: 平行平板 1%, 抵抗 1%, 同軸 5%)
sh data/sample/rlc_check.sh bin/ofe bin/ofe_post /tmp/rlc-check
```

## 移植性の絶対規則 (OpenFDTD が Windows CI で踏んだもの)

- **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
  フラット配列を使う。
- **OpenMP のループ変数は `int`**。MSVC の OpenMP 2.0 は 64bit のループ変数を
  受け付けない。節点数が `INT_MAX` 未満であることは `setup()` で確認済みなので、
  節点ループは `const int n = (int)A->n;` として `int` で回す。
- **float\*/double\* の取り違え禁止**: 配列の実型と読み出しポインタ型の
  不一致は Windows で 0xC0000005 クラッシュ。本リポジトリは全て `double`。
- libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
- MSVC フラグは CMakeLists の既存ブロックに従う
  (`/utf-8`, `_USE_MATH_DEFINES`, `/STACK:16777216`)。
- 数学定数は `PI` / `EPS0` / `MU0` / `C0` の既存マクロを使う。

## 数値まわりの規則

- 剛性行列は CRS (27 点ステンシル)。構造格子なので格納位置は
  `crs_offset()` で探索なしに計算する。パターンを変えるときは
  `crs_alloc()` と `crs_offset()` を必ず同時に直す。
- Dirichlet 条件は行列を書き換えず、`fix[]` マスクで恒等行として扱う
  (`crs_spmv`)。解いた後に**元の行列**から反作用 Σ(Kφ) を取って
  電荷・電流にするため、行列を破壊してはいけない。
- 要素行列は 2x2x2 Gauss 積分 (直方体要素に対して厳密)。材料係数は
  `assemble.c` の `material_coef()` に集約する (mode 0..3)。
- 静磁場 (`analysis = M`) は純 Neumann 系。正味電流が 0 なので可解だが特異なので、
  1 節点だけ固定して解く。L は相互エネルギー ∫J・A から取るので定数分に依らない。
  **右辺 b は L の計算にそのまま使うので、ソルバーに渡すコピー側だけ
  固定節点を 0 にすること**。
- 渦電流 (`analysis = F`) は複素対称系 K + jωM を COCG で解く。内積に共役を
  取らない双一次形式を使うこと (共役を取ると収束しない)。境界は自然境界条件
  (磁気壁) のままにする。**Dirichlet で切ると帰路電流が境界に逃げて
  Z_loop が壊れる**。基準導体は V'=0 に固定し、ポートのみ励振して
  Y (np x np) を作り、Z_loop = inv(Y) とする。
- 複素数は実部・虚部の配列を別に持つ (C99 complex は MSVC 互換性の都合で使わない)。
- 非線形磁性体 (`bh`) は Newton-Raphson。ν の逐次代入は飽和領域で発散するので
  使わない。残差とヤコビアンで ν の評価点 (Gauss 点) を揃えること。
- 収束判定は相対残差。`solver` キーの既定は `10000 100 1e-9`。

## 機能追加の規則

- 入力キー追加は `sol/input_data.c` に、既定値は「キー省略時に従来動作と
  完全一致」になるよう初期化する (後方互換)。
- 新機能には data/sample/ の検証ケース (できれば解析解付き) と
  `data/sample/rlc_check.sh` への追加、CI スモーク (3 OS) を必ず付ける。
- `ofe_post` が `ofe.out` の新フィールドを無害に扱えることを確認する
  (フォーマットを変えたら先頭のマジック `OFEOUT01` も上げる)。
- CPU 実装にだけ追加した機能は README に対応状況を明記する。

## 設定ファイル (.claude/)

詳細な規則は `.claude/rules/` に分けてあり、該当ファイルを編集するときだけ
自動で読み込まれます (`portability.md` / `numerics.md` / `validation.md`)。
`/build` `/verify` `/add-analysis` のスラッシュコマンドと、MSVC で落ちる書き方を
編集直後に検出するフックも入っています。詳しくは `.claude/README.md`。

## CI

`.github/workflows/ci.yml`: Linux / macOS (libomp) / Windows (MSVC + Ninja)。
外部ライブラリ依存は無い。Windows の検証ステップは Git for Windows の bash で
`rlc_check.sh` を実行する。タグ `v*` push で Release にバイナリを添付。
