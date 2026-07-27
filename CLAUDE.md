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
# 相対パスは script 内で cd するので絶対パスで渡すこと
sh data/sample/rlc_check.sh "$PWD/bin/ofe" "$PWD/bin/ofe_post" /tmp/rlc-check
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
- 分散材料 (`debye`/`lorentz`/`drude`/`colecole`) は時間因子 e^{jωt}、ε = ε' − jε''。
  `er` に ε'、`ei` に **+ε''** (損失が正) を足す。極を足すときはこの符号規約に従う。
  Cole-Cole は α = 0 で Debye に厳密一致するので、その恒等式を検証に使う。
- 温度依存は σ の読み出し**手前で一度だけ**掛ける (`input_data.c` の末尾)。σ の
  読み出しは `Material[].sigma` と `CondSigma[]` の 2 系統あり、下流で個別に
  補正すると必ず漏れる。
- 3 次元渦電流 (`analysis = A`) は A (辺) と φ (節点) の連成。**節点側の式を
  jω で割ると複素対称**になるので COCG がそのまま使える。この形にしないと
  (1,2) = TG と (2,1) = jωGᵀT が食い違って対称でなくなる。
- A-φ 系は (A, φ) → (A + Gψ, φ − jωψ) で不変なので特異。右辺を Dirichlet の
  持ち上げ b = −A x_D で作れば b は必ず値域に入るので COCG は収束し、
  Z はゲージ不変なので**ゲージ固定しない方が速い** (実測で反復 1/6、Z は 7 桁一致)。
  A・φ 自体が要るときだけ `gauge = 1` (tree-cotree)。その場合 awall の辺を
  優先して木に入れないと境界条件と端子電圧が壊れる。
- **1 次 Nedelec 要素の誤差は要素の最大寸法で決まる** (場が変化する方向の
  刻みではない)。場が一様な方向を粗くすると誤差が出る (実測で 13%)。
  辺要素を使う解析では要素を等方的にとること。

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
