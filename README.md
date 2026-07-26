# OpenFEM

準静的 **有限要素法 (FEM)** による、電気回路パラメータ (R / L / C) 抽出ソルバー。

配線・基板パターン・伝送線路などの 3 次元形状から、静電界・定常電流界を
FEM で解いて **容量行列 C・インダクタンス行列 L・抵抗/コンダクタンス行列 R, G**
を求め、SPICE の等価回路 (サブサーキット) として書き出します。

[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) と同じ構成 (C99 + CMake +
OpenMP、`mesh`/`material`/`geometry` の入力書式、`ofd.log` 相当のログ、
ポストプロセッサ分離) を踏襲した別リポジトリです。FDTD が「時間領域の
電磁界」を解くのに対し、こちらは「回路として見たときの定数」を求めます。

## 何を解くか

いずれも節点スカラーポテンシャル φ に対する準静的な弱形式

```
∫ c ∇w・∇φ dV = 0
```

を、構造格子上の 8 節点 3 次補間 (trilinear) 6 面体要素で離散化します
(2x2x2 Gauss 積分。直方体要素なので厳密積分)。

| 解析 | 係数 c | 得られる量 |
|---|---|---|
| `C` | ε0 εr | 短絡容量行列 (Maxwell 行列) [F] または [F/m] |
| `L` | ε0 (真空) | L = μ0 ε0 inv(C0) : TEM 近似のインダクタンス行列 [H/m] |
| `R` | σ | コンダクタンス行列 G と R = inv(G) |

導体 (電極) 上の節点を Dirichlet 境界とし、ポート k のみを `voltage` [V]、
他の導体と基準導体 (id=0) を 0 [V] として解きます。各導体上の反作用
Σ (K φ)_i から電荷 Q_j (静電界) / 電流 I_j (電流界) を求め、
`C[k][j] = Q_j / V`、`G[k][j] = I_j / V` として行列を組み立てます。

外側境界は自然境界条件 (∂φ/∂n = 0、磁気壁) です。開放問題では
解析領域を対象より十分広く取ってください。

線形ソルバーは対角スケーリング前処理付き共役勾配法 (Jacobi-PCG、OpenMP 並列)
です。剛性行列は CRS (27 点ステンシル) で保持します。

## 処理部の構成

| ディレクトリ | 役割 |
|---|---|
| `src/sol_Main.c` | ソルバー `ofe` のエントリ。入力読込 → `setup()` → `solve()` → 出力 |
| `post/post_Main.c` | ポストプロセッサ `ofe_post`。`ofe.out` → `rlc.csv` / `ofe_circuit.sp` |
| `sol/input_data.c` | `.ofe` 入力の解釈 |
| `sol/setup.c` | 格子展開、セル材料・導体節点の割り当て |
| `sol/ingeometry.c` | 形状内外判定 (形状番号は OpenFDTD の `geometry` と共通) |
| `sol/crs.c` | CRS 疎行列 (構造格子なので格納位置を探索なしで計算) |
| `sol/assemble.c` | 6 面体要素の要素行列と全体行列の作成 |
| `sol/solver_cg.c` | Jacobi-PCG |
| `sol/solve.c` | 各解析の駆動と回路パラメータの抽出 |
| `sol/matutil.c` | 小行列の逆行列 (Gauss-Jordan) |
| `sol/outputRLC.c` | ログ出力と `ofe.out` の書き出し |
| `include/` | 共有ヘッダ (`fem.h`, `fem_prototype.h`) |

## ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
# -> bin/ofe, bin/ofe_post
```

依存は C コンパイラと CMake のみ (OpenMP は任意、`-DWITH_OPENMP=OFF` で無効化)。
Linux / macOS (AppleClang + libomp) / Windows (MSVC) で CI を回しています。

## 使い方

```bash
bin/ofe -n 4 data/sample/microstrip.ofe   # 解析 (-n : スレッド数)
bin/ofe_post                              # ofe.out -> rlc.csv, ofe_circuit.sp
```

出力:

- `ofe.log` — 実行ログ。収束履歴と回路パラメータ。正常終了で `=== normal end ===`
- `ofe.out` — ポスト処理用バイナリ
- `rlc.csv` — 行列 (C, L, G, R) と伝送線路定数の一覧
- `ofe_circuit.sp` — SPICE サブサーキット (`.SUBCKT OFE_CIRCUIT ...`)

## 入力ファイル (.ofe)

1 行 1 キーの `key = value ...` 形式。`#` 以降はコメント、`end` で終端。
先頭行はプログラム名とバージョン。`mesh` / `material` / `geometry` の書式は
OpenFDTD の `.ofd` と共通なので、既存モデルの記述をそのまま流用できます。

```
OpenFEM 1 0
title = microstrip line
xmesh = -4e-3 30 -1e-3 80 1e-3 30 4e-3
ymesh = 0 16 0.4e-3 2 0.435e-3 60 4e-3
zmesh = 0 1 1e-4
material = 4.4 0
geometry = 2 1 -4e-3 4e-3 0 0.4e-3 0 1e-4
conductor = 0 1 -4e-3     4e-3     0      0        0 1e-4
conductor = 1 1 -0.375e-3 0.375e-3 0.4e-3 0.435e-3 0 1e-4
analysis = C L
tline = Z
end
```

| キー | 書式 | 既定値 | 説明 |
|---|---|---|---|
| `title` | 文字列 | 空 | 表題 |
| `xmesh` / `ymesh` / `zmesh` | `x0 n1 x1 n2 x2 ...` | (必須) | 領域境界と分割数 |
| `material` | `<epsr> <sigma>` | — | 材料。番号は 2 から順に振られる (0=真空, 1=PEC 予約)。`.ofd` の 4/5 値形式も受け付け、透磁率は無視する |
| `geometry` | `<material_id> <shape> <g0>..<g5>` | — | 材料を割り当てる形状 (セル中心で判定、後の行が優先) |
| `conductor` | `<id> <shape> <g0>..<g5>` | (必須) | 電極。`id=0` が基準導体 (グランド)、`id>=1` がポート。`id=-1` は取り消し (くり抜き)。節点で判定するので厚さ 0 の指定は「面」として扱える |
| `analysis` | `C` / `L` / `R` の組 | `C` | 実行する解析 |
| `tline` | `X` / `Y` / `Z` | なし | 伝送線路軸。指定すると単位長あたりの値で出力する (`L` は必須) |
| `linelength` | 実数 [m] | 解析領域長 | 等価回路を作るときの線路長 |
| `nsection` | 整数 | 1 | 等価回路 (梯子) の段数 |
| `voltage` | 実数 [V] | 1 | 励振電圧 |
| `solver` | `<maxiter> <nout> <converg>` | `10000 100 1e-9` | 反復解法の最大回数・出力間隔・収束判定 (相対残差) |

`shape` は OpenFDTD と共通で、`1`:直方体、`2`:楕円体、`11`/`12`/`13`:X/Y/Z 円柱
(いずれも引数は外接直方体 `xmin xmax ymin ymax zmin zmax`)。

## サンプルと検証

`data/sample/` :

| ファイル | 内容 | 解析解との比較 |
|---|---|---|
| `parallel_plate.ofe` | 平行平板コンデンサ | C = ε0 εr A/d と一致 (誤差 < 0.01%) |
| `resistor_bar.ofe` | 直方体抵抗 | R = d/(σA) と一致 (誤差 < 0.01%) |
| `coax.ofe` | 同軸線路 (単位長) | C', L', Z0 とも誤差 ~2% (円形導体の階段近似による) |
| `microstrip.ofe` | マイクロストリップ線路 | Z0 = 49.3 Ω, ε_eff = 3.30 (Hammerstad の目安と整合) |
| `coupled_microstrip.ofe` | 結合線路 (2 ポート) | C/L 行列と結合係数、4 段梯子の SPICE 出力 |

検証スクリプト:

```bash
sh data/sample/rlc_check.sh bin/ofe bin/ofe_post /tmp/rlc-check
```

## 制限・未対応

- 実装は CPU (OpenMP) のみ。**CUDA / MPI は未対応** (OpenFDTD にある GPU・
  マルチプロセス版に相当するものは無い)。
- 格子は直交構造格子。斜め・曲面の導体は階段近似になる (同軸で ~2%)。
- `L` は導体を完全導体とみなす **TEM 近似** (`L = μ0 ε0 inv(C_vacuum)`)。
  磁性体 (μr ≠ 1)、内部インダクタンス、表皮効果は考慮しない。
- 導体損 (表皮効果による直列抵抗) は求めない。`R` 解析は電極の置き方に
  応じて「直列 DC 抵抗」または「並列コンダクタンス」のいずれかになる。
- 周波数依存材料、渦電流 (A-φ 定式化)、静磁場の直接解析は未対応。

## OpenFDTD との関係

- 入力書式 (`mesh` / `material` / `geometry` の書き方、形状番号)、
  ビルド構成、移植性の規則 (C99 VLA 禁止・`MATH_LIB` 経由の libm リンク等) を
  OpenFDTD に合わせています。`sol/ingeometry.c` の形状定義は OpenFDTD の
  同名ファイルに由来します。
- 一方で解法は独立で、FDTD のコードは共有していません。

## ライセンス

MIT License (`LICENSE` を参照)。

`sol/ingeometry.c` の形状定義 (形状番号と引数の意味) は OpenFDTD の
同名ファイルに由来します。上流 OpenFDTD (e-em.co.jp) の配布条件も
併せて確認してください。
