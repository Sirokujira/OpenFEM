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

いずれも節点スカラー量 u (電位 φ またはベクトルポテンシャル成分 Az) に対する
準静的な弱形式

```
∫ c ∇w・∇u dV = ∫ w f dV
```

を、構造格子上の 8 節点 3 次補間 (trilinear) 6 面体要素で離散化します
(2x2x2 Gauss 積分。直方体要素なので厳密積分)。係数 c は方向毎に持てるので
対角テンソルの異方性材料 (`anisoeps` / `anisomur`) を扱えます。

| 解析 | 方程式 / 係数 | 得られる量 |
|---|---|---|
| `C` | ∇・(ε0 εr ∇φ) = 0 | 短絡容量行列 (Maxwell 行列) [F] または [F/m] |
| `L` | ∇・(ε0 ∇φ) = 0 (真空) | L = μ0 ε0 inv(C0) : TEM 近似 (外部インダクタンスのみ) [H/m] |
| `R` | ∇・((σ + ω ε0 εr tanδ) ∇φ) = 0 | 並列コンダクタンス行列 G と R = inv(G)。導電損 + 誘電損 |
| `M` | ∇・(ν ∇Az) = −Jz | DC インダクタンス行列 [H/m]。**内部インダクタンス込み、μr 対応**。B-H 曲線があれば**非線形** |
| `F` | ∇・(∇Az/(μ0 μr)) − jωσAz + σV' = 0 | 周波数 f での直列 R(f), L(f) [Ω/m, H/m]。**表皮効果・近接効果込み** |

`M` (静磁場) はベクトルポテンシャル Az の 2 次元断面定式化です。ポート k に +I、
基準導体に −I を一様電流密度で流して解き、相互エネルギー
∫J^(j)・A^(k) dV = L_kj I_k I_j から L 行列を得ます。導体内部の電流分布を
含むので内部インダクタンスが自動的に入ります (正味電流が 0 なので純 Neumann 系は
可解。1 節点を固定して解いた解を使い、結果は定数分に依りません)。

材料に B-H 曲線 (`bh`) を与えると `M` は **非線形静磁場解析** になります。
ν = H(|B|)/|B| を Gauss 点毎に評価し、**Newton-Raphson 法**

```
R(A) = K(ν(|∇A|)) A − b = 0
J_ij = Σ_g w [ ν (∇N_i・∇N_j) + 2 (dν/dB²) (∇N_i・∇A)(∇N_j・∇A) ] detJ
```

で解きます (H(B) が単調なら J は正定値なので Jacobi-PCG が使えます)。
得られる L は与えた電流での**磁束鎖交 (割線) インダクタンス**で、電流を上げると
飽和して下がります。重ね合わせが成り立たないので単一ポート専用です。
ν の逐次代入は飽和領域で振動するため使っていません。

`F` (時間調和渦電流) は複素 Az の 2 次元断面定式化です。基準導体を V'=0 に
固定してポート j を V'=1 [V/m] で励振し、各ポートの電流
I_k = ∫_k σ(V'_k − jωAz) dS から単位長あたりのアドミタンス行列 Y を作り、
Z_loop = inv(Y) から R(f) = Re Z、L(f) = Im Z/ω を得ます。導体内部の電流分布を
解くので**表皮効果・近接効果**が入ります (複素対称行列なので COCG 法で解く)。
境界条件は自然境界条件 (磁気壁) で、これは「帰路電流はモデル内の導体を流れる」
という伝送線路の前提そのものです。格子が表皮深さを刻めていない場合は警告します。

導体に `conductorsigma` を与えると、断面積から **DC 直列抵抗**
Rs = 1/(σA) [Ω/m] を求めます (帰路は基準導体が共有するものとして
Rs[k][j] = R0 + (k=j のとき Rk))。これで RLGC が揃います。

導体 (電極) 上の節点を Dirichlet 境界とし、ポート k のみを `voltage` [V]、
他の導体と基準導体 (id=0) を 0 [V] として解きます。各導体上の反作用
Σ (K φ)_i から電荷 Q_j (静電界) / 電流 I_j (電流界) を求め、
`C[k][j] = Q_j / V`、`G[k][j] = I_j / V` として行列を組み立てます。

外側境界はいずれの解析も自然境界条件 (∂u/∂n = 0) です。静電界では電気力線が
境界に平行に、静磁場では磁束線が境界に垂直になります。閉じたシールド構造では
これが厳密な条件になり、開放問題では解析領域を対象より十分広く取ってください。

線形ソルバーは対角スケーリング前処理付き共役勾配法 (実対称は Jacobi-PCG、
渦電流の複素対称系は COCG、いずれも OpenMP 並列) です。剛性行列は
CRS (27 点ステンシル) で保持します。

## 処理部の構成

| ディレクトリ | 役割 |
|---|---|
| `src/sol_Main.c` | ソルバー `ofe` のエントリ。入力読込 → `setup()` → `solve()` → 出力 |
| `post/post_Main.c` | ポストプロセッサ `ofe_post`。`ofe.out` → `rlc.csv` / `ofe_circuit.sp` |
| `sol/input_data.c` | `.ofe` 入力の解釈 |
| `sol/setup.c` | 格子展開、セル材料・導体セル/節点の割り当て、導体断面積と DC 直列抵抗 |
| `sol/ingeometry.c` | 形状内外判定 (形状番号は OpenFDTD の `geometry` と共通) |
| `sol/crs.c` | CRS 疎行列 (構造格子なので格納位置を探索なしで計算) |
| `sol/assemble.c` | 6 面体要素の要素行列と全体行列の作成 |
| `sol/solver_cg.c` | Jacobi-PCG (実対称) |
| `sol/solver_cocg.c` | COCG (複素対称、渦電流解析用) |
| `sol/solve.c` | 各解析の駆動、回路パラメータの抽出、静磁場解析 |
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
- `rlc.csv` — 行列 (C, L, Ldc, G, R, Rs) と伝送線路定数の一覧
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
| `mur` | `<material_id> <mur>` | 1 | 比透磁率 (静磁場解析)。`material` の後に書く |
| `bh` | `<material_id> <H> <B>` | なし | B-H 曲線の点 (複数行)。H, B とも正で単調増加、原点は書かない。与えると `M` が非線形解析になる |
| `nlsolver` | `<maxiter> <tol> <damping>` | `50 1e-5 1.0` | 非線形 (B-H) Newton 反復の設定 |
| `tand` | `<material_id> <tand>` | 0 | 誘電正接。`frequency` と併用して誘電損 G を求める |
| `debye` | `<material_id> <eps_inf> (<deps> <tau>)...` | なし | Debye 分散材料 (多極可)。`frequency` の値で εr と tanδ に展開する |
| `lorentz` | `<material_id> <eps_inf> (<deps> <f0> <delta>)...` | なし | Lorentz 分散材料 (多極可)。`debye` と混在させると極が足し合わされる |
| `anisoeps` | `<material_id> <ex> <ey> <ez>` | 等方性 | 異方性の比誘電率 (対角テンソル) |
| `anisomur` | `<material_id> <mx> <my> <mz>` | 等方性 | 異方性の比透磁率 (対角テンソル) |
| `conductorsigma` | `<conductor_id> <sigma>` | なし | 導体の導電率。DC 直列抵抗 Rs の計算に使う |
| `frequency` | 実数 [Hz] | 0 | 誘電損 (tanδ)・分散材料・渦電流 (`F`) を評価する周波数。0 なら tanδ は無視 |
| `current` | 実数 [A] | 1 | 静磁場解析の励振電流 |
| `analysis` | `C` / `L` / `R` / `M` / `F` の組 | `C` | 実行する解析 |
| `tline` | `X` / `Y` / `Z` | なし | 伝送線路軸。指定すると単位長あたりの値で出力する (`L`, `M`, `F` は必須) |
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
| `plate_line_dc.ofe` | 平行平板線路の静磁場 | L'dc = μ0(d + 2t/3)/W と一致 (誤差 0.01%)、Rs' も厳密 |
| `coax_loss.ofe` | 損失同軸 (C/L/G/L_dc/Rs) | L'dc 誤差 0.22%、Rs' 誤差 0.53%、G' = ωC'tanδ を厳密に再現 |
| `plate_line_ac.ofe` | 平行平板線路の渦電流 | 1 次元厳密解 Z = 2γcoth(γt)/(σW) + jωμ0d/W と比較 (1kHz: R 誤差 0.00% / 10MHz: R 誤差 0.12%) |
| `plate_line_bh.ofe` | 非線形磁性体を挟んだ平行平板線路 | L(I) = B(I/W)d/I + 2μ0t/(3W) と 4 電流で比較 (いずれも誤差 0.00%、Newton 2〜3 回) |
| `dispersive_plate.ofe` | 多極分散 (Debye + Lorentz) | 1GHz での C, G が展開式と厳密一致 (誤差 0.00%) |
| `aniso_plate.ofe` | 異方性誘電体 (εx,εy,εz = 10,5,2) | C = ε0 εz A/d と厳密一致 (軸の対応の検査) |

検証スクリプト:

```bash
sh data/sample/rlc_check.sh bin/ofe bin/ofe_post /tmp/rlc-check
```

## 制限・未対応

- 実装は CPU (OpenMP) のみ。**CUDA / MPI は未対応** (OpenFDTD にある GPU・
  マルチプロセス版に相当するものは無い)。
- 格子は直交構造格子。斜め・曲面の導体は階段近似になる (同軸で ~2%)。
  非構造格子・適合格子は未対応。
- `L` (TEM) は外部インダクタンスのみ。内部インダクタンスと μr が要るときは
  `M` (静磁場) を使う。`M` は 2 次元断面の定式化なので `tline` が要る。
- 表皮効果は `F` (渦電流) で扱う。`Rs` と `M` の L は DC 値 (電流一様) なので、
  高周波では `F` を使うこと。格子が表皮深さを刻めていないと R(f) を
  過小評価するため、その場合は警告を出す。
- 分散材料は Debye / Lorentz の多極 (最大 8 極)。Drude・Cole-Cole は未対応。
- 異方性は**対角テンソル**のみ (主軸が格子軸に一致する場合)。非対角成分
  (任意方向の主軸) と、異方性 + 非線形 (B-H) の併用は未対応。
- 非線形磁性体は**単調な B-H 曲線のみ**。ヒステリシス (履歴・残留磁束)、
  異方性、非線形と渦電流 (`F`) の同時解析は未対応。
- `M` / `F` は 2 次元断面 (伝送線路) 専用。3 次元の渦電流 (A-φ 定式化)、
  変位電流を含む全波解析は未対応 (全波は OpenFDTD 側の担当)。

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
