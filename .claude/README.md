# .claude — Claude Code 設定

このディレクトリは [Claude Code](https://claude.com/claude-code) 用の
プロジェクト設定です。git に入れてチームで共有します。

構成の考え方は次の 2 つを参考にしています:

- [shanraisshan/claude-code-best-practice](https://github.com/shanraisshan/claude-code-best-practice)
  — CLAUDE.md は 200 行未満に保ち、大きくなったら `rules/` に分割する。
  決まりきった挙動は指示文ではなく `settings.json` で決定論的に縛る。
  日常的に繰り返す作業はスラッシュコマンドにする。
- [affaan-m/ECC (Everything Claude Code)](https://github.com/affaan-m/ECC)
  — `agents/` `commands/` `rules/` `hooks/` の分け方と、
  エージェントの frontmatter (`name` / `description` / `tools` / `model`) 形式。

## 中身

| パス | 役割 |
|---|---|
| `settings.json` | 権限 (許可/拒否) と PostToolUse フックの登録 |
| `rules/portability.md` | C/ヘッダー/CMakeLists 編集時に自動で読まれる移植性の絶対規則 |
| `rules/numerics.md` | `sol/` 編集時に読まれる行列・境界条件・解法の不変条件 |
| `rules/validation.md` | サンプル/CI/入出力互換の規則 |
| `commands/build.md` | `/build` — ビルドと警告確認 |
| `commands/verify.md` | `/verify` — 解析解との比較検証 |
| `commands/add-analysis.md` | `/add-analysis` — 新しい解析モードの追加手順 |
| `agents/fem-reviewer.md` | 移植性 + 数値の観点で差分をレビューするサブエージェント |
| `hooks/check-portability.sh` | 編集直後に MSVC で落ちる書き方を検出する (exit 2 で差し戻す) |

`rules/*.md` は frontmatter の `paths:` にマッチするファイルを触るときだけ
読み込まれます (常時読み込みではないので CLAUDE.md を膨らませずに済みます)。

## フック

`Edit` / `Write` / `MultiEdit` の直後に `hooks/check-portability.sh` が走り、
**Windows CI で実際に落ちた**書き方を検出します:

1. `#pragma omp parallel for` が支配する for 文でのループ変数宣言 (MSVC `error C3015`)
2. `<complex.h>` の使用
3. `target_link_libraries` での libm 直接指定

いずれも正規表現で誤検出なく判定できるものだけに絞ってあります。
検出したら exit 2 で差し戻すので、そのまま直してください。

動作確認:

```bash
echo '{"tool_input":{"file_path":"sol/crs.c"}}' | sh .claude/hooks/check-portability.sh
echo "exit=$?"   # 0 なら問題なし
```

## CLAUDE.md との使い分け

- `CLAUDE.md` … 毎回必要な最小限 (ビルド/テストの叩き方、機能追加の手順、CI の構成)
- `.claude/rules/` … 該当ファイルを触るときだけ必要な詳細
- `.claude/settings.json` … 指示ではなく仕組みで縛れるもの (権限・フック)
