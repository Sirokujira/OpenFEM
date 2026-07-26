---
paths:
  - "data/sample/**"
  - ".github/workflows/**"
  - "post/**"
description: 検証ケース・CI・入出力互換の規則
---

# 検証と入出力の規則

## 新機能には必ず検証ケースを付ける

- `data/sample/` に**解析解と比較できる**ケースを追加し、
  `data/sample/rlc_check.sh` に比較を足して CI (3 OS) で回します。
- 許容誤差は「測った値」ではなく**物理的な根拠**で決めること
  (1 次元厳密 → 1%、円形導体の階段近似 → 5% など)。ケースのコメントに
  期待値の導出を書きます。
- 期待値は awk で計算できる形に落とすか、`.ofe` のコメントに導出を残して
  スクリプトへ定数で埋め込みます。

## 入力キーの追加

- `sol/input_data.c` に追加し、既定値は「**キー省略時に従来動作と完全一致**」
  になるよう初期化します (後方互換)。
- 材料番号を参照するキー (`mur` / `tand` / `debye` など) は `material` 行より
  後に書く必要がある旨を README に明記します。

## 出力の互換

- `ofe.out` のフォーマットを変えたら先頭のマジック (`OFEOUT0n`) を上げ、
  `post/post_Main.c` を必ず追従させます。
- `ofe_post` が新しいフィールドを無害に扱えることを確認すること。

## CI

`.github/workflows/ci.yml` は Linux / macOS (libomp) / Windows (MSVC + Ninja) の 3 OS。
外部ライブラリ依存はありません。Windows の検証ステップは Git for Windows の
bash で `rlc_check.sh` を実行します。タグ `v*` push で Release にバイナリを添付。

CI を落としたまま次の作業に進まないこと。Windows だけ落ちた場合は
`.claude/rules/portability.md` を先に読み直します。
