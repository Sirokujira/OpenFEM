---
description: 新しい解析モード (analysis = X) を追加するときの手順
argument-hint: [解析の名前と支配方程式]
---

新しい解析モード「$ARGUMENTS」を追加します。既存の C / L / R / M / F と同じ
組み立てに従ってください。

1. **定式化を先に書く** — 弱形式、係数、境界条件、抽出する量 (何から行列を作るか)
   を `sol/solve.c` の関数コメントに日本語で残す。ここが曖昧なまま実装しない。
2. `include/fem.h` に `ANALYSIS_X` ビットと結果行列・`HaveX` を追加。
3. `sol/input_data.c` の `analysis` 解釈に文字を追加。必要な前提
   (`tline` が要る / `frequency` が要る 等) は入力段でエラーにする。
4. 係数が新しいなら `sol/assemble.c` の `material_coef()` に mode を足す。
5. `sol/solve.c` に `solve_xxx()` を追加。既存の CRS・マスク・PCG/COCG を再利用する
   (`.claude/rules/numerics.md` の不変条件を守ること)。
6. `sol/setup.c` で行列を確保・解放。
7. `sol/outputRLC.c` にログ出力と `ofe.out` への書き出しを追加し、
   **マジックを 1 つ上げて** `post/post_Main.c` を追従させる。
8. `sol/utils.c` の `monitor2()` に解析名を追加。
9. `data/sample/` に**解析解と比較できる**ケースを追加し、
   `data/sample/rlc_check.sh` に比較を足す。
10. README の解析一覧・入力キー表・サンプル表・「制限・未対応」を更新。
    CLAUDE.md にも不変条件があれば追記。

最後に `/verify` を実行し、既存項目が全て通ることを確認してから報告してください。
