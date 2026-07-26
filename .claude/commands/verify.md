---
description: 解析解との比較検証 (rlc_check.sh) を実行して結果を報告する
---

回路パラメータが解析解と一致するか検証します。**コードを変更したら必ず実行**すること。

```bash
cmake --build build -j"$(nproc)"
sh data/sample/rlc_check.sh "$PWD/bin/ofe" "$PWD/bin/ofe_post" /tmp/rlc-check
```

確認すること:

- 全項目が `OK` であること (1 つでも `NG` なら原因を特定してから次へ)
- 誤差の傾向が変わっていないこと。厳密解のケース (平行平板 / 抵抗 /
  plate_line_dc / plate_line_ac の 1kHz) は誤差 0.0x% を維持しているはず。
  ここが悪化したら離散化ではなく実装の誤りを疑う。

報告は「どの項目が何 % ずれたか」を表で簡潔に。全部 OK なら 1 行で十分です。
