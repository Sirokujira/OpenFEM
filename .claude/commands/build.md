---
description: OpenFEM を Release でビルドし、警告ゼロを確認する
---

`ofe` / `ofe_post` をビルドします。

1. `cmake -B build -DCMAKE_BUILD_TYPE=Release` (既に構成済みならスキップ)
2. `cmake --build build -j"$(nproc)"`
3. 警告が出ていないか確認する。出ていたら直してから次に進む。

厳しめの確認をするときは OpenMP を切って警告を全部出す:

```bash
cmake -B /tmp/bw -DCMAKE_BUILD_TYPE=Release -DWITH_OPENMP=OFF \
  -DCMAKE_C_FLAGS="-Wall -Wextra -Wpedantic"
cmake --build /tmp/bw -j4
```

ビルドが通ったら報告は簡潔に。バイナリは `bin/` に出ます。
