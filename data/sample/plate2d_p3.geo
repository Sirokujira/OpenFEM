// 断面 2 次元の平行平板線路 (plate2d.msh と同じ形状・同じ物理タグ) を gmsh で 3 次三角形に
// 切る。tline = X なので断面は y-z 面 (x = 0)。
//   導体 0 : z = -t .. 0、導体 1 : z = d .. d+t (幅 W)、間隙 z = 0 .. d は真空
// 物理タグ 1 = 真空、10 = 導体 0、11 = 導体 1 (plate2d_dc.ofe と同じ)
// 構造的 (transfinite) に切るので、-order 2 / 3 は同じ頂点の格子になる (次数の比較用):
//   gmsh plate2d_p3.geo -2 -order 3 -format msh22 -o p3raw.msh
//   gmsh p3raw.msh -0 -format msh22 -o plate2d_p3.msh
//   gmsh p3raw.msh -0 -format msh41 -o plate2d_p3_41.msh
//   gmsh plate2d_p3.geo -2 -order 2 -format msh22 -o plate2d_p3o2.msh
// (gmsh は変換のたびに節点を振り直すので、2.2 と 4.1 の一致を見る対は同じ変換で出す)
W = 1e-3; T = 0.05e-3; D = 0.2e-3;
If (!Exists(NW))
  NW = 10;	// 幅の分割
EndIf
If (!Exists(NG))
  NG = 3;	// 間隙の分割 (静磁場の過渡解析の積層板 lam2d.msh は NW 4 / NG 40 / NT 2 の 1 次 :
		//   gmsh plate2d_p3.geo -2 -setnumber NW 4 -setnumber NG 40 -setnumber NT 2 -setnumber GB 0.5 -format msh22 -o lam2d.msh。
		//   2 次の lam2d_p2.msh は NG 10 -order 2、導体の J-A の plate2d_nt8 / _p2 / _p3.msh は NT 8 -order 1 / 2 / 3 で、
		//   どちらも p3raw.msh と同じく gmsh で一度切ってから -0 で 2.2 に書き直す)
EndIf
If (!Exists(NT))
  NT = 4;	// 導体内の分割 (gmsh -setnumber NT n で変える : 収束次数の検査用)
EndIf
z[] = {-T, 0, D, D + T};
For k In {0:3}
  Point(2 * k + 1) = {0, -W / 2, z[k]};
  Point(2 * k + 2) = {0,  W / 2, z[k]};
EndFor
// 横の線 (y 方向) と縦の線 (z 方向)
For k In {0:3}
  Line(10 + k) = {2 * k + 1, 2 * k + 2};
  Transfinite Curve{10 + k} = NW + 1;
EndFor
For k In {0:2}
  Line(20 + k) = {2 * k + 1, 2 * k + 3};
  Line(30 + k) = {2 * k + 2, 2 * k + 4};
EndFor
// 導体の中は両面に向けて細かく (表皮効果)、間隙は一様
Transfinite Curve{20, 30, 22, 32} = NT + 1 Using Bump 0.5;
If (!Exists(GB))
  Transfinite Curve{21, 31} = NG + 1;
Else
  Transfinite Curve{21, 31} = NG + 1 Using Bump GB;	// 間隙も両面に向けて細かく (lam2d.msh : 面積の揃わない三角形)
EndIf
For k In {0:2}
  Curve Loop(40 + k) = {10 + k, 30 + k, -(11 + k), -(20 + k)};
  Plane Surface(50 + k) = {40 + k};
  Transfinite Surface{50 + k};
EndFor
Physical Surface(10) = {50};
Physical Surface(1) = {51};
Physical Surface(11) = {52};
