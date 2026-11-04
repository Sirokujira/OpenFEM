// 円柱の導体 (半径 R = 0.5mm、長さ L = 0.5mm、z 軸) の四面体格子 (熱連成の曲がった面の検証)
//   gmsh cyl3.geo -3 -order 2 -format msh22 -o cyl_p2.msh   (-order 1 / 3 も同じ geo)
// 物理タグ 1 = 体積、10 = z = 0 の端面 (電極 0)、11 = z = L の端面 (電極 1)、20 = 側面
SetFactory("OpenCASCADE");
R = 0.5e-3; L = 0.5e-3;
If (!Exists(h))
  h = 0.25e-3;
EndIf
Cylinder(1) = {0, 0, 0, 0, 0, L, R};
Physical Volume(1) = {1};
Physical Surface(10) = {3};
Physical Surface(11) = {2};
Physical Surface(20) = {1};
Mesh.MeshSizeMin = h; Mesh.MeshSizeMax = h;
Mesh.HighOrderOptimize = 0;
