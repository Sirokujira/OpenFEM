// gmsh coax3.geo -setnumber h 0.4e-3 -3 -order 3 -format msh22 -o coax_p3.msh (h = 0.2e-3 で coax_p3f.msh)
// 同軸 (内径 a = 0.5mm, 外径 b = 1.5mm, 長さ 0.1mm) の四面体格子 (高次は -order で上げる)
SetFactory("OpenCASCADE");
a = 0.5e-3; b = 1.5e-3; lz = 0.1e-3;
Disk(1) = {0, 0, 0, b};
Disk(2) = {0, 0, 0, a};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
out[] = Extrude {0, 0, lz} { Surface{3}; Layers{1}; };
Physical Volume(1) = {out[1]};
Physical Surface(10) = {out[2]};
Physical Surface(11) = {out[3]};
Mesh.MeshSizeMin = h; Mesh.MeshSizeMax = h;
Mesh.HighOrderOptimize = 0;
