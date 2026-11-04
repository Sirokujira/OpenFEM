// 3 次四面体 (tet20 / tri10) の検証用の箱 : gmsh box_p3.geo -3 -order 3 -format msh22 -o box_p3.msh
lx = 1e-3; ly = 1e-3; lz = 0.2e-3;
nx = 3; ny = 3; nz = 2;
Point(1) = {0,0,0,1}; Point(2) = {lx,0,0,1}; Point(3) = {lx,ly,0,1}; Point(4) = {0,ly,0,1};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
Transfinite Curve{1,3} = nx+1; Transfinite Curve{2,4} = ny+1;
Transfinite Surface{1};
out[] = Extrude {0,0,lz} { Surface{1}; Layers{nz}; };
Physical Volume(1) = {out[1]};
Physical Surface(10) = {1};
Physical Surface(11) = {out[0]};
