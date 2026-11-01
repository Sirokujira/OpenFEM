// 断面 2 次元の 3 次三角形 (tri10) : 入力エラーになること (rlc_check.sh) の検証用
// gmsh tri10_2d.geo -2 -order 3 -format msh22 -o tri10_2d.msh
Point(1) = {0, 0, 0, 0.5e-3}; Point(2) = {1e-3, 0, 0, 0.5e-3};
Point(3) = {1e-3, 1e-3, 0, 0.5e-3}; Point(4) = {0, 1e-3, 0, 0.5e-3};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4}; Plane Surface(1) = {1};
Physical Surface(1) = {1};
