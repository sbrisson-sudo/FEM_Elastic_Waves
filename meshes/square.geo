L = 100.0;
dx = 20.0;

Point(1) = {0.0, 0.0, 0.0, dx};
Point(2) = {  L, 0.0, 0.0, dx};
Point(3) = {  L,   L, 0.0, dx};
Point(4) = {0.0,   L, 0.0, dx};

Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};

Curve Loop(21) = {11, 12, 13, 14};

Plane Surface(111) = {21};

Transfinite Surface {111};
Recombine Surface {111};

Mesh.ElementOrder = 1;


Physical Curve("bottom",    31) = {11};
Physical Curve("right",     32) = {12};
Physical Curve("top",       33) = {13};
Physical Curve("left",      34) = {14};

Physical Surface("main",    35) = {111};

Mesh 2;
Save "square1.msh";