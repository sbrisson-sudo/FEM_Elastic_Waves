// GEOMETRY

// mesh parameters
h1 = 10.0;
h2 = Sqrt(2)/2*h1;

// points definition
Point(1) = {0.0     , 0.0       , 0.0};
Point(2) = {h1      , 0.0       , 0.0};
Point(3) = {h1+h2   , h2        , 0.0};
Point(4) = {h1+h2   , h1+h2     , 0.0};
Point(5) = {h1      , h1+2*h2   , 0.0};
Point(6) = {0.0     , h1+2*h2   , 0.0};
Point(7) = {-h2     , h1+h2     , 0.0};
Point(8) = {-h2     , h2        , 0.0};

// curves definition (Line, Curve, Spline...)
Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 5};
Line(15) = {5, 6};
Line(16) = {6, 7};
Line(17) = {7, 8};
Line(18) = {8, 1};
Curve Loop(21) = {11,12,13,14,15,16,17,18};

// surfaces definition
Plane Surface(111) = {21};
Recombine Surface {111};    // for quadrangles
Mesh.ElementOrder = 4;


// physical curves
Physical Curve("bottom",21) = {11};
Physical Curve("rigth",22) = {12,13,14};
Physical Curve("top",23) = {15};
Physical Curve("left",24) = {16,17,18};

// physical surfaces
Physical Surface("main",211) = {111};

// MESHING
Mesh.MeshSizeFactor=0.5; // densify the mesh
// Mesh.MshFileVersion = 4;
Mesh 2;                  // mesh in 2D
Save "octogon4.msh";      // save file


