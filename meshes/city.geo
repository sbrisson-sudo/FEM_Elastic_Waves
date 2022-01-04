Mesh.ElementOrder = 4;

Lx = 1000.0;
Ly = 300.0;
dx = 25.0;

// main block

Point(1) = {0.0,        0.0,        0.0, dx};
Point(2) = {Lx,   0.0,        0.0, dx};
Point(3) = {Lx,   Ly,   0.0, dx};
Point(4) = {0.0,        Ly,   0.0, dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

// buildings blocks

// number of buildings
Nb = 4;                 
// left base positions
Xb = {600.0, 100.0, 225.0, 750.0};
// widths
Wb = {100.0, 100.0, 150.0, 200.0};    
// heights
Hb = {200.0, 100.0, 300.0, 100.0};    

For k In {0:Nb-1:1}

    Printf("k = %f", k );

    w = Wb[k];
    h = Hb[k];
    x = Xb[k];

    Point(5+4*k) = {x,  Ly, 0.0, dx};
    Point(6+4*k) = {x + w, Ly, 0.0, dx};
    Point(7+4*k) = {x + w, Ly + h, 0.0, dx};
    Point(8+4*k) = {x, Ly + h, 0.0, dx};

    Line(5+4*k) = {5+4*k, 6+4*k};
    Line(6+4*k) = {6+4*k, 7+4*k};
    Line(7+4*k) = {7+4*k, 8+4*k};
    Line(8+4*k) = {8+4*k, 5+4*k};

    Curve Loop(2+k) = {5+4*k, 6+4*k, 7+4*k, 8+4*k};

    Plane Surface(2+k) = {2+k};
    Transfinite Surface {2+k};
    Recombine Surface {2+k};

 EndFor


Mesh 2;
Coherence Mesh;

Save "city4.msh";
