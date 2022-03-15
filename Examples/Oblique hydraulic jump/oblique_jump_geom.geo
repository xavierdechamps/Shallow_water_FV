
bottom = 5;
hexit = 10;
angle = 10*Pi/180;

quad_mesh = 1;  // [0] Triangular mesh or [1] rectangular mesh

x_ramp = 15;
Point(1) = {0, 0, 0};
Point(2) = {bottom, 0, 0};
Point(3) = {bottom+x_ramp, x_ramp*Tan(angle), 0};
Point(4) = {bottom+x_ramp, hexit, 0};
Point(5) = {0, hexit, 0};
Point(6) = {bottom, hexit, 0};

Line(10) = {1,2};
Line(20) = {2, 3};
Line(30) = {3, 4};
Line(40) = {4, 6};
Line(50) = {6, 5};
Line(60) = {5, 1};
Line(70) = {6, 2};

Line Loop(1) = {10, -70, 50, 60};
Plane Surface(1) = {1};
Line Loop(2) = {20, 30, 40, 70};
Plane Surface(2) = {2};

Vp = 20;
Hp1 = Vp/2;
Hp2 = Hp1*3;
Transfinite Line {10, 50} = Hp1 Using Progression 1; // horizontal bottom+top
Transfinite Line {20, 40} = Hp2 Using Progression 1; // slope bottom+top
Transfinite Line {30,60,70} = Vp Using Progression 1; // vertical inflow+outflow
Transfinite Surface {1} = {1, 2, 6, 5} Alternated;
Transfinite Surface {2} = {3, 4, 6,2} Alternated;

If (quad_mesh == 1)
    Mesh.RecombinationAlgorithm = 0; // 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad
    Recombine Surface {1};
EndIf

Physical Line(20) = {10,20,40,50};
Physical Line(21) = {30};
Physical Line(22) = {60};
Physical Surface(1) = {1,2};
//Physical Surface(13) = {2};
