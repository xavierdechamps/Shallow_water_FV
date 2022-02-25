
Lu = 1000 ;   // Length upstream
Ld = 75 ;     // Length downstream
Lb = 10 ;     // Length of narrowest section
Le = 20 ;     // Length of expansion
b1 = 40 ;     // Largest width
b2 = 26.5 ;   // Smallest width

quad_mesh = 1;  // [0] Triangular mesh or [1] rectangular mesh

y1 = (b1-b2)*0.5;

// Box definition
Point(1) = {0 ,           0 ,   0};
Point(2) = {Lu-0.5*Lb-Le, 0 ,   0};
Point(3) = {Lu-0.5*Lb,    y1,   0};
Point(4) = {Lu+0.5*Lb,    y1,   0};
Point(5) = {Lu+0.5*Lb+Le, 0,    0};
Point(6) = {Lu+Ld,        0,    0};
Point(7) = {Lu+Ld,        b1,    0};
Point(8) = {Lu+0.5*Lb+Le, b1,    0};
Point(9) = {Lu+0.5*Lb,    b1-y1,   0};
Point(10) = {Lu-0.5*Lb,    b1-y1,   0};
Point(11) = {Lu-0.5*Lb-Le, b1 ,   0};
Point(12) = {0 ,           b1 ,   0};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 5};
Line(5) = {5 , 6};
Line(6) = {6 , 7};
Line(7) = {7 , 8};
Line(8) = {8 , 9};
Line(9) = {9 , 10};
Line(10) = {10 , 11};
Line(11) = {11 , 12};
Line(12) = {12 , 1};

Line(20) = {2 , 11};
Line(21) = {3 , 10};
Line(22) = {4 , 9};
Line(23) = {5 , 8};

// Left
Line Loop(1) = {1 , 20, 11, 12};
Plane Surface(1) = {1};

// Contraction left
Line Loop(2) = {2 , 21, 10, -20};
Plane Surface(2) = {2};

// Central area
Line Loop(3) = {3, 22, 9, -21};
Plane Surface(3) = {3};

// Shock
angle_exp = Atan( 0.5 * (b1-b2) / Le );
angle_shock1 = Pi*0.5 - angle_exp * 1.3; 
angle_shock2 = -Pi*0.5 + angle_exp ; 
x1 = Lu+0.5*Lb + 16.5 ;
y1 = 0 ;
y2 = b1*0.5;
x2 = (y2-y1)/Tan(angle_shock1) + x1 ;
Point(13) = {x1, y1, 0, 1.0};
Point(14) = {x2, y2, 0, 1.0};
Point(15) = {x1, b1, 0, 1.0};
Line(24) = {13 , 14};
Line(25) = {15 , 14};

Extrude {0, 0, 1} { Curve{24,25}; } // surface 29-33
Intersect Curve {4,23} Surface {29} // Points 19 (bottom) and 20 center
Intersect Curve {8,23} Surface {33} // Points 21 (top) and 22 center
Recursive Delete { Surface{29}; Surface{33};}
Recursive Delete { Curve{23}; Curve{4}; Curve{8};}

Line(24) = {19 , 20};
Line(25) = {21 , 22};
Line(26) = {4 , 19};
Line(27) = {19 , 5};
Line(28) = {8 , 21};
Line(29) = {21 , 9};
Point(15) = {x1-9, y2, 0, 1.0};
Circle(30) = {20, 15, 22};

// Diverging area
Line Loop(4) = {26,24,30,-25,29,-22};
Plane Surface(4) = {4};

// Right
Line Loop(5) = {27,5,6,7,28,25,-30,-24};
Plane Surface(5) = {5};

// lc = 3.;
// Vp = 31; // Vertical
// Hleft = 51; // Horizontal left
// Hconv = 21; // converging area
// Hcenter = 11 ; // Center
// Hdiv = 21; // diverging area
// Hright = 51; // Horizontal right

size_far_away = 40;
size_middle = 4;
size_near = 1;
If (quad_mesh == 0) // Triangles
    raff1 = 2;
    raff2 = 4;
    raff3 = 8;
Else               // Rectangles
    raff1 = 3;
    raff2 = 6;
    raff3 = 10;
EndIf

//--   Characteristic lengths at key points
Characteristic Length {1, 12, 6,7} = size_far_away;
// Characteristic Length {2,3,4,5, 8,9,10,11} = size_near;

//--   Merge a post-processing view containing the target anisotropic mesh sizes
Include "channel_contraction_Cueto_gradient_Ld75.pos";

y1 = (b1-b2)*0.5;
// Attractor Line 21 
Field[1] = Distance;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {21};

// Attractor Shock
Field[3] = Distance;
Field[3].NNodesByEdge = 200;
Field[3].EdgesList = {24,30,25};

Field[11] = Threshold;
Field[11].IField = 1;
Field[11].LcMin = size_near / raff1;
Field[11].LcMax = size_near / raff1;
Field[11].DistMin = 0.01;
Field[11].DistMax = 5;
Field[11].StopAtDistMax = 1;

Field[13] = Threshold;
Field[13].IField = 3;
Field[13].LcMin = size_near / raff3;
Field[13].LcMax = size_near / raff2;
Field[13].DistMin = 0.0001;
Field[13].DistMax = 2;
Field[13].StopAtDistMax = 1;

Field[14] = Box; // Refine mesh in the inflow
Field[14].VIn = size_far_away/4;
Field[14].VOut = size_far_away;
Field[14].XMin = 700;
Field[14].XMax = Lu+Lb/2+Le+4000;
Field[14].YMax =  10+b1;
Field[14].YMin = -10;
Field[14].ZMax = 1;
Field[14].ZMin = -1;

Field[15] = Box; // Global refinement of central area
Field[15].VIn = size_middle;
Field[15].VOut = size_far_away;
Field[15].XMin = 950;
Field[15].XMax = 1500;
Field[15].YMax =  10+b1;
Field[15].YMin = -10;
Field[15].ZMax = 1;
Field[15].ZMin = -1;

Field[16] = Box; // Finer refinement of central area
Field[16].VIn = size_near;
Field[16].VOut = size_far_away;
Field[16].XMin = Lu-0.5*Lb;
Field[16].XMax = Lu+0.5*Lb+Le*1.2;
Field[16].YMax =  10+b1;
Field[16].YMin = -10;
Field[16].ZMax = 1;
Field[16].ZMin = -1;

Field[17] = Box; // Corners left bottom
Field[17].VIn  = size_near / raff2;
Field[17].VOut = size_far_away;
Field[17].XMin = Lu-0.5*Lb;
Field[17].XMax = Lu;
Field[17].YMax = y1+b2/4;
Field[17].YMin = y1;
Field[17].ZMax = 1;
Field[17].ZMin = -1;

Field[18] = Box; // Corners left top
Field[18].VIn  = size_near / raff2;
Field[18].VOut = size_far_away;
Field[18].XMin = Lu-0.5*Lb;
Field[18].XMax = Lu;
Field[18].YMax = b1-y1;
Field[18].YMin = b1-y1-b2/4;
Field[18].ZMax = 1;
Field[18].ZMin = -1;

Field[19] = Box; // Corners right bottom
Field[19].VIn  = size_near / raff2;
Field[19].VOut = size_far_away;
Field[19].XMin = Lu+0.5*Lb;
Field[19].XMax = Lu+0.5*Lb+0.5*Le;
Field[19].YMax = y1+b2/4;
Field[19].YMin = 0;
Field[19].ZMax = 1;
Field[19].ZMin = -1;

Field[20] = Box; // Corners right top
Field[20].VIn  = size_near / raff2;
Field[20].VOut = size_far_away;
Field[20].XMin = Lu+0.5*Lb;
Field[20].XMax = Lu+0.5*Lb+0.5*Le;
Field[20].YMax = b1;
Field[20].YMin = b1-y1-b2/4;
Field[20].ZMax = 1;
Field[20].ZMin = -1;

Field[50] = Min;
Field[50].FieldsList = {11,13,14,15,16,17,18,19,20};
Background Field = 50;
// Don't extend the elements sizes from the boundary inside the domain
// Mesh.CharacteristicLengthExtendFromBoundary = 0;

Mesh.Algorithm = 5 ; // Delaunay will handle complex mesh size
// fields better - in particular size fields with large element size gradients

If (quad_mesh == 1)
    Mesh.RecombinationAlgorithm = 0; // 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad
    Recombine Surface {1,2,3,4,5,6};
EndIf

Physical Line(3) = {1,2,3,26,27,5,7,28,29,9,10,11}; // Walls
Physical Line(2) = {6}; // Outlet
Physical Line(1) = {12}; // Inlet
Physical Surface(1)={1,2,3,4,5}; 

