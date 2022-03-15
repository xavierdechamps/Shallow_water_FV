
Lu = 1000 ;   // Length upstream
Ld = 4000 ;    // Length downstream
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

// Diverging area right
Line Loop(4) = {4,23,8,-22};
Plane Surface(4) = {4};

// Right
Line Loop(5) = {5, 6, 7, -23};
Plane Surface(5) = {5};

// lc = 3.;
Vp = 31; // Vertical
Hleft = 51; // Horizontal left
Hconv = 21; // converging area
Hcenter = 11 ; // Center
Hdiv = 21; // diverging area
Hright = 51; // Horizontal right

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
Include "channel_contraction_Cueto_gradient_Ld4000.pos";


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
Field[15].XMin = Lu-50;
Field[15].XMax = Lu+50;
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
Field[19].XMax = Lu+0.5*Lb+0.4*Le;
Field[19].YMax = y1+b2/4;
Field[19].YMin = 0;
Field[19].ZMax = 1;
Field[19].ZMin = -1;

Field[20] = Box; // Corners right top
Field[20].VIn  = size_near / raff2;
Field[20].VOut = size_far_away;
Field[20].XMin = Lu+0.5*Lb;
Field[20].XMax = Lu+0.5*Lb+0.4*Le;
Field[20].YMax = b1;
Field[20].YMin = b1-y1-b2/4;
Field[20].ZMax = 1;
Field[20].ZMin = -1;

Field[50] = Min;
Field[50].FieldsList = {14,15,16,17,18,19,20};
Background Field = 50;
// Don't extend the elements sizes from the boundary inside the domain
// Mesh.CharacteristicLengthExtendFromBoundary = 0;

Mesh.Algorithm = 5 ; // Delaunay will handle complex mesh size
// fields better - in particular size fields with large element size gradients

If (quad_mesh == 1)
    Mesh.RecombinationAlgorithm = 0; // 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad
    Recombine Surface {1,2,3,4,5,6};
EndIf

Physical Line(3) = {1,2,3,4,5,7,8,9,10,11}; // Walls
Physical Line(2) = {6}; // Outlet
Physical Line(1) = {12}; // Inlet
Physical Surface(1)={1,2,3,4,5,6};