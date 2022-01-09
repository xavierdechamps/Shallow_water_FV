L2 = 22.234;
L1 = 25-L2;
L3 = 5.0;
H1 = 20.0;
H2 = 10.548;

y1 = (H1-H2)/2. ;

// Box definition
Point(1) = {0 ,       0 ,   0};
Point(2) = {L1,       0 ,   0};
Point(3) = {L1+L2,    y1,   0};
Point(4) = {L1+L2+L3, y1,   0};
Point(5) = {L1+L2+L3, y1+H2,0};
Point(6) = {L1+L2,    y1+H2,0};
Point(7) = {L1,       H1 ,  0};
Point(8) = {0,        H1 ,  0};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {3 , 4};
Line(4) = {4 , 5};
Line(5) = {5 , 6};
Line(6) = {6 , 7};
Line(7) = {7 , 8};
Line(8) = {8 , 1};
Line(9) = {6 , 3};
Line(10) = {7 , 2};

// Left
Line Loop(1) = {1 , -10, 7, 8};
Plane Surface(1) = {1};

// Middle
Line Loop(2) = {2 , -9, 6, 10};
Plane Surface(2) = {2};

// Right
Line Loop(3) = {3, 4, 5, 9};
Plane Surface(3) = {3};

lc = 3.;
Vp = 81; // Vertical
Hp = 80; // Horizontal middle
Hp2 = 21; // Horizontal left+right

uniform = 1 ;

If (uniform == 0)
    //--   Merge a post-processing view containing the target anisotropic mesh sizes
    Include "channel_contraction_Lai_sol.pos";

    Plugin(MathEval).Expression0 = "v0/60" ;
    Plugin(MathEval).View = 0 ;
    Plugin(MathEval).Run ;

    //--   Apply the view as the current background mesh
    Background Mesh View[1];
    //--   Characteristic lengths at key points
    Characteristic Length {1, 2, 3, 4, 5, 6, 7, 8} = lc;
Else
    Transfinite Line {4,8,9,10} = Vp Using Progression 1; // Vertical lines
    Transfinite Line {2,6} = Hp Using Progression 1; // Horizontal middle
    Transfinite Line {1,7} = Hp2 Using Progression 1; // Horizontal left
    Transfinite Line {3,5} = Hp2*4 Using Progression 1; // Horizontal right

    Transfinite Surface {1} = {2,7,8,1} Alternated;
    Transfinite Surface {2} = {2,3,6,7} Alternated;
    Transfinite Surface {3} = {4,5,6,3} Alternated;
EndIf

Physical Line(3) = {1,2,3,5,6,7}; // Walls
Physical Line(2) = {4}; // Outlet
Physical Line(1) = {8};
Physical Surface(1)={1,2,3};

