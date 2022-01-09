L=1000;
H=2;

// Box definition
Point(1) = {0, H, 0};
Point(2) = {0, 0, 0};
Point(3) = {L, H, 0};
Point(4) = {L, 0, 0};

Line(13) = {3, 1};
Line(14) = {2, 4};
Line(15) = {4, 3};
Line(16) = {1, 2};
Line Loop(17) = {13, 16, 14, 15};
Plane Surface(18) = {17};

lc = 0.4;
Lp = 1001;
Hp = 3;
Characteristic Length {1, 2, 4, 3} = lc;
Transfinite Line {13, 14} = Lp Using Progression 1;
Transfinite Line {16, 15} = Hp Using Progression 1;
Transfinite Surface {18} = {3, 4, 2, 1} Alternated;


Physical Line(1) = {16};  // Inflow
Physical Line(2) = {15}; // Outflow
Physical Line(3) = {13, 14}; // Wall
Physical Surface(1)={18};

