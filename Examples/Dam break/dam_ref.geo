
// Geometry definition
Point(1) = {0, 0, 0}; 
Point(2) = {95, 0, 0}; 
Point(3) = {95, 95, 0}; 
Point(4) = {105, 95, 0}; 
Point(5) = {105, 0, 0}; 
Point(6) = {200, 0, 0};
Point(7) = {200, 200, 0};
Point(8) = {105, 200, 0};
Point(9) = {105, 170, 0};
Point(10) = {95, 170, 0};
Point(11) = {95, 200, 0};
Point(12) = {0, 200, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};
Line(13) = {3, 10};
Line(14) = {4, 9};

Line Loop(15) = {1, 2, 13, 10, 11, 12};
Plane Surface(16) = {15};
Line Loop(17) = {3,14,9,-13};
Plane Surface(18) = {17};
Line Loop(19) = {4,5,6,7,8,-14};
Plane Surface(20) = {19};

// Number of elements per 5m (gap in the middle)
Vp = 1;
Transfinite Line {1,5,7,11} = Vp*95/5+1 Using Progression 1; // horizontal bottom+top
Transfinite Line {6,12} = Vp*200/5+1 Using Progression 1; // vertical left/right
Transfinite Line {2,4} = Vp*95/5+1 Using Progression 1; // vertical middle bottom
Transfinite Line {8,10} = Vp*30/5+1 Using Progression 1; // vertical middle top
Transfinite Line {13,14} = Vp*75/5+1 Using Progression 1; // vertical gap
Transfinite Line {3,9} = Vp*10/5+1 Using Progression 1; // horizontal gap

Transfinite Surface {16} = {1, 2, 11,12} Alternated;
Transfinite Surface {18} = {3,4,9,10} Alternated;
Transfinite Surface {20} = {6,7,8,5} Alternated;

Physical Line(19) = {1,2,3,4,5,6,7,8,9,10,11,12};
Physical Surface(22) = {18,20};
Physical Surface(23) = {16};