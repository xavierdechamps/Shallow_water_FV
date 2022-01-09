
bottom = 5;
hexit = 10;
angle = 10*Pi/180;

x_ramp = 15;
Point(1) = {0, 0, 0};
Point(2) = {bottom, 0, 0};
Point(3) = {bottom+x_ramp, x_ramp*Tan(angle), 0};
Point(4) = {bottom+x_ramp, hexit, 0};
Point(5) = {0, hexit, 0};
Point(6) = {bottom, hexit, 0};

Line(1) = {1,2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 6};
Line(5) = {6, 5};
Line(6) = {5, 1};
Line(7) = {6, 2};

Line Loop(1) = {1, -7, 5, 6};
Plane Surface(1) = {1};
Line Loop(2) = {2, 3, 4, 7};
Plane Surface(2) = {2};

Vp = 20;
Hp1 = Vp/2;
Hp2 = Hp1*3;
Transfinite Line {1, 5} = Hp1 Using Progression 1; // horizontal bottom+top
Transfinite Line {2, 4} = Hp2 Using Progression 1; // slope bottom+top
Transfinite Line {3,6,7} = Vp Using Progression 1; // vertical inflow+outflow
Transfinite Surface {1} = {1, 2, 6, 5} Alternated;
Transfinite Surface {2} = {3, 4, 6,2} Alternated;

Physical Line(20) = {1,2,4,5};
Physical Line(21) = {3};
Physical Line(22) = {6};
Physical Surface(1) = {1,2};
//Physical Surface(13) = {2};
