// points for the square
side_length = 1.0;

Point(1) = {0,0,0,1.0};
Point(2) = {side_length,0,0,1};
Point(3) = {side_length,side_length,0,1};
Point(4) = {0,side_length,0,1};

// construct different lengths
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Line(50) = {1,3};
Physical Line(101) = {4};
Physical Line(102) = {2};

Recombine Surface{1};