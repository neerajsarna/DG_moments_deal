// same geometry as box_cylinder.geo but in 2D
radius = 0.5;

//half of the edge length of the bounding box
half_length = 4.0;
char_length = 0.1;

// vertices of the square
Point(1) = {-half_length,-half_length,0,1.0};
Point(2) = {half_length,-half_length,0,1.0};
Point(3) = {half_length,half_length,0,1.0};
Point(4) = {-half_length,half_length,0,1.0};

// points for the Circle
Point(5) = {radius,0,0,1.0};
Point(6) = {-radius,0,0,1.0};
Point(7) = {0,0,0,1.0};

// Lines for the Square
Line(1) = {1,2};	// bottom
Line(2) = {2,3};	// right
Line(3) = {3,4};	// top
Line(4) = {4,1};	// left

// Circles
Circle(5) = {5,7,6};
Circle(6) = {6,7,5};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6};

Plane Surface(1) = {1,2};

// all the physcial lines and surfaces
Physical Surface(1) = {1};

// defining the boundaries of the domain
Physical Line(101) = {4};
Physical Line(102) = {2};
Physical Line(50) = {1,3};
//Physical Line(1) = {3};
Physical Line(2) = {5,6};

Recombine Surface{1};
