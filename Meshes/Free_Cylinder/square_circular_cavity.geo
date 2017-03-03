// points on the circle
Point(1) = {0.25, 0.0, 0, 1.0};
Point(2) = {0.0, 0.0, 0, 1.0};
Point(3) = {-0.25, 0.0, 0, 1.0};

//points on the square
Point(4) = {-2.0, -0.5, 0, 1.0};
Point(5) = {2.0, -0.5, 0, 1.0};
Point(6) = {2.0, 0.5, 0, 1.0};
Point(7) = {-2.0, 0.5, 0, 1.0};

//circular loops
Circle(1) = {1, 2, 3};
Circle(2) = {3, 2, 1};

// straight lines
Line(3) = {7, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};

// line loops
Line Loop(7) = {6, 3, 4, 5};
Line Loop(8) = {1, 2};

Plane Surface(9) = {7, 8};

//physical entities
//inflow
Physical Line(101) = {3};
// specular wall
Physical Line(50) = {4,6};
// outflow
Physical Line(102) = {5};
//inner circle
Physical Line(0) = {1,2};

Physical Surface(14) = {9};

Recombine Surface {9};
