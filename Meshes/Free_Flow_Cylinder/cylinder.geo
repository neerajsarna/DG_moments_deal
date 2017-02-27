
Point(1) = {0.5, 0.0, 0, 1.0};
Point(2) = {0.0, 0.0, 0, 1.0};
Point(3) = {-0.5, 0.0, 0, 1.0};
Point(4) = {0.0, 5.0, 0, 1.0};
Point(5) = {0.0, -5.0, 0, 1.0};
Circle(1) = {1, 2, 3};
Circle(2) = {3, 2, 1};
Circle(3) = {4, 2, 5};
Circle(4) = {5, 2, 4};
Line Loop(5) = {3, 4};
Line Loop(6) = {1, 2};
Plane Surface(7) = {5, 6};
Recombine Surface {7};
Physical Surface(0) = {7};
Physical Line(101) = {3};
Physical Line(102) = {4};
Physical Line(0) = {1, 2};
