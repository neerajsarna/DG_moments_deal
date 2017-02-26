// Gmsh project created on Sun Feb 26 16:28:57 2017
Point(1) = {0.25, 0.0, 0, 1.0};
Point(2) = {0.0, 0.0, 0, 1.0};
Point(3) = {-0.25, 0.0, 0, 1.0};

Point(4) = {-10.0, -10.0, 0, 1.0};
Point(5) = {10.0, -10.0, 0, 1.0};
Point(6) = {10.0, 10.0, 0, 1.0};
Point(7) = {-10.0, 10.0, 0, 1.0};

Line(1) = {7, 4};
Line(2) = {4, 5};
Line(3) = {5, 6};
Line(4) = {6, 7};
Circle(5) = {1, 2, 3};
Circle(6) = {3, 2, 1};
Line Loop(7) = {1, 2, 3, 4};
Line Loop(8) = {5, 6};

Plane Surface(9) = {7, 8};


Physical Line(101) = {1};
Physical Line(0) = {2};
Physical Line(102) = {3};
Physical Line(1) = {4};
Physical Line(2) = {5, 6};

Recombine Surface {9};
