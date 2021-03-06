	// points on the circle
char_len = 0.5;
Point(1) = {0.25, 0.0, 0, char_len};
Point(2) = {0.0, 0.0, 0, char_len};
Point(3) = {-0.25, 0.0, 0, char_len};

//points on the square
Point(4) = {-2.0, -0.5, 0, char_len};
Point(5) = {2.0, -0.5, 0, char_len};
Point(6) = {2.0, 0.5, 0, char_len};
Point(7) = {-2.0, 0.5, 0, char_len};

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


Recombine Surface {9};

Physical Line(101) = {3};
Physical Line(0) = {4};
Physical Line(102) = {5};
Physical Line(1) = {6};
Physical Line(2) = {1, 2};
Physical Surface(15) = {9};
