radius = 0.5;

//half of the edge length of the bounding box
half_length = 4.0;
char_length = 0.1;

// points for the bounding cube
Point(1) = {-half_length, -half_length, 0, 1.0};
Point(2) = {half_length, -half_length, 0, 1.0};
Point(3) = {half_length, half_length, 0, 1.0};
Point(4) = {-half_length, half_length, 0, 1.0};
Point(5) = {0.0, 0.0, 0, 1.0};
Point(6) = {0.0, half_length, 0, 1.0};
Point(7) = {half_length, 0.0, 0, 1.0};
Point(8) = {-half_length, 0.0, 0, 1.0};
Point(9) = {0.0, -half_length, 0, 1.0};

// points for the circumference of the cylinder
Point(10) = {radius, 0.0, 0, char_length};
Point(11) = {0.0, radius, 0, char_length};
Point(12) = {-radius, 0.0, 0, char_length};
Point(13) = {0.0, -radius, 0, char_length};

// circular loops
Circle(1) = {10, 5, 11};
Circle(2) = {11, 5, 12};
Circle(3) = {12, 5, 13};
Circle(4) = {13, 5, 10};

Line(5) = {7, 3};
Line(6) = {3, 6};
Line(7) = {6, 4};
Line(8) = {4, 8};
Line(9) = {8, 1};
Line(10) = {1, 9};
Line(11) = {9, 2};
Line(12) = {2, 7};
Line(13) = {6, 11};
Line(14) = {10, 7};
Line(15) = {12, 8};
Line(16) = {13, 9};
Line Loop(17) = {13, -1, 14, 5, 6};
Plane Surface(18) = {17};
Line Loop(19) = {7, 8, -15, -2, -13};
Plane Surface(20) = {19};
Line Loop(21) = {9, 10, -16, -3, 15};
Plane Surface(22) = {21};
Line Loop(23) = {4, 14, -12, -11, -16};
Plane Surface(24) = {23};

// recombination to construct quad meshes
Recombine Surface {20, 18, 24, 22};

// extrusion in the z direction of length 1
Extrude {0, 0, 3} {
  Surface{20, 18, 24, 22};
  Layers{3};
  Recombine;
}

// definition of physical surfaces and volumes
Physical Surface(133) = {115, 38, 132, 51, 34, 20, 78, 105, 77, 73, 96, 46, 127, 65, 88, 119, 100, 18, 24, 22};
Physical Volume(134) = {4, 3, 2, 1};
