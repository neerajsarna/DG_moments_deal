// points for the square
side_length = 1.0;
char_length = 0.5;

Point(1) = {0,0,0,char_length};
Point(2) = {side_length,0,0,char_length};
Point(3) = {side_length,side_length,0,char_length};
Point(4) = {0,side_length,0,char_length};

// construct different lengths
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Recombine Surface{1};

Extrude{0,0,15}{
	Surface{1};
	Layers{6};
	Recombine;
}

Physical Surface(27) = {25, 26, 21, 17, 13, 1};
Physical Volume(28) = {1};
