Point(1) = {-2.0, -0.3, 0, 1.0};
Point(2) = {-0.5, -0.3, 0, 1.0};
Point(3) = {1.5, -0.3, 0, 1.0};
Point(4) = {1.5, 0.3, 0, 1.0};
Point(5) = {-0.5, 0.3, 0, 1.0};
Point(6) = {-2.0, 0.3, 0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(1) = {1,2,3,4,5,6};

Plane Surface(1) = {1};

// inflow
Physical Line(101) = {6};

// lower wall
Physical Line(0) = {1,2};

// outflow
Physical Line(102) = {3};

// upper wall
Physical Line(1) = {4,5};


Physical Surface(13) = {1};
Recombine Surface {1};
