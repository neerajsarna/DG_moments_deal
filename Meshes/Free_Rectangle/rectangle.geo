left_edge = -2.0;
right_edge = 2.0;
bottom_edge = -0.5;
top_edge = 0.5;

// left bottom
Point(1) = {left_edge,bottom_edge,0,1.0};

// right bottom
Point(2) = {right_edge,bottom_edge,0,1.0};

// right top
Point(3) = {right_edge,top_edge,0,1.0};

// left top
Point(4) = {left_edge,top_edge,0,1.0};

// bottom edge
Line(1) = {1,2};

// right edge
Line(2) = {2,3};

// top edge 
Line(3) = {3,4};

// left edge
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

// inflow with a prescribed velocity
Physical Line(103) = {4};

// lower wall
Physical Line(50) = {1,3};

// outflow
Physical Line(102) = {2};


Physical Surface(13) = {1};
Recombine Surface {1};
Transfinite Line {1,3} = 100 Using Progression 1;
Transfinite Line {2,4} = 100 Using Progression 1;
Transfinite Surface {1};
