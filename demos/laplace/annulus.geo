// Gmsh Geometry file for meshing an annulus

r = 1;
R = 2;
h = 0.1;

Point(1) = {0, 0, 0};
Point(2) = { r, 0, 0, h};
Point(3) = {-r, 0, 0, h};
Point(4) = { R, 0, 0, h};
Point(5) = {-R, 0, 0, h};
Circle(1) = {2,1,3};
Circle(2) = {3,1,2};
Circle(3) = {4,1,5};
Circle(4) = {5,1,4};
Curve Loop(4) = {1,2};
Curve Loop(5) = {3,4};
Physical Point(1) = {2,3};
Physical Point(2) = {4,5};
Physical Curve(1) = {1,2};
Physical Curve(2) = {3,4};
Plane Surface(5) = {4,-5};
Physical Surface(1) = {5};
