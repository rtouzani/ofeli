// Gmsh Geometry file for meshing a unit square

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(5) = {1, 2, 3, 4};
Plane Surface(14) = {5};
Physical Point(1) = {1,2,3,4};
Physical Line(1) = {1,2,3,4};
Physical Surface(1) = {14};
