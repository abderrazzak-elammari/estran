// Thacker (1981) paraboloid, 4 m x 4 m, closed basin.
DefineConstant[ lc = 0.04 ];

Point(1) = {0, 0, 0, lc};
Point(2) = {4, 0, 0, lc};
Point(3) = {4, 4, 0, lc};
Point(4) = {0, 4, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Curve("wall") = {1, 2, 3, 4};
Physical Surface("water") = {1};
