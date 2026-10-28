// LeVeque (1998) bump: 2 m x 1 m, still water at 1 m over a 0.8 m bump.
// Open ends at x = 0 and x = 2, walls along y.
DefineConstant[ lc = 0.02 ];

Point(1) = {0, 0, 0, lc};
Point(2) = {2, 0, 0, lc};
Point(3) = {2, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Curve("wall") = {1, 3};
Physical Curve("open") = {2, 4};
Physical Surface("water") = {1};
