// Dam break in a flat channel, 50 m x 2 m, dam at x = 25 m. The dam is a
// mesh line so the initial jump sits on cell edges.
DefineConstant[ lc = 0.25 ];
L = 50; W = 2; xd = 25;

Point(1) = {0, 0, 0, lc};
Point(2) = {xd, 0, 0, lc};
Point(3) = {L, 0, 0, lc};
Point(4) = {L, W, 0, lc};
Point(5) = {xd, W, 0, lc};
Point(6) = {0, W, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};
Curve Loop(1) = {1, 7, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 4, -7};
Plane Surface(2) = {2};

Physical Curve("wall") = {1, 2, 4, 5, 6};
Physical Curve("open") = {3};
Physical Surface("water") = {1, 2};
