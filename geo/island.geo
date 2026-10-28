// Briggs et al. (1995) conical island, basin coordinates of the experiment.
// 25 m x 28.2 m, 0.32 m of water, island centred at (12.96, 13.80): toe
// radius 3.6 m, crest radius 1.1 m, 0.625 m high. The toe, the still
// waterline and the crest are embedded so the mesh follows them.
DefineConstant[ lc = 0.25 ];   // offshore
DefineConstant[ li = 0.07 ];   // on and around the island
xc = 12.96; yc = 13.80;

Point(1) = {0, 0, 0, lc};
Point(2) = {25, 0, 0, lc};
Point(3) = {25, 28.2, 0, lc};
Point(4) = {0, 28.2, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Point(10) = {xc, yc, 0, li};
r[] = {3.6, 2.32, 1.1};
For i In {0:2}
  p1 = newp; Point(p1) = {xc + r[i], yc, 0, li};
  p2 = newp; Point(p2) = {xc - r[i], yc, 0, li};
  c1 = newc; Circle(c1) = {p1, 10, p2};
  c2 = newc; Circle(c2) = {p2, 10, p1};
  Curve{c1, c2} In Surface{1};
EndFor
Point{10} In Surface{1};

Field[1] = Distance;
Field[1].PointsList = {10};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = li;
Field[2].SizeMax = lc;
Field[2].DistMin = 4.2;
Field[2].DistMax = 8.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;

Physical Curve("wave") = {4};
Physical Curve("open") = {1, 2, 3};
Physical Surface("water") = {1};
