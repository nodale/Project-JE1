SetFactory("OpenCASCADE");
Geometry.OCCParallel = 0;
Geometry.OCCBooleanPreserveNumbering = 0;

B = 3;           // Number of blades
theta = 2*Pi/B;   // Sector angle in degrees
H = 0.14;         // Cylinder height
R = 0.17;         // Cylinder outer radius
Rhub = 0.048;      // Hub radius

Mesh.Algorithm = 2;
Mesh.Algorithm3D = 10;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 1;
Mesh.MinimumElementsPerTwoPi = 100;
Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.ElementOrder = 1;
Mesh.HighOrderOptimize = 2;
Mesh.MshFileVersion = 2.2;
Mesh.RecombineAll = 0;

Cylinder(1) = {-0.012,0,-H/2, 0,0,H, R, theta};
Cylinder(2) = {-0.012,0,-H/2, 0,0,H, Rhub, theta};

sector() = BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete;};

Rotate	{{0,0,1}, {0,0,0}, -theta/2} { Volume{sector()};}

Physical Surface("cyclic_1") = {5};
Physical Surface("cyclic_2") = {3};
Physical Surface("air") = {6, 2, 4, 1};

b() = ShapeFromFile("blade.brep");
Surface Loop(10) = {b()};
b[] = {b()};
Physical Surface("blade") = {b()};
Volume(10) = {10};
fin() = BooleanDifference{ Volume{sector()}; Delete; }{ Volume{10}; Delete;};
Physical Volume("all") = {fin()};

Field[10] = BoundaryLayer;
Field[10].CurvesList = {b[]};
Field[10].hwall_n = 1e-5;
Field[10].thickness = 1e-3;
Field[10].ratio = 1.2;
Field[10].Quads = 0; // IMPORTANT
Background Field = 10;


Point(10000) = {0,0,0, 1.0}; // reference point
Field[1] = Distance;
Field[1].NodesList = {10000};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.10; 
Field[2].LcMax = 0.18;  
Field[2].DistMin = 0.01;
Field[2].DistMax = 0.04;
Background Field = 2;

Transfinite Line{2, 4, 7, 8, 5, 10} = 30;
Transfinite Curve{3, 1} = 40;
Transfinite Line{6, 9, 11, 12} = 40;

Mesh 2;
//RecombineMesh;
Mesh 3;





