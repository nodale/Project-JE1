SetFactory("OpenCASCADE");
Geometry.OCCParallel = 0;

B = 3;           // Number of blades
theta = 2*Pi/B;   // Sector angle in degrees
H = 0.2;         // Cylinder height
R = 0.4;         // Cylinder outer radius
Rhub = 0.032;      // Hub radius

b() = ShapeFromFile("blade.brep");
Surface Loop(10) = {b()};
Volume(10) = {10};

Cylinder(1) = {0,0,-H/2, 0,0,H, R, theta};

Cylinder(2) = {0,0,-H/2, 0,0,H, Rhub, theta};

sector() = BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete;};

Rotate	{{0,0,1}, {0,0,0}, -theta/2} { Volume{sector()};}

fin() = BooleanDifference{ Volume{sector()}; Delete; }{ Volume{10}; Delete;};

//Mesh 3;


