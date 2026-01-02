// Gmsh geometry generated from STL only

R  = 0.2;
Rh = 0.032;
H  = 0.1;
A  = 1.5708;
Ab = 0.785398;

Merge "blade.stl";

DefineConstant[
  angle = {40, Min 20, Max 120, Step 1, Name "Angle for surface detection"},
  forceParametrizablePatches = {1, Choices{0,1}, Name "Force creation of parametrizable patches"},
  includeBoundary = 1,
  curveAngle = 180
];

ClassifySurfaces{angle*Pi/180, includeBoundary, forceParametrizablePatches, curveAngle*Pi/180};
CreateGeometry;

// Push blade slightly to avoid numerical touching
Translate {0.0001, 0, 0} { Surface{:}; }

Surface Loop(1) = Surface{:};
Volume(1) = {1};

Physical Volume("blade") = {1};
Physical Surface("inlet") = Surface In BoundingBox{-R,-R,-H/2-1e-6, R,R,-H/2+1e-6};
Physical Surface("outlet") = Surface In BoundingBox{-R,-R,H/2-1e-6, R,R,H/2+1e-6};

