SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 4.1;
Mesh.Renumber = 1;
Mesh.SaveAll = 0;


radius = 12.5;

Sphere (1) = {0,0,-radius,radius,0,Pi/2};

Physical Volume("Rigid_dome") = {1};
Physical Surface("Dome") = {1};
Physical Surface("dome_fixed") = {2};

lc = 1;
Mesh.CharacteristicLengthFromPoints = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.MeshSizeFromCurvature = 0;


Mesh.Algorithm = 6;
Mesh 3;
Save "rigid_dome.msh";





