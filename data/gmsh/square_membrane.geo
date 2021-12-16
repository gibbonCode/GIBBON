
SetFactory("OpenCASCADE");
Mesh.MshFileVersion = 4.1;
Mesh.Renumber = 1;
Mesh.SaveAll = 0;


membrane_L = 50;
membrane_W = 50;
membrane_H = .5;

Rectangle (1) = {-membrane_L/2, -membrane_W/2, 0, membrane_L,membrane_W};
Extrude { 0, 0, membrane_H  }{ Surface{1};}

Physical Volume("Membrane") = {1};
Physical Surface("membrane_bottom") = {1};
Physical Surface("membrane_top") = {6};
Physical Surface("membrane_fixed") = {2,3,4,5};

Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.MeshSizeFromCurvature = 0;

Transfinite Curve{1,2,3,4,7,9,11,12} = 36;
Transfinite Curve{5,6,8,10} = 4;
Transfinite Surface{1,2,3,5,6};
Transfinite Surface{4} Right;
Transfinite Volume{1};

Mesh.Algorithm = 6;
Mesh 3;
Save "square_membrane.msh";
