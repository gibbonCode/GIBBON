function [E,V]=rhombicDodecahedronHex(r)

%% Get rhombic dodecahedron

[F,V]=rhombicDodecahedron(r);

%% Add centre point

V(end+1,:)=zeros(1,3);

%% Define the elements

E1=[fliplr(F(1,:)) fliplr([size(V,1) F(5,3) F(5,4) F(8,1)]) ];
E2=[fliplr(F(2,:))  fliplr([F(9,2) F(10,3) size(V,1) F(9,1)]) ];
E3=[fliplr(F(3,:)) fliplr([size(V,1)  F(7,3) F(6,4) F(6,1)]) ];
E4=[fliplr((F(4,:))) fliplr([F(12,2) F(12,3) size(V,1) F(11,1)])  ];
E=[E1;E2;E3;E4];

