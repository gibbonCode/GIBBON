function [F,V]=rhombicDodecahedron(r)

%% Construct rhombic dodecahedron using cube and octahedron

%Get cube
[Vc,~]=platonic_solid(2,1);
Vc=Vc./(ones(size(Vc,1),1)*max(Vc,[],1)); %Scale coordinates

%Get octahedron
[Vo,~]=platonic_solid(3,1);
Vo=Vo.*2; %Scale coordinates
[R,~]=euler2DCM([0 0 0.25*pi]); %Define rotation
Vo=Vo*R; %Rotate coordinates

%Compose VERTICES of rhombic dodecahedron
V=[Vc;Vo];
V=r.*(V./(ones(size(V,1),1)*max(V,[],1))); %Scale radius

%Compose FACES of rhombic dodecahedron
sideFormat=[1 size(Vc,1)+1 5 size(Vc,1)+2];
topFormat=[5 size(Vc,1)+2 6 size(Vc,1)+6];
bottomFormat=[1 size(Vc,1)+2 2 size(Vc,1)+5];
F=[sideFormat;sideFormat+1;sideFormat+2; [sideFormat(1:3)+3 sideFormat(4)-1];...
   topFormat; [topFormat(1:3)+1 size(Vc,1)+6]; [topFormat(1:3)+2 size(Vc,1)+6]; [topFormat(1)+3 topFormat(2)-1 topFormat(3)-1 size(Vc,1)+6];...
   bottomFormat; [bottomFormat(1:3)+1 size(Vc,1)+5]; [bottomFormat(1:3)+2 size(Vc,1)+5]; [bottomFormat(1)+3 bottomFormat(2)-1 bottomFormat(3)-1 size(Vc,1)+5];];

%Fix face normals 
F(1:4,:)=F(1:4,[1 4 3 2]);
F(end-3:end,:)=F(end-3:end,[1 4 3 2]);

%TO DO: 
% 1)Remove need for fixing of face normals by just specifying the
% coordinates properly. 
% 2)Similarly remove need for scaling

