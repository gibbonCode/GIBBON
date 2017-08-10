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

 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
