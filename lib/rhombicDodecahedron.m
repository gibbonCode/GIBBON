function [F,V]=rhombicDodecahedron(r)

% function [F,V]=rhombicDodecahedron(r)
% ------------------------------------------------------------------------
% This function creates a rhombic dodecahedron with a radius r. 
%
% ------------------------------------------------------------------------


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

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
