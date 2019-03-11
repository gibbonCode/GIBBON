function [Et,Vt,Ft]=hexahedral_hexagon_beam(r,ne,hz,nz,XYZ_centre)

%%


%% CREATING 2D SECTION MESH

%Basic regular quad mesh
[X,Y]=meshgrid(linspace(0,1,ne+1));
Z=zeros(size(X));
V=[X(:) Y(:) Z(:)];
[F,V] = surf2patch(X,Y,Z);

%Creating 3 sheared quad meshes to construct hexagon

%V1 
V1=V;
SH=eye(3,3); SH(1,2)=0; SH(2,1)=0.5;
V1=V1*SH; 
S=eye(3,3); S(1,1)=1; S(2,2)=sqrt(3/4);
V1=V1*S; 

%V2
a=pi/3;
R=[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 0;];
V2=V1*R; 

%V3
V3=V;
SH=eye(3,3); SH(1,2)=0; SH(2,1)=-0.5;
V3=V3*SH; 
S=eye(3,3); S(1,1)=1; S(2,2)=sqrt(3/4);
V3=V3*S; 
V3(:,1)=V3(:,1)+0.5; V3(:,2)=V3(:,2)+sqrt(3/4); 

%Composing hexagon
Vu=[V1;V2;V3];
Fu=[F;F+size(V1,1);F+2.*size(V1,1)];

%Removing double points
[Fu,Vu,~,~,~,~]=unique_patch(Fu,Vu,[]);

%Scaling radius
[THETA,R] = cart2pol(Vu(:,1),Vu(:,2));
[Vu(:,1),Vu(:,2)] = pol2cart(THETA,r.*R);

%% CREATING 3D EXTRUDED MESH

F=Fu; V=Vu; 

z_range=linspace(0,hz,nz+1);
Vt=repmat(V,numel(z_range),1);
Z_add=ones(size(V,1),1)*z_range; Z_add=Z_add(:);
Vt(:,3)=Vt(:,3)+Z_add;
Vt=Vt-ones(size(Vt,1),1)*mean(Vt,1); %centering around mean
Vt=Vt+ones(size(Vt,1),1)*XYZ_centre; %Translate to desired centre

Et=[];
for iz=1:1:numel(z_range)-1
    %       bottom                top 
    Et=[Et; F+(size(V,1).*(iz-1)) F+(size(V,1).*iz)];
end

% f_order=[4 3 2 1 8 7 6 5];
% Et=Et(:,f_order);

[Ft]=hex2patch(Et);

end
 
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
