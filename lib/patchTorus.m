function [F,V]=patchTorus(r,nr,rc,nc,patchType)

pRange=linspace(0,2*pi,nr+1);
tRange=linspace(0,2*pi,nc+1);
[P,T]=ndgrid(pRange,tRange);
[Ip,~]=ndgrid(1:1:nr,1:1:nc);
P=P(1:end-1,1:end-1); T=T(1:end-1,1:end-1); %Crop away double points

if strcmp('tri',patchType) || strcmp('honey',patchType) 
    %Offset mesh in T direction to obtain aproximate equilateral triangular mesh
    pointSpacing=abs(tRange(1)-tRange(2));
    T(1:2:end,:)=T(1:2:end,:)+(pointSpacing/2);
end

%Derive coordinates
X=(rc+(r*cos(P))).*cos(T);
Y=(rc+(r*cos(P))).*sin(T);
Z=r*sin(P);
C=Ip;

%Convert to quadrilateral patch data
[F,V,C] = surf2patch(X,Y,Z,C);
[C]=vertexToFaceMeasure(F,C);
C=round(C);

%Merge ends
I=[(1:size(Z,1)-1)' (1:size(Z,1)-1)' (2:size(Z,1))' (2:size(Z,1))'];
J=[size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1)];
F_sub1=sub2ind(size(Z),I,J);
F=[F;F_sub1];
C=[C;(min(C(:)):1:max(C(:)))'];

%Merge side
I=[ones(size(Z,2)-1,1) size(Z,1).*ones(size(Z,2)-1,1) size(Z,1).*ones(size(Z,2)-1,1) ones(size(Z,2)-1,1)];
I(end+1,:)=I(end,:);
J=[(2:size(Z,2))' (2:size(Z,2))' (1:size(Z,2)-1)' (1:size(Z,2)-1)'];
J(end+1,:)=[J(1,3:4) J(end,1:2)];
F_sub2=sub2ind(size(Z),I,J);
F=[F;fliplr(F_sub2)];
if iseven(nr)
    C=[C;(min(C(:))-2)*ones(size(F_sub2,1),1)];
else
     C=[C;(min(C(:))-1)*ones(size(F_sub2,1),1)];
end

switch patchType
    case 'quad' %Quadrilaterial elements
        
    case {'tri','tri_q','honey'}
        logicSlashType=iseven(C);
        F1=[F(logicSlashType,2) F(logicSlashType,3) F(logicSlashType,1); F(logicSlashType,3) F(logicSlashType,4) F(logicSlashType,1)];
        F2=[F(~logicSlashType,2) F(~logicSlashType,3) F(~logicSlashType,4); F(~logicSlashType,4) F(~logicSlashType,1) F(~logicSlashType,2)];
        F=[F1;F2];
    otherwise
        error('wrong input for argument ptype, valid inputs are quad and tri');
end

if strcmp('honey',patchType) 
    % GET DUAL FOR HONEY-COMB
    [V,Fd]=patch_dual(V,F);
    numVert=cellfun(@(x) size(x,2),Fd);
    F=Fd{numVert==6};
end
    
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
