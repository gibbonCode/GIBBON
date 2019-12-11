function [varargout]=rhombicDodecahedronMesh(varargin)

% function [E,V,C,F,FC]=rhombicDodecahedronMesh(r,nCopies)
% ------------------------------------------------------------------------
% Creates a rhombic dodecahedron mesh where r sets the radias and nCopies
% (a 1x3 vector) sets the number of copies in the x, y, and z direction.
% The output consists of:
%
% Fc_Q, Fc_T: the quadrilateral and triangular face cell arrays (1 cell
% entry per element).
%
% Ft_Q, Ft,T: the quadrilateral and triangular face arrays
%
% Ct_Q, Ct,T: color/label data for the face arrays
%
% Vt: the vertex array
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2019/02/08 Created
% 2019/10/13 Changed orientation for gridding to create more regular grid
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        r=varargin{1};
        nCopies=2;
    case 2
        r=varargin{1};
        nCopies=varargin{2};
end

if isempty(nCopies)
    nCopies=2;
end

if numel(nCopies)==1
    nCopies=nCopies*ones(1,3);
end

%%

% Get single rhombic dodecahedron
[~,Vs]=rhombicDodecahedron(r);
V_offsets=2*r*eye(3,3); %Offset vectors

% Copy over rhombic dodecahedron in grid like fashion 4 times and shift
% some of the copies to nest them together. 
nCopies1=nCopies;
[Vg1]=gridClone(Vs,nCopies1,V_offsets);

nCopies2=nCopies;
nCopies2([1,2])=nCopies2([1,2])-1;
[Vg2]=gridClone(Vs,nCopies2,V_offsets);

nCopies3=nCopies;
nCopies3([1,3])=nCopies3([1,3])-1;
[Vg3]=gridClone(Vs,nCopies3,V_offsets);

nCopies4=nCopies;
nCopies4([2,3])=nCopies4([2,3])-1;
[Vg4]=gridClone(Vs,nCopies4,V_offsets);

%%

Vg2(:,1)=Vg2(:,1)+r;
Vg2(:,2)=Vg2(:,2)+r;

Vg3(:,1)=Vg3(:,1)+r;
Vg3(:,3)=Vg3(:,3)+r;

Vg4(:,2)=Vg4(:,2)+r;
Vg4(:,3)=Vg4(:,3)+r;

V=[Vg1;Vg2;Vg3;Vg4];

%% Create element and face arrays
E=reshape((1:1:size(V,1)),14,size(V,1)/14)';
C=(1:1:size(E,1))'; %Element colors
[F,CF]=element2patch(E,C);

%% Merging point and fix element and face indices
[F,V,~,indFix]=mergeVertices(F,V);
E=indFix(E);

%% Collect output
varargout{1}=E;
varargout{2}=V;
varargout{3}=C;
varargout{4}=F;
varargout{5}=CF;

end

%%

function [Vg]=gridClone(Vs,nCopies,V_offsets)

nTotal=prod(nCopies); %Total number of copies
indC=ones(size(Vs,1),1)*(1:1:nTotal);
indC=indC(:);
[I,J,K] = ind2sub(nCopies,indC);
I=I-1; J=J-1; K=K-1;

%Defining offsets
D1=I*V_offsets(1,:);
D2=J*V_offsets(2,:);
D3=K*V_offsets(3,:);

%Defining vertices matrix
Vg=repmat(Vs,nTotal,1)+(D1+D2+D3);

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
