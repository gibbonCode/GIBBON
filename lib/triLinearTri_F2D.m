function F_cell=triLinearTri_F2D(VX,Vx,TRI)

F_cell=cell(size(TRI,1),1);
for qc=1:1:numel(F_cell);
    indNow=TRI(qc,:);
    X=VX(indNow,:); %Initial coordinates
    x=Vx(indNow,:); %Current coordinates
    F=triLinearTri_subF2D(X,x); %The deformation gradient tensor
    F_cell{qc}=F; %Store in cell array
end

end

function F=triLinearTri_subF2D(X,x)

% Define shape function set N=[1-r-s; r; s;];

%Hardcoded derivative
dN_dRST =[-1 -1;...
           1  0;...
           0  1];

% Compute derivatives of initial position vectors with respect to shape functions
dX_dRST=zeros(2,2);
for q=1:1:size(dN_dRST,1);
    dX_dRST=dX_dRST+(X(q,:)'*dN_dRST(q,:));
end

% Compute derivatives of shape functions with respect to initial position vectors
dN_dX=zeros(2,2);
for q=1:1:size(dN_dRST,1);
    dN_dX(q,:)=(dX_dRST'\dN_dRST(q,:)')';
end

%% DERIVE THE DEFORMATION GRADIENT TENSOR
F_2D=zeros(2,2);
for q=1:1:size(x,1)
    F_2D=F_2D+(x(q,:)'*dN_dX(q,:));
end
F=eye(3,3);
F(1:2,1:2)=F_2D;

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
