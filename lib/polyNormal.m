function [varargout]=polyNormal(V_poly)

% [N,Vn,Nv]=polyNormal(V_poly)
% ------------------------------------------------------------------------
% Normals directions are derived based on cross product of curve segments
% with the outware z-direction. A horizontal and clockwise curve segment
% therefore results in an upward normal. The normals are provided for each
% segment or for each point on the curve. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2017/03/16 Created
%------------------------------------------------------------------------

%%

zDir=[0 0 1];

if size(V_poly,2)==2 %Cope with 2D
    V_poly(:,3)=0;
end

[Vn,U]=polyDerivative(V_poly,1);
N=cross(zDir(ones(size(U,1),1),:),U);
N=vecnormalize(N);

varargout{1}=N(:,[1 2]); 
varargout{2}=Vn(:,[1 2]); 

if nargout==3
    [~,U]=polyDerivative(V_poly,3);
    Nv=cross(zDir(ones(size(U,1),1),:),U);
    Nv=vecnormalize(Nv);
    varargout{3}=Nv(:,[1 2]);
end

end

function [V,U]=polyDerivative(Vg,dirOpt)

switch dirOpt
    case 1 %Forward
        V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
        U=diff(Vg,1,1);        
    case 2 %Backward
        V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
        U=flipud(diff(flipud(Vg),1,1));        
    case 3 %Central
        V=Vg;
        Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
        Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
        U=Uf;
        U(:,:,2)=Ub;
        U=gnanmean(U,3);        
end
U=vecnormalize(U);

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
