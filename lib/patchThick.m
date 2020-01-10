function [varargout]=patchThick(varargin)

% function [E,V,Fp1,Fp2]=patchThick(Fp1,Vp1,dirSet,layerThickness,numSteps)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 
% 2019/04/18 Created based on quadThick
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        Fp1=varargin{1};
        Vp1=varargin{2};
        dirSet=1;
        layerThickness=mean(patchEdgeLengths(Fp1,Vp1));
        numSteps=1;
    case 3
        Fp1=varargin{1};
        Vp1=varargin{2};
        dirSet=varargin{3};
        layerThickness=mean(patchEdgeLengths(Fp1,Vp1));
        numSteps=1;
    case 4
        Fp1=varargin{1};
        Vp1=varargin{2};
        dirSet=varargin{3};
        layerThickness=varargin{4};
        numSteps=1;
    case 5
        Fp1=varargin{1};
        Vp1=varargin{2};
        dirSet=varargin{3};
        layerThickness=varargin{4};
        numSteps=varargin{5};
end

if size(dirSet,2)==3 %single vector or vector set for offsetting
    if size(dirSet,1)==1 %If one vector is given copy for all points
        Nv=dirSet(ones(size(Vp1,1),1),:);
    else %assume vectors are provided for all points
        Nv=dirSet;
        if size(dirSet,1)~=size(Vp1,1)
            error('Number of normals should match number of points');
        end
    end
else
    %Get vertex normals
    [~,~,Nv]=patchNormal(Fp1,Vp1);
    
    %Flip if desired
    Nv=Nv.*dirSet;
end

%% Compute coordinates

%Thicken inwards
Vp2=Vp1+layerThickness.*Nv;

%Get coordinates
X=linspacen(Vp1(:,1),Vp2(:,1),numSteps+1);
Y=linspacen(Vp1(:,2),Vp2(:,2),numSteps+1);
Z=linspacen(Vp1(:,3),Vp2(:,3),numSteps+1);

%Collect node set
V=[X(:) Y(:) Z(:)];

%% Create element sets 

if isa(Fp1,'cell')
    E=Fp1;
    Fp1n=Fp1;
    Fp2n=Fp1;
    for q=1:1:numel(Fp1)
        [E{q},Fp1n{q},Fp2n{q}]=createElementSet(Fp1{q},Vp1,numSteps);
    end
else     
    [E,Fp1n,Fp2n]=createElementSet(Fp1,Vp1,numSteps);    
end

%% Collect output
varargout{1}=E;
varargout{2}=V;
varargout{3}=Fp1n;
varargout{4}=Fp2n;
 
end

function [E,Fp1,Fp2]=createElementSet(Fp1,Vp1,numSteps)

%Create element matrix
E=repmat(Fp1,[numSteps,2]);
E_add=0:size(Vp1,1):size(Vp1,1)*(numSteps-1);
E_add=E_add(ones(size(Fp1,1),1),:);
E_add=E_add(:);
E_add=E_add(:,ones(size(Fp1,2),1));
E_add=[E_add E_add+size(Vp1,1)];
E=E+E_add;

%Create top and bottom face set
Fp1=E(1:size(Fp1,1),1:size(Fp1,2));
Fp2=E(1+(end-size(Fp1,1)):end,size(Fp1,2)+1:end);

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
