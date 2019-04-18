function [varargout]=quadThick(varargin)

% function [E,V,Fq1,Fq2]=quadThick(Fq,Vq,dirSet,layerThickness,numSteps)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 
% 2017/08/31 Added function description to top of function
% 2019/04/18 Updated to use patchThick (generalized form of old quadThick)
% instead 
% ------------------------------------------------------------------------

%%

[E,V,Fp1,Fp2]=patchThick(varargin{:});
 
%Collect output
varargout{1}=E;
varargout{2}=V;
varargout{3}=Fp1;
varargout{4}=Fp2;

%%
% Old version: 
% switch nargin
%     case 2
%         Fq=varargin{1};
%         Vq=varargin{2};
%         dirSet=1;
%         layerThickness=mean(patchEdgeLengths(Fq,Vq));
%         numSteps=1;
%     case 3
%         Fq=varargin{1};
%         Vq=varargin{2};
%         dirSet=varargin{3};
%         layerThickness=mean(patchEdgeLengths(Fq,Vq));
%         numSteps=1;
%     case 4
%         Fq=varargin{1};
%         Vq=varargin{2};
%         dirSet=varargin{3};
%         layerThickness=varargin{4};
%         numSteps=1;
%     case 5
%         Fq=varargin{1};
%         Vq=varargin{2};
%         dirSet=varargin{3};
%         layerThickness=varargin{4};
%         numSteps=varargin{5};
% end
% 
% %%
% 
% if size(dirSet,2)==3 %single vector or vector set for offsetting
%     if size(dirSet,1)==1 %If one vector is given copy for all points
%         Nv=dirSet(ones(size(Vq,1),1),:);
%     else %assume vectors are provided for all points
%         if size(dirSet,1)~=size(Vq,1)
%             error('Number of normals should match number of points');
%         end
%     end
% else
%     %Get vertex normals
%     [~,~,Nv]=patchNormal(Fq,Vq);
%     
%     %Flip if desired
%     Nv=Nv.*dirSet;
% end
% 
% %Thicken inwards
% Vq2=Vq+layerThickness.*Nv;
% 
% %Get coordinates
% X=linspacen(Vq(:,1),Vq2(:,1),numSteps+1);
% Y=linspacen(Vq(:,2),Vq2(:,2),numSteps+1);
% Z=linspacen(Vq(:,3),Vq2(:,3),numSteps+1);
% 
% %Collect node set
% V=[X(:) Y(:) Z(:)];
% 
% %Create element matrix
% E=repmat(Fq,[numSteps,2]);
% E_add=0:size(Vq,1):size(Vq,1)*(numSteps-1); 
% E_add=E_add(ones(size(Fq,1),1),:);
% E_add=E_add(:);
% E_add=E_add(:,ones(4,1));
% E_add=[E_add E_add+size(Vq,1)];
% E=E+E_add;
% 
% %Create boundary face set
% Fq1=E(1:size(Fq,1),1:4);
% Fq2=E(1+(end-size(Fq,1)):end,5:end);
% 
% %Collect output
% varargout{1}=E;
% varargout{2}=V;
% varargout{3}=Fq1;
% varargout{4}=Fq2;
 
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
