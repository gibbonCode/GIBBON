function [varargout]=patchAnnotate(F,V,d,varargin)

% function [hv,hf]=patchAnnotate(F,V,d,varargin)
% ------------------------------------------------------------------------
%
% Change log: 
% 2019/08/06 Created
% 
% ------------------------------------------------------------------------

%% Parse input

% Force 3D coordinates
if size(V,2)==2
    V(:,3)=0; 
end

%% Get node indices 

nodeIndices=1:1:size(V,1);
faceIndices=1:1:size(F,1);

%% Get coordinates for text segments

if isempty(d) 
    d=mean(patchEdgeLengths(F,V))/10; %
end

[NF,VF,NV]=patchNormal(F,V);
VV=V+NV.*d;
VF=VF+NF.*d;

%% Annotate nodes and faces

hv=pointAnnotate(VV,nodeIndices,varargin{:},'FontWeight','Normal','FontAngle','italic'); %Annotate Nodes
hf=pointAnnotate(VF,faceIndices,varargin{:},'FontWeight','Bold'); %Annotate faces

%% Collect output

varargout{1}=hv;
varargout{2}=hf;

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
