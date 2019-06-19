function [varargout]=tri2rhombi(varargin)

% function [FQ,VQ,CQ]=tri2rhombi(F,V,C)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        C=[];
    case 3        
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
end

%%

conStruct=patchConnectivity(F,V,{'ef','ev'});
EF=conStruct.edge.face;
E=conStruct.edge.vertex;

logicBoundary=any(EF==0,2);
indF=EF(logicBoundary,1);

VF=patchCentre(F,V);
VE=patchCentre(E,V);
VC=VF;
VC(indF,:)=VE(logicBoundary,:);

e=E(~logicBoundary,:);
ef=EF(~logicBoundary,:);

VQ=[V;VC];
FQ=[e(:,1) ef(:,1)+size(V,1) e(:,2) ef(:,2)+size(V,1)];

%Derive color data for refined set
if ~isempty(C) && nargout>2
    if size(C,1)==size(F,1) %Face color data        
        CQ=(C(ef(:,1),:)+C(ef(:,2),:))./2;
    elseif size(C,1)==size(V,1) %Vertex color data
        CC=vertexToFaceMeasure(F,C);
        CQ=[C;CC];
    else 
        error('Color data should be nxq in size whereby n is the number of faces or the number of vertices');
    end
else
    CQ=[];
end

%% Collect output
varargout{1}=FQ; %Faces
varargout{2}=VQ; %Vertices
varargout{3}=CQ; %New color data for faces or vertices

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
