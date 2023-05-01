function [varargout]=patchEdgeCollapse(varargin)

% function [F,V,logicValid,indFix2]=patchEdgeCollapse(F,V,E,logicKeep,meanOption)
%-------------------------------------------------------------------------
%
%
%-------------------------------------------------------------------------

%% Parse input
switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        logicKeep=[];
        meanOption=1;
    case 4
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        logicKeep=varargin{4};
        meanOption=1;
    case 5
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        logicKeep=varargin{4};
        meanOption=varargin{5};
end

if isempty(logicKeep)
   logicKeep=false(size(E));
   logicKeep(:,1)=1;
end
%% 
%

numEdges=size(E,1);
for q=1:1:numEdges
    edgeCollapseNow=E(q,:); %Current edge to collaps
    logicKeepNow=logicKeep(q,:);
    if all(logicKeepNow==0) || all(logicKeepNow==1)
        error('Invalid logic for keeping edge points provided, cannot remove or keep both points');
    end
    if nnz(ismember(edgeCollapseNow,F))==2 && (edgeCollapseNow(1)~=edgeCollapseNow(2))
        
        F(F==edgeCollapseNow(logicKeepNow==0))=edgeCollapseNow(logicKeepNow==1); %Replace index in face array
        E(E==edgeCollapseNow(logicKeepNow==0))=edgeCollapseNow(logicKeepNow==1); %Replace index in edge array
        if meanOption==1
            V(edgeCollapseNow(1),:)=(V(edgeCollapseNow(1),:)+V(edgeCollapseNow(2),:))/2; %Replace coordinate by mean of the two edge points
        end
        
    end
end

logicValid=~any(diff(sort(F,2),[],2)==0,2);
F=F(logicValid,:); %Keep only valid faces

%Remote unused points and update faces matrix
[F,V,indFix2]=patchCleanUnused(F,V);

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=logicValid;
varargout{4}=indFix2;

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
