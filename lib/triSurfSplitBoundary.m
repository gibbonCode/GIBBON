function [varargout]=triSurfSplitBoundary(varargin)

% function [F,V,E,CF,CV]=triSurfSplitBoundary(F,V,E,n,CF,CV)
% ------------------------------------------------------------------------
%
% 
% To do: Use for loop and recursion instead
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=size(E,1)+1;
        CF=[];
        CV=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=varargin{4};
        CF=[];
        CV=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=varargin{4};
        CF=varargin{5};
        CV=[];
    case 6
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=varargin{4};
        CF=varargin{5};
        CV=varargin{6};
end

if isempty(CF)
    CF=ones(size(F,1),1);
end

%%
if isempty(E)
    error('Edge set is empty');
end

if size(E,1)==n
    %Do nothing: Input = output
elseif size(E,1)>n
    error('Current number of edges exceeds number of desired edges');
else
    while size(E,1)~=n
        [~,indMax]=max(edgeLengths(E,V));
        numNodesPre=size(V,1); %Number of nodes before split
        [F,V,~,CF,CV]=triEdgeSplit(F,V,E(indMax,:),CF,CV);
        numNodesPost=size(V,1); %Number of nodes after split
        [Eb]=patchBoundaryLabelEdges(F,V,CF); %All boundary edges
        logicKeep=all(ismember(Eb,E) |...
            ismember(Eb,(numNodesPre+1):numNodesPost),2);
        E=Eb(logicKeep,:);
    end
end

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=E;
varargout{4}=CF;
varargout{5}=CV;

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
