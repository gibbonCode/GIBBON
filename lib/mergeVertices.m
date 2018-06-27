function [varargout]=mergeVertices(varargin)

% function [F,V,ind1,ind2]=mergeVertices(F,V,numDigitsMerge)

%% Parse input

switch nargin    
    case 2
        F=varargin{1};
        V=varargin{2};
        numDigitsMerge=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        numDigitsMerge=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

if isempty(numDigitsMerge)    
    D=patchEdgeLengths(F,V);    
    numDigitsMerge=6-numOrder(mean(D));
end

%% Merge nodes

[~,ind1,ind2]=unique(pround(V,numDigitsMerge),'rows');
V=V(ind1,:);

if size(F,1)==1
    F=ind2(F)'; %Fix indices in F
else
    F=ind2(F); %Fix indices in F
end

%% Collect output

varargout{1}=F;
varargout{2}=V;
varargout{3}=ind1;
varargout{4}=ind2;

%%
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
