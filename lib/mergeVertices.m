function [varargout]=mergeVertices(varargin)

% function [F,V,ind1,ind2]=mergeVertices(F,V,nKeep)

%% Parse input

switch nargin    
    case 2
        F=varargin{1};
        V=varargin{2};
        nKeep=5;
    case 3
        F=varargin{1};
        V=varargin{2};
        nKeep=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

%% Merge nodes

[~,ind1,ind2]=unique(pround(V,nKeep),'rows');
V=V(ind1,:);
F=ind2(F);

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
