function [varargout]=scatterV(V,varargin)

% function [hp]=scatterV(V,varargin)
% ------------------------------------------------------------------------
% This function is similar to scatter3 except that the coordinate set can
% be specified in a single nx3 array. In fact the function simply does: 
%
%       hp=scatter3(V(:,1),V(:,2),V(:,3),varargin{:});
%
% Kevin Mattheus Moerman
%
% ------------------------------------------------------------------------

%% 
if ~isempty(V)
    nDims=size(V,2);
    
    %Add zeros if input is 2D
    if nDims==2
        V(:,3)=0;
    end
    
    hp=scatter3(V(:,1),V(:,2),V(:,3),varargin{:});
    
else
    hp=[];
end

switch nargout
    case 1
        varargout{1}=hp;
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
