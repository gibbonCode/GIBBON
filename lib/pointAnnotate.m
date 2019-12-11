function [varargout]=pointAnnotate(V,nodeIndices,varargin)

% function [ht]=pointAnnotate(V,nodeIndices,varargin)
% ------------------------------------------------------------------------
% This function annotates the points definced by V using the indices
% nodeIndices. Optional additional inputs are those associated with
% MATLAB's text function. 
% 
% Change log: 
% 2019/08/06 Created
% ------------------------------------------------------------------------

%% Get node indices 
if isempty(nodeIndices)
    nodeIndices=1:1:size(V,1);
end

%% Create text data 
t = sprintfc('%i',nodeIndices); %Text cell for node id's

%% Plot text at coordinates
varargout{1}=text(V(:,1),V(:,2),V(:,3),t,varargin{:}); %Plot text at nodes

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
