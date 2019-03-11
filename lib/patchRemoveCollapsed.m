function [varargout]=patchRemoveCollapsed(F)

% function [F,logicKeep]=patchRemoveCollapsed(F)
% ------------------------------------------------------------------------
% 
% ------------------------------------------------------------------------

%% Detect repeated indices for all faces
% Create logic defining collapsed faces

F_sort=sort(F,2); %Sort faces in 2nd direction so 2 1 2 -> 1 2 2
d=diff(F_sort,[],2); %Difference in 2nd direction so 1 2 2 -> 1 0
logicKeep=~any(d==0,2); %Logic for faces without zeros in difference measure of sorted faces

%% Create a new face set
F=F(logicKeep,:); %Selecting faces without repeated indices

%% Collect output

varargout{1}=F;
varargout{2}=logicKeep;

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
