function [varargout]=coarsepatch(F,V,f)

% function [F,V,Eb]=coarsepatch(F,V,f)
%-------------------------------------------------------------------------
% 
% 
% 
% Change log: 
% 2019/06/25 Added variable output handling as well as outputting of
% boundary edges
%-------------------------------------------------------------------------
%%
% Run reducepatch
[F,V] = reducepatch(F,V,f,'fast');

%%
% Avoid/fix reducepatch bugs
[F,V]=mergeVertices(F,V); %Merge vertices
F=patchRemoveCollapsed(F); %remove collapsed (edges)
[F,V]=triSurfRemoveThreeConnect(F,V); 
[F,V]=patchCleanUnused(F,V); %Remove unused points
F=uniqueIntegerRow(F); %Removing double faces

%%

varargout{1}=F;
varargout{2}=V;

if nargout==3
    %Also compute boundary
    varargout{3}=patchBoundary(F,V);
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
