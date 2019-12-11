function [E]=curveToEdgeList(N)
% function [E]=curveToEdgeList(N)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%%

if numel(N)==1 %the size of the list is specified
    indList=(1:1:N)';
elseif ismatrix(N) %ordered vertices are provided
    indList=(1:1:size(N,1))';
else %The indices are provided
    indList=N(:);
end

%%

E=[indList(1:end-1) indList(2:end)];


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
