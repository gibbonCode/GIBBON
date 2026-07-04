function [pp,t]=cs3dPath(V,p,w)

% function [pp,t]=cs3dPath(V,p,w)
% ------------------------------------------------------------------------
% This function generates the pp-form for a cubic-smoothing spline to the
% 2D or 3D curve defined by V. 
% 
% Equivalent to CSCVN performance if p=1 and w=ones(size(V,1),1)
%
% See also: |cscvn|, |csaps|
% ------------------------------------------------------------------------

%%

dt = sum((diff(V,[],1).^2),2)'; %Point spacing measure
t = cumsum([0,dt.^(1/4)]); %Curve length measure

pp = csaps(t',V',p,[],w); %Smoothened ppform
 
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