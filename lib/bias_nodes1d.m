function [x]=bias_nodes1d(x,f)

% function [x]=bias_nodes1d(x,f)
% ------------------------------------------------------------------------
% This function biases the spacing for the entries in x using the factor f
% such that the output array goes from x(1) to x(end) but biased using x^f
% (i.e. the spacing decreasing depending on the magnitude of x). 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

min_x=min(x,[],2)*ones(1,size(x,2)); 
x=x-min_x; 
max_x=max(x,[],2)*ones(1,size(x,2)); 
x=max_x.*((x.^f)./(max_x.^f)); 
x=x+min_x;
 
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
