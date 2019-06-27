function [cMap]=iwarmcold(varargin)

% function [cMap]=iwarmcold(n)
% ------------------------------------------------------------------------
% Creates the colormap data for n levels for the  inverse warm-cold
% colormap. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/03/21
%------------------------------------------------------------------------

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cf=[213 15  37; 238 178 17;]./255;
c=rot90(cf,2);
h=cf;

h=resampleColormap(h,3);
c=resampleColormap(c,3);

cMap=[c; [0 0 0]; h];

[cMap]=resampleColormap(cMap,n);
 
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
