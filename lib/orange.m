function [c]=orange(varargin)

% function [c]=orange(m)
% ------------------------------------------------------------------------
% This function outputs the rgb color values (between 0 and 1) for an
% orange color, i.e. [255 123 21]/255
% If the optional output m is provided it is used to copy the color m
% times, such than an mX3 array is returned. 
% ------------------------------------------------------------------------

%%

%Orange color
c=[255 123 21]/255;

%Replicate if needed
switch nargin 
    %case 0
    %c=[255 123 21]/255;
    case 1
        m=varargin{1};
        c=repmat(c,m,1);
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
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
