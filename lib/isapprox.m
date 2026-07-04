function [L]=isapprox(varargin)

% function [L]=isapprox(A,B,tolLevel)
% ------------------------------------------------------------------------
% This function returns if A is approximately equal to B (to within
% tolLevel) in the form of the logical L. 
%
% Change log 
% 2023/08/31 KMM: Updated description/documentation 
% ------------------------------------------------------------------------
%%

switch nargin
    case 1
        A=varargin{1};
        B=zeros(size(A));
        tolLevel=eps(1);
    case 2
        A=varargin{1};
        B=varargin{2};
        tolLevel=eps(1);
    case 3
        A=varargin{1};
        B=varargin{2};
        tolLevel=varargin{3};
end

if numel(B)==1
    B=B.*ones(size(A));
end

L=abs(A-B)<tolLevel;

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
