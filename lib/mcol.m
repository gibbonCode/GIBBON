function [b]=mcol(a)

% function [b]=mcol(a)
% ------------------------------------------------------------------------
% This function converts the input matrix (or vector) a to a column vector. 
%
%
% Change log: 
% 2023/03/10: Kevin Moerman Created
%
% ------------------------------------------------------------------------
%%

if ismatrix(a)
    if ~isvector(a) %Check if it is a vector already
        b=a(:);
    else %Vector so either row already or a column
        if ~iscolumn(a)
            b=a'; %Transpose
        else
            b=a; %Keep
        end
    end
else
    error('Input is not of matrix type. The mrow function is only defined for matrix arrays, see also ismatrix');
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
