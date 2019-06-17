function str=strwrap(varargin)

% function str=strwrap(str,n,pattern)
%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        str=varargin{1};
        n=[];
        pattern=[];
    case 2
        str=varargin{1};
        n=varargin{2};
        pattern=[];
    case 3
        str=varargin{1};
        n=varargin{2};
        pattern=varargin{3};
end

if isempty(n)
    n=1;
end

if isempty(pattern)
    pattern=', ';
end

%%
try %Use count as introduced in R2016b
    N = count(str,pattern); 
catch %Use old approach
    N = numel(strfind(str,pattern));
end

rangeSteps=n:n:N;
for q=1:1:numel(rangeSteps)
    str=regexprep(str,pattern,'\n',rangeSteps(q)-(q-1));
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
