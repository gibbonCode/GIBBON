function d=divisors(varargin)

% function d=divisors(x,arg)
% -----------------------------------------------------------------------
% Determine divisors of input x (which needs to be natural number)
% By default the output provides a vector containing all positive divisors
% but this can be altered depending on arg. 
%
% Change log: 
% 2019/06/27 Added variable input handling, improved error handling
% -----------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        x=varargin{1};
        outputType='pos';
    case 2
        x=varargin{1};
        outputType=varargin{2};
end

%%

x=abs(x); %Force absolute
X=0:1:x; %The candidate range [0,x]

d=X(ismember(X,x./X)); %Keep entries which apear as integers after division

%Check if negative output is requested
switch outputType
    case 'pos'
        
    case 'neg'
        d=-fliplr(d);
    case 'all'
        d=[-fliplr(d) d];
    otherwise
        error('Incorrect outputType provided, valid options are pos, neg or all');
end

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
