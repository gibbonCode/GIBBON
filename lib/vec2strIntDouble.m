function [t]=vec2strIntDouble(varargin)

%%
switch nargin
    case 1
        n=varargin{1};
        formatLong='%6.7e';
    case 2
        n=varargin{1};
        formatLong=varargin{2};
end

%%

if isnumeric(n) %If it is numeric
    n=double(n);
    if isrounded(n) %If it looks like an integer
        t_form=repmat('%d, ',1,numel(n)); 
    else %Not an integer
        t_form=repmat([formatLong,', '],1,numel(n)); 
    end
    t_form=t_form(1:end-2); %Take away last space and comma
    
    %Convert to string
    t=sprintf(t_form,n);
elseif char(n)
    t=n;
else
    error('Input should be numeric')
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
