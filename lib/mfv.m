function mfv(varargin)

% function mfv(H)
% %Makes figure or figures visible depending on the input handles. If none
% are provide all are searched and made visible. 


%% Parse input

switch nargin
    case 0 %No handles provided so assume all
        H=findall(0,'type','figure');
    case 1
        H=varargin{1};
    otherwise
        error('Wrong number of input arguments');
end

%% Make figure(s) visible

for q=1:1:numel(H)
    h=H(q);    
    if verLessThan('matlab', '8.4.0.150421 (R2014b)')
        set(h,'Visible','On');
    else 
        h.Visible='On';
    end
    drawnow;
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
