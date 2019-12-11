function [cMap]=kvw(varargin)

% function [cMap]=kvw(n)
% ------------------------------------------------------------------------
% Creates the colormap data for the black-viridis-white colormap. This
% colormap start black, continues as the viridis colormap and ends white. 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2019/06/27
%------------------------------------------------------------------------

%%
switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%%
ns=75;
v=viridis(250);
k=linspacen([0 0 0],v(1,:),ns)';
w=linspacen(v(end,:),[1 1 1],ns)';

cMap=[k(1:end-1,:);v;w(2:end,:)];
cMap=resampleColormap(cMap,n);

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
