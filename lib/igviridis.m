function [cMap]=igviridis(varargin)

% function [cMap]=igviridis(n)
% ------------------------------------------------------------------------
% Creates the colormap data for the inver gray-viridis colormap. This
% colormap combines gray scale colors with the viridis colors. The colormap
% is circular in that it starts and ends with white. 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log: 
% 2019/03/28
% 2019/06/27 Increased color resolution
%------------------------------------------------------------------------

%%
switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

%%

Ng=250;
Nw=75;
Nk=50;
Nv=Ng-Nk-Nw;

g=linspacen([0 0 0],[1 1 1],Ng)';
v=flipud(viridis(Nv));
k=linspacen(v(end,:),[0 0 0],Nk)';
w=linspacen([1 1 1],v(1,:),Nw)';

cMap=[g(1:end-1,:); [1 1 1]; w; v; k];
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
