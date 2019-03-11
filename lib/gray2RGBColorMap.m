function [C]=gray2RGBColorMap(varargin)

% function [C]=gray2RGBColorMap(G,cMap,cLim)
% ------------------------------------------------------------------------
% This function maps the monochrome data vector G, containing for instance 
% gray scale values, to a colormapped RGB array C using the colormap cMap
% and the color limits cLim. If C is not provided the jet colormap is used.
% If cLim is not provided or empty the colorlimits are adjusted according
% to the extreme values in G. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        G=varargin{1};
        cMap=jet(250);
        cLim=[];
    case 2
        G=varargin{1};
        cMap=varargin{2};
        cLim=[];
    case 3
        G=varargin{1};
        cMap=varargin{2};
        cLim=varargin{3};
end

if isempty(cLim)
   cLim=[min(G(:)) max(G(:))]; 
end

G=G(:); %G in column form

%% Create colormapped RGB array

%Snap to color limits
G(G<cLim(1))=cLim(1);
G(G>cLim(2))=cLim(2);

%Normalized to range [0-1]
G=G-min(G(:)); 
max_G=max(G(:));
if max_G>0
    G=G./max(G(:));
end
        
numColorLevels=size(cMap,1); %Number of color levels
linearIndex_G=1+round(G.*(numColorLevels-1));
C=cMap(linearIndex_G,:);

 
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
