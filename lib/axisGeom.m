function axisGeom(varargin)

% function axisGeom(h,fontSize)
%-------------------------------------------------------------------------
% This function sets axis properties commonly used for viewing geometry in
% 3D, e.g. axis equal; axis vis3d; axis tight; The optional inputs include
% the axis handle and the desired font size. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2018/04/03 Fixed bug in relation to MATLAB R2017a-MATLAB R2018a where
% axes(h) triggers a drawnow event.  axes(h) was therefore suppressed. 
%------------------------------------------------------------------------

%%
% Parse input
switch nargin
    case 0
        h=gca;
        fontSize=15;
    case 1
        h=varargin{1};
        fontSize=15;
    case 2
        h=varargin{1};
        fontSize=varargin{2};
end

if isempty(h)
    h=gca;
end

%%

% axes(h); %Causes drawnow event in MATLAB R2017a-MATLAB R2018a

view(3); %Switch to 3D view

%Set axis labels
xlabel('X','FontSize',fontSize);
ylabel('Y','FontSize',fontSize); 
zlabel('Z','FontSize',fontSize);

%Set axis properties
axis equal; axis vis3d; axis tight; %Set axis properties for geometry viewing
grid on; %Turn on the grid
box on; %Turn on a box around axis
h.FontSize=fontSize; %Set font size
h.Clipping='off'; %Turn off clipping in relation to zooming

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
