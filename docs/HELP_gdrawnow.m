%% gdrawnow
% Below is a demonstration of the features of the |gdrawnow| function

%% Syntax
% |gdrawnow;|

%% Description
% The |gdrawnow| function is similar to the |drawnow| command but also
% activates the |vcw| widget if present for the current figure window.
% Note that activation of |vcw| also hides the axis interactive toolbars. 
% 
% See also: |drawnow|

%% Examples
clear; close all; clc; 

%% Calling |gdrawnow| equivalent to |drawnow|

%% 
% Some example data
[X,Y,Z]=peaks(25);

%%
figure; %Opens default MATLAB figure without vcw button
surf(X,Y,Z); %Visualize something in an axis
axisGeom; %Set axis options for geometry viewing
gdrawnow; %Fully equivalent to drawnow

%% Calling |gdrawnow| to drawnow and also activate existing |vcw| button

%%

cFigure; %Opens a cFigure window which contains the vcw button
surf(X,Y,Z); %Visualize something in an axis
axisGeom; %Set axis options for geometry viewing
gdrawnow; %drawnow + vcw activation

%% Calling |gdrawnow| to drawnow and activate |vcw| repeatedly

%%
% This example shows repeated calls to |gdrawnow| for a figure window e.g.
% when creating subplots. 

cFigure;
for i=1:1:2
    for j=1:1:2
        q=sub2ind([2,2],i,j);
        subplot(2,2,q); hold on;
        surf(X,Y,Z);
        colorbar;
        axisGeom;        
        gdrawnow; %drawnow + vcw activation
    end
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
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
