%% gtitle
% Below is a demonstration of the features of the |gtitle| function

%%
clear; close all; clc;

%% Syntax
% |hf=gtitle(textString,);|

%% Description 
% The |gtitle| function can be used to create a title at the top of a
% figure window. The title is tied to the figure rather than any axes and
% can therefore be combined with figures with multiple subplots. 

%% Examples 

%%
% Settings 
titleString='Title on top'; %The title string
fontSize=25; %The desired font size

%%
% Visualize something with a title
hf=cFigure; 

hText=gtitle(titleString,fontSize,hf); %Create the title

surf(peaks(25)); 
colormap(gjet(250));
axisGeom; 
camlight headlight; 
drawnow; 

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
