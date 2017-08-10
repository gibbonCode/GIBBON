%% vcw
% Below is a demonstration of the features of the |vcw| function

%% Syntax
% |vcw(hf,buttonOpt);|

%% Description
% The vcw function, the view control widget, allows the user to rotate, pan
% and zoom a figure using key presses and mouse gestures. This makes 3D
% view manipulation feel like what one would expect of a CAD package. 
%
% The vcw widget is loaded by default when using cFigure but can be loaded
% for any other figure by using: 
% figure; vcw;
% You can load vcw either before or after creating axes. 
%
% This code was inspired by the fcw function by Oliver Woodford (which was
% based on Torsten Vogel's view3d function, which was in turn inspired by
% rotate3d from The MathWorks, Inc.).
% VCW is inspired by and similar to fcw function by Oliver Woodford.
% However vcw offers the following added features (including toolbar
% button):  
%%
%
% # handing of colorbars (bug in fcw when view(2) is used combined with
% panning which induced zooming and panning)   
% # overobj based axes selection so that the current axes is determined
% based on mouse pointer location for most functions   
% # A toggle button for activation and deactivation in the figure toolbar 
% # ability to start vcw before objects are plotted 
% # "proper" closure of the vcw widget, in fcw the q button did not exit
% the keyDown functions such as panning etc. Now the quit action 
% deactivates the widget   
% # Upon activation of the vcw widget the plotting and default view
% manipulation tools and buttons are disabled 
% (to avoid interference with vcw) 
% # Added "linked" mode by using ALT button to alter views for all axes in
% figure uppon keypress  
% # Altered keypress functions and behaviour with SHIFT (which now negates
% directions)  
% # Added i key option to display help information. 
% # Allows for saving of a new default view

%% 
% Once loaded users may press v key to activate widget (or press the
% toolbar button). Press i to show the following help information:
% 
% <<vcw_help_show.png>>

%% Examples

%% Using the preloaded view control widget (vcw) in |cFigure|
% A |cFigure| window with the view control widget loaded by default. To
% activate widget press v or click on vcw button in the toolbar. 
% 
% <<vcw_icon_menubar.png>>

%% Use with |figure| 
% The default MATLAB |figure| does not contain the view control widget. To
% load it here the user must enter |vcw;| either before or after axes are
% created. e.g.: 
%%
% |figure; surf(peaks(25)); axis equal; axis tight; vcw;|
%%
% which is equivalent to: 
%%
% |figure; vcw; surf(peaks(25)); axis equal; axis tight;| 

%% Tip
% You can create your own figure function that simply contains |figure;
% vcw;| to create a standard MATLAB figure containing the view control
% widget. 

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
