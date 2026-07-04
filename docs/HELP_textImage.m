%% textImage
% Below is a demonstration of the features of the |textImage| function

%%
clear; close all; clc;

%% Syntax
% |[glyphImage]=textImage(textString,FontName,FontSize,padAmount);|

%% Description
% This creates image data containing text content with a desired size and
% font. 

%% Examples

%%
% Plot settings
cMap=gray(250); %Colormap

%% Example: Using |textImage| to create an image with text
% You can use |listTrueTypeFonts| to check what fonts are available on your
% machine

%%
% Settings
textString={'Text'}; %String to be rendered in image

%%
% Creating the text image
[M]=textImage(textString);

%%
% Visualizing the text image
cFigure;
imagesc(M);
axis equal; axis tight; grid off;  
colormap(cMap); caxis([min(M(:)) max(M(:))]); colorbar;
drawnow;

%% Example: Using full set of inputs
% You can use |listTrueTypeFonts| to check what fonts are available on your
% machine

%%
% Settings
textString={'Text'}; %String to be rendered in image
FontName='Arial'; %Font name 
FontSize=25; %Font height in pixels
padAmount=20; %Pixels padded around string

%%
% Creating the text image
[M]=textImage(textString,FontName,FontSize,padAmount);

%%
% Visualizing the text image
cFigure;
imagesc(M);
axis equal; axis tight; grid off;  
colormap(cMap); caxis([min(M(:)) max(M(:))]); colorbar;
drawnow;

%% Example: Use with cell inputs and text across multiple lines
% You can use |listTrueTypeFonts| to check what fonts are available on your
% machine

%%
% Settings
textString={'Lorem ipsum dolor sit amet',...
            'consectetur adipiscing elit, sed',...
            'do eiusmod tempor incididunt ut',... 
            'labore et dolore magna aliqua.'}; %String to be rendered in image
FontName='Arial'; %Font name 
FontSize=50; %Font height in pixels
padAmount=25; %Pixels padded around string

%%
% Creating the text image
[M]=textImage(textString,FontName,FontSize,padAmount);

%%
% Visualizing the text image
cFigure;
imagesc(M);
axis equal; axis tight; grid off;  
colormap(cMap); caxis([min(M(:)) max(M(:))]); colorbar;
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
