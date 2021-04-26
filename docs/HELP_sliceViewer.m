%% sliceViewer
% Below is a demonstration of the features of the |sliceViewer| function

%%
clear; close all; clc;

%% Syntax
% |hf=sliceViewer(varargin);|

%% Description 
% The |sliceViewer| function provides a figure window based GUI for 3D image
% visualization

%% Examples 
%
%% Example: Visualizing MRI data

%% 
% Example image data
load mri;
M=squeeze(D); %example image data set
v=2./[1,1,.4]; %example voxel size

%%
% Visualise using viewerType=1 featuring clickable navigation in mutually
% orthogonal 2D views
viewerType=1; %1 is default
sliceViewer(M,v,viewerType);

%%
% Visualise using viewerType=2 featuring a single 3D view and GUI sliders
% for navigation and thresholding
viewerType=2; 
sliceViewer(M,v,viewerType);

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
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
