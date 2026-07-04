%% cunique
% Below is a demonstration of the features of the |cunique| function

%%
clear; close all; clc;

%% Syntax
% |[A_uni,ind1,ind2,Ac]=cunique(A);|

%% Description 
% The |imx| function provides a figure window based GUI for 3D image
% segmentation

%% Examples

%%
% Plot settings
fontSize=20; 

%% Example 1: Getting unique entries and occurance counts for 1xN arrays

n=15;
A=round(25*rand(1,n)); %Rounded random set in range 0-25
A(1)=A(end); %Force at least one double occurance for this example
A

%Get unique set and counts
[A_uni,ind1,ind2,Ac]=cunique(A)

%% Example 2: Getting unique entries and occurance counts for NxM arrays

n=5;
m=6;
A=round(25*rand(n,m)); %Rounded random set in range 0-25
A(1)=A(end); %Force at least one double occurance for this example
A

%Get unique set and counts
[A_uni,ind1,ind2,Ac]=cunique(A)

%%
% Visualizing input array and occurange counts

cFigure; 
subplot(1,2,1); 
title('The input array')
hold on;
imagesc(A);
image_numeric(A,[],0,fontSize);
axis tight; axis equal; 
colormap(gca,gjet(max(A(:))));
icolorbar; 

subplot(1,2,2); 
title('The occurance counts')
hold on;
imagesc(Ac);
image_numeric(Ac,[],0,fontSize);
axis tight; axis equal; 
colormap(gca,gjet(max(Ac(:))));
icolorbar; 
drawnow;

%% Example 3: Getting unique entries and occurance counts for NxMx... arrays

n=3;
m=4;
l=2;

A=round(25*rand(n,m,l)); %Rounded random set in range 0-25
A(1)=A(end); %Force at least one double occurance for this example
A

%Get unique set and counts
[A_uni,ind1,ind2,Ac]=cunique(A)

%% Example 4: Using 'rows' option

n=5;
m=3;

A=round(25*rand(n,m)); %Rounded random set in range 0-25
A(1,:)=A(end,:); %Force at least one double row for this example
A

%Get unique set and counts
[A_uni,ind1,ind2,Ac]=cunique(A,'rows')

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
