%% MRcart2im
% Below is a demonstration of the features of the |MRcart2im| function

%%
clear; close all; clc;

%% Syntax
% |[I,J,K]=MRcart2im(Xc,Yc,Zc,v,OR,r,c);|

%% Description 

%% Examples 
% 

defaultFolder = fileparts(fileparts(mfilename('fullpath'))); %Set main folder
pathName=fullfile(defaultFolder,'data','DICOM','0001_human_calf');
loadName=fullfile(pathName,'IMDAT','IMDAT.mat');

IMDAT_struct=load(loadName); %The image data structure
G = IMDAT_struct.G; %Geometric/spatial information
OR = G.OR; %Origin location
v = G.v; %The voxel size
r = G.r; %Row direction
c = G.c; %Column direction
M = IMDAT_struct.type_1; %The image data
siz = IMDAT_struct.ImageSize;

%%


%%

L_plot=false(size(M));
L_plot(round(siz(1)/2),:,:)=1;
L_plot(:,round(siz(2)/2),:)=1;
L_plot(:,:,1)=1;
L_plot(:,:,end)=1;
% L_plot(:,:,round(siz(3)/2))=1;

[F,V_XYZ_reg,C]=im2patch(M,L_plot,'vb',v);

V_IJK=cart2im(V_XYZ_reg(:,[1 2 3]),v);

%%
I=V_IJK(:,1);
J=V_IJK(:,2);
K=V_IJK(:,3);

[Xc,Yc,Zc]=im2MRcart(I,J,K,v,OR,r,c);

V_XYZ_real=[Xc Yc Zc];

s=cross(c',r')'; %Determine s => slice direction vector

%%

[Ic,Jc,Kc]=MRcart2im(Xc,Yc,Zc,v,OR,r,c);

V_reg_check=im2cart(Ic,Jc,Kc,v);

%%

cFigure; 
subplot(1,2,1); hold on; 
title('Regular system (no tilt/shift)');
gpatch(F,V_XYZ_reg,C,'none'); 
plotV(V_reg_check,'r.','MarkerSize',5);
colormap gray; 
axisGeom; camlight headlight; 

subplot(1,2,2); hold on; 
title('Real world scanner coordinates');
gpatch(F,V_XYZ_real,C,'none'); 
plotV(OR([2 1 3])','y.','MarkerSize',50)
quiverVec(OR([2 1 3])',r([2 1 3])',50,'r');
quiverVec(OR([2 1 3])',c([2 1 3])',50,'g');
quiverVec(OR([2 1 3])',s([2 1 3])',50,'b');
colormap gray; 
axisGeom; camlight headlight; 

gdrawnow; 

%%


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
