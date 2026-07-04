%% patchTorus
% Below is a demonstration of the features of the |patchTorus| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=patchTorus(r,nr,rc,nc,patchType);|

%% Description
% This function generates patch for a torus with desired dimensions and
% mesh type

%% Examples

%% 
% Plot settings
fontSize=10;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=2;
markerSize=25;

%% Building a mesh of a torus

%Torus parameters
r=1; %Sphere radius
rc=2.5; %Central radius
nr=16;
nc=25;
patchTypes={'quad','tri_q','tri','honey'};

% Open figure for plotting
cFigure; 

%Plot the various mesh types
pColors=gjet(4);
for q=1:1:4
    [F,V]=patchTorus(r,nr,rc,nc,patchTypes{q});
    subplot(2,2,q); hold on;
    title([patchTypes{q}],'FontSize',fontSize,'Interpreter','none');
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hp=patch('Faces',F,'Vertices',V);
%     [hn]=patchNormPlot(F,V,0.3);
    set(hp,'FaceColor',pColors(q,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    camlight headlight; 
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  grid on;
end
drawnow;

%% Simulate a doughnut

%Torus parameters
r=1; %Sphere radius
rc=2.5; %Central radius
n=3; 
nr=64*n;
nc=100*n;
patchType='tri';

[F,V]=patchTorus(r,nr,rc,nc,patchType);
C_bread=[254 191 78]./255;
C_red=pColors(end,:);
C_blue=pColors(1,:);
C=C_bread(ones(size(V,1),1),:);

X=V(:,3);
L=X>mean(X);
C(L,:)=ones(nnz(L),3);

X=V(:,1);
Lr=X>-2 & X<-1.8 & L;
C(Lr,:)=C_red(ones(nnz(Lr),1),:);

Lr=X>-0.1 & X<0.1 & L;
C(Lr,:)=C_red(ones(nnz(Lr),1),:);

Lr=X>1.8 & X<2 & L;
C(Lr,:)=C_red(ones(nnz(Lr),1),:);

d=0.25;
Lr=X>-2-d & X<-1.8-d & L;
C(Lr,:)=C_blue(ones(nnz(Lr),1),:);

Lr=X>-0.1-d & X<0.1-d & L;
C(Lr,:)=C_blue(ones(nnz(Lr),1),:);

Lr=X>1.8-d & X<2-d & L;
C(Lr,:)=C_blue(ones(nnz(Lr),1),:);

%%
% Visualize a doughnut
cFigure; hold on;
title('doughnut','FontSize',fontSize); 
hp=gpatch(F,V,C,'none');
hp.FaceColor='interp';
axisGeom(gca,fontSize);
camlight headlight; lighting phong;
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
