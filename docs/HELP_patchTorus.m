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
np=100;
r=1; %Sphere radius
rc=2.5; %Central radius
nr=16;
nc=25;
patchTypes={'quad','tri_q','tri','honey'};

% Open figure for plotting
cFigure; 

%Plot the various mesh types
pColors=gjet(4);
for q=1:1:4;
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

%% Simulate your lunch

%Torus parameters
np=100;
r=1; %Sphere radius
rc=2.5; %Central radius
n=2; 
nr=64*n;
nc=100*n;
patchType='tri';

[F,V]=patchTorus(r,nr,rc,nc,patchType);
C_bread=[254 191 78]./255;
C_red=pColors(end,:);
C_blue=pColors(1,:);
C=C_bread(ones(size(F,1),1),:);

X=V(:,3);
XF=mean(X(F),2);
L=XF>mean(XF);
C(L,:)=ones(nnz(L),3);

X=V(:,1);
XF=mean(X(F),2);
Lr=XF>-2 & XF<-1.8 & L;
C(Lr,:)=C_red(ones(nnz(Lr),1),:);

Lr=XF>-0.1 & XF<0.1 & L;
C(Lr,:)=C_red(ones(nnz(Lr),1),:);

Lr=XF>1.8 & XF<2 & L;
C(Lr,:)=C_red(ones(nnz(Lr),1),:);

d=0.25;
X=V(:,1);
XF=mean(X(F),2);
Lr=XF>-2-d & XF<-1.8-d & L;
C(Lr,:)=C_blue(ones(nnz(Lr),1),:);

Lr=XF>-0.1-d & XF<0.1-d & L;
C(Lr,:)=C_blue(ones(nnz(Lr),1),:);

Lr=XF>1.8-d & XF<2-d & L;
C(Lr,:)=C_blue(ones(nnz(Lr),1),:);

cFigure; hold on;
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hp=patch('Faces',F,'Vertices',V);
set(hp,'FaceColor','flat','FaceVertexCData',C,'edgeColor','none');
camlight headlight; lighting phong;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on; axis off; 
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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
