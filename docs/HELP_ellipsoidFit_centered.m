%% HELP_ellipsoidFit_centered
% Below is a demonstration of the features of the |ellipsoidFit_centered| function

%% Syntax
% |[M,ellipStretch,R,MU]=ellipsoidFit_centered(X,MU);|

%% Description 
% The |ellipsoidFit_centered| function fits an ellipsoid to data when the
% ellipsoid centre is known. If the centre is not provided the mean of the
% input point set will be assumed to be the centre. 

%% Examples

%% 
clear; close all; clc;

%%
% Plot settings

figColor='w';
figColorDef='white';
fontSize=11;

%% Example: Using |ellipsoidFit_centered| to fit an ellipsoid to a point cloud with known centre

%%
% Simulating an ellipsoid with known directions

% Ellipsoid axis stretch factors
ellipStretchTrue=[pi 2 0.5];
MU_true=[1 6 pi];

% Create ellipsoid patch data
[F,X,~]=geoSphere(3,1);
x=X(:,1); 
FX=mean(x(F),2);
logicKeep=FX>0;
F=F(logicKeep,:);
indKeep=unique(F(:));
indFix=nan(size(X,1),1);
indFix(indKeep)=1:numel(indKeep);
X=X(indKeep,:);
F=indFix(F); 
X=X.*ellipStretchTrue(ones(size(X,1),1),:);

%Create Euler angles to set directions
E=[0.25*pi 0.25*pi -0.25*pi];
[R_true,~]=euler2DCM(E); %The true directions for X, Y and Z axis
X=(R_true*X')'; %Rotate polyhedron

X=X+MU_true(ones(size(X,1),1),:); %Centre points around mean

%Add noise
n_std=0.2;  %Standard deviation
Xn=X+n_std.*randn(size(X));

%%
% This is the true axis system
R_true

%%
% These are the true stretch factors
ellipStretchTrue

%%

[M,ellipStretchFit,R_fit,MU]=ellipsoidFit_centered(Xn,MU_true);

%%
% This is the fitted axis system. The system axes should be colinear with
% the true axes but can be oposite in direction. 
R_fit=R_fit(1:3,1:3)

%%
% These are the fitted stretch factors
ellipStretchFit

%%
% Building a fitted (clean) ellipsoid for visualization

%Create sphere
[F_fit,V_fit,~]=geoSphere(4,1);

%Transforming sphere to ellipsoid
V_fit_t=V_fit;
V_fit_t(:,end+1)=1;
V_fit_t=(M*V_fit_t')'; %Rotate polyhedron
V_fit=V_fit_t(:,1:end-1);

%%
% Visualizing results

MU=mean(X,1); %Origin for vectors
a=[7 7]; %Vector size

[Fq,Vq,Cq]=quiver3Dpatch(MU(1)*ones(1,3),MU(2)*ones(1,3),MU(3)*ones(1,3),R_fit(1,:),R_fit(2,:),R_fit(3,:),eye(3,3),a); %Fitted vectors
[Fq2,Vq2,Cq2]=quiver3Dpatch(MU(1)*ones(1,3),MU(2)*ones(1,3),MU(3)*ones(1,3),R_true(1,:),R_true(2,:),R_true(3,:),eye(3,3),a); %True vectors

figuremax(figColor,figColorDef);
title('The true (green) and fitted ellipsoid (red) and axis directions (solid, transparant respectively)','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(Xn,'k.','MarkerSize',15);

hp=patch('Faces',F,'Vertices',X);
set(hp,'FaceColor','g','FaceAlpha',1,'EdgeColor','k');

hp=patch('Faces',F_fit,'Vertices',V_fit);
set(hp,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none');

patch('Faces',Fq,'Vertices',Vq,'FaceColor','flat','FaceVertexCData',Cq,'FaceAlpha',1);
patch('Faces',Fq2,'Vertices',Vq2,'FaceColor','flat','FaceVertexCData',Cq2,'FaceAlpha',0.2,'EdgeColor','none');
camlight('headlight');
axis equal; view(3); axis vis3d; axis tight;  grid on;
set(gca,'FontSize',fontSize);
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
