%% cap_patchCylinder
% Below is a demonstration of the features of the |cap_patchCylinder| function

%% Syntax
% |[Fm,Vm]=cap_patchcylinder(F,V,F,V2,nr,nz);|

%% Description 
% The |cap_patchCylinder| function closes (or caps) the patch data from a
% set of cylinders. 

%% Examples

%%
clear; close all; clc;

%%
% Plot settings
figColor='w'; 
figColorDef='white';
fontSize=10;
faceAlpha1=0.4;

%% Example: Creating input geometry
% Construct regular cylinder
r=10; 
nr=40; 
h=35;
nz=35;
ptype='tri';
[F,V]=patchcylinder(r,nr,h,nz,ptype);

%%
% Plottting model
hf1=figuremax(figColor,figColorDef);
title('Regular cylinder from patchcylinder','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
hpm=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','g','FaceAlpha',faceAlpha1);
axis equal; view(3); axis tight;  grid on; 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% 
% Deforming regular cylinder to create outer wall
[T,R,Z]=cart2pol(V(:,1),V(:,2),V(:,3));
R_add=2.*cos(2*pi.*V(:,3).*2.*(1/(h)));
R=R+R_add;
[V(:,1),V(:,2),V(:,3)]=pol2cart(T,R,Z);
X_add=3.*sin(2*pi.*V(:,3).*1.*(1/(h)));
Y_add=3.*cos(2*pi.*V(:,3).*1.*(1/(h)));
V(:,1)=V(:,1)+X_add;
V(:,2)=V(:,2)+Y_add;

%%
% Creating inner wall
[T,R,Z]=cart2pol(V(:,1),V(:,2),V(:,3));
a=0.5; b=2;
R_sub=b+(a+a.*cos(2*pi.*V(:,3).*3.*(1/(h))));
R2=R-R_sub;
V2=V;
[V2(:,1),V2(:,2),V2(:,3)]=pol2cart(T,R2,Z);

%% 
% Plotting walls
hf2=figuremax(figColor,figColorDef);
title('Two seperate walls','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
hpm1=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','r','FaceAlpha',0.5);
hpm2=patch('Faces',F,'Vertices',V2,'EdgeColor','k','FaceColor','b','FaceAlpha',0.5);
axis equal; view(3); axis tight;  grid on; 
set(gca,'FontSize',fontSize);
camlight headlight;
drawnow;

%% Example: Capping the cylinder geometry

[Fm,Vm]=cap_patchcylinder(F,V,F,V2,nr,nz);

hf3=figuremax(figColor,figColorDef);
title('Closed (capped) vessel model','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
hpm3=patch('Faces',Fm,'Vertices',Vm,'EdgeColor','k','FaceColor','r','FaceAlpha',1);
axis equal; view(3); axis tight;  grid on; set(gca,'FontSize',fontSize);
lighting flat; camlight('headlight');
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
