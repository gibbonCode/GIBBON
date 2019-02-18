%% truncatedOctahedronMesh
% Below is a demonstration of the features of the |truncatedOctahedronMesh| function

%%
clear; close all; clc;

%% Syntax
% |[F1c,F2c,F1,F2,C1,C2,VT]=truncatedOctahedronMesh(r,nCopies)|

%% Description
% Creates a truncated octahedron mesh where r sets the radias and nCopies
% (a 1x3 vector) sets the number of copies in the x, y, and z direction.
% The output consists of:
%
% F1c, F2c: the hexgonal and quadrilateral face cell arrays (1 cell entry
% per element).  
%
% F1, F2: the hexgonal and quadrilateral face arrays
%
% C1, C2: color/label data for the face arrays
%
% VT: the vertex array

%% 
% Plot settings
fontSize=15;
faceAlpha1=1;

%%

r=sqrt(5)/4; %Radii, results in a width of 1

n=5; %Desired number of copies in each direction 

%The actual input 
nCopies=[n n n+ceil((n+1)/2)+1]; %Number of offset copies


%%
[F1c,F2c,F1,F2,C1,C2,VT]=truncatedOctahedronMesh(r,nCopies);

%%
% Plotting results
cFigure; hold on;
suptitle('A mesh of truncated octahedra');
gpatch(F1,VT,C1,'k',faceAlpha1);
gpatch(F2,VT,C2,'k',faceAlpha1);
colormap(gjet);
axisGeom(gca,fontSize);
camlight('headlight'); 
drawnow; 
