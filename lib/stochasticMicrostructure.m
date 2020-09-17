%%%% Stochastic Bicontinuous Microstructures %%%%
%%%% Based on: Soyarslan et al. "3D stochastic bicontinuous
%%%% microstructures: Generation, topology and elasticity" %%%%
%%%% Author: Sebastien Callens
%%%% September 2020

clear; close all; clc;
addpath('C:\Users\sjpcallens\Documents\Matlab\Gibbon\lib');
%% Input
L = 1; % characteristic length
Ns = 80; % number of sampling points
Nw = 120; % number of waves
q0 = 55; % wave number
relD = 0.55; % relative density

%% Generation of points
[U,V,W] = meshgrid(linspace(0,L,Ns),linspace(0,L,Ns),linspace(0,L,Ns));
x = [reshape(U,[],1),reshape(V,[],1),reshape(W,[],1)];

%% Generation of random directions and phases
% Set anisotropy in x,y,z direction. If value=1: the entire range is
% sampled (no anisotropy)
xrange = 1; % If it is 1, we sample from the "whole possible range"
yrange = 1;
zrange = 0.4;

qi = [-xrange+2*xrange*randn(Nw,1),-yrange+2*yrange*randn(Nw,1),-zrange+2*zrange*randn(Nw,1)]; % create u,v,w components of qi in random direction (between -1 and 1)
qi = q0*normr(qi); % normalize rows and multiply with wave number
phii = 2*pi*randn(Nw,1);


%% Generation of function value at points
for i = 1:size(x,1)
    funx(i) = sqrt(2/Nw)*sum(cos(qi*x(i,:)'+phii));
end
funx_grid = reshape(funx,size(U));

%% Controlling relative density
ksi = sqrt(2)*erfinv(2*relD-1);

%% Generation of level sets and exporting
fv = isosurface(U,V,W,funx_grid,ksi);
p1 = patch(fv);
p1.FaceColor = [65,180,177]/255;
p2 = patch(isocaps(U,V,W,funx_grid,ksi));
p2.FaceColor = [65,180,177]/255; 
p2.FaceColor = p1.FaceColor;
p1.EdgeColor = 'none';
p2.EdgeColor = p1.EdgeColor;
axisGeom
view([45 30])
axis equal
axis off
camlight
material shiny

F1 = p1.Faces;
F2 = p2.Faces+max(reshape(F1,1,[]));
V1 = p1.Vertices;
V2 = p2.Vertices;

Ftot = [F1;F2];
Vtot = [V1;V2];
stlStruct.solidNames = {'test'};
stlStruct.solidVertices = {Vtot};
stlStruct.solidFaces = {Ftot};
stlStruct.solidNormals = {[]};

% export_STL_txt('Stochastic4.stl',stlStruct);
