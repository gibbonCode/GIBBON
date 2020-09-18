function [f,v,c]=stochasticMicrostructure(inputStruct)

% -----------------------------------------------------------------------
% This function generates Stochastic Bicontinuous Microstructures
%
% Input structure and default values:
%
%   inputStruct.L=1; % characteristic length
%   inputStruct.Ns=80; % number of sampling points
%   inputStruct.Nw=120; % number of waves
%   inputStruct.q0=55; % wave number
%   inputStruct.relD=0.5; % relative density
%   inputStruct.anisotropyFactors=[1 1 1]; %Anisotropy factors
%   inputStruct.isocap=1; %Option to cap the isosurface
%
% Based on: Soyarslan et al. "3D stochastic bicontinuous
% microstructures: Generation, topology and elasticity"
% https://doi.org/10.1016/j.actamat.2018.01.005 
%
% Original author: Sebastien Callens, September 2020
%
% Change log: 
% 2020/09/17 Sebastien Callens: Created original code
% 2020/09/18 Kevin M. Moerman: Altered to be function
% 2020/09/18 Kevin M. Moerman: Vectorized code (removed need for for loop)
% 2020/09/18 Kevin M. Moerman: Added mesh union, node merging, mesh cleanup
% 2020/09/18 Kevin M. Moerman: Added face color data
% 2020/09/18 Kevin M. Moerman: Added anisotropy factors
% 2020/09/18 Kevin M. Moerman: Inverted face orientation
% -----------------------------------------------------------------------

%% Parse input

%Create default structure
defaultInputStruct.L=1; % characteristic length
defaultInputStruct.Ns=80; % number of sampling points
defaultInputStruct.Nw=120; % number of waves
defaultInputStruct.q0=55; % wave number
defaultInputStruct.relD=0.55; % relative density
defaultInputStruct.anisotropyFactors=[1 1 1]; 
defaultInputStruct.isocap=1; %Option to cap the isosurface

%Complete input with default if incomplete
[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

%Get parameters from input structure
L = inputStruct.L; % characteristic length
Ns = inputStruct.Ns; % number of sampling points
Nw = inputStruct.Nw; % number of waves
q0 = inputStruct.q0; % wave number
relD = inputStruct.relD; % relative density
anisotropyFactors = inputStruct.anisotropyFactors; 
isocap= inputStruct.isocap; 

%% Generation of points
[U,V,W] = meshgrid(linspace(0,L,Ns),linspace(0,L,Ns),linspace(0,L,Ns));
x = [reshape(U,[],1),reshape(V,[],1),reshape(W,[],1)];

%% Generation of random directions and phases
% Set anisotropy in x,y,z direction. If value=1: the entire range is
% sampled (no anisotropy)
xrange = anisotropyFactors(1); % If it is 1, we sample from the "whole possible range"
yrange = anisotropyFactors(2);
zrange = anisotropyFactors(3);

qi = [-xrange+2*xrange*randn(Nw,1),-yrange+2*yrange*randn(Nw,1),-zrange+2*zrange*randn(Nw,1)]; % create u,v,w components of qi in random direction (between -1 and 1)
qi = q0*vecnormalize(qi); % normalize rows and multiply with wave number
phii = 2*pi*randn(Nw,1);

%% Generation of function value at points

try
    funx=sqrt(2/Nw)*sum(cos(qi*x'+phii),1);
catch
    funx=sqrt(2/Nw)*sum(cos(qi*x'+phii(:,ones(size(x,1),1))),1);
end

% Old for-loop approach
% funx=zeros(1,size(x,1));
% for i = 1:size(x,1)
%     funx(i) = sqrt(2/Nw)*sum(cos(qi*x(i,:)'+phii));
% end

funx_grid = reshape(funx,size(U));

%% Controlling relative density
ksi = sqrt(2)*erfinv(2*relD-1);

%% Generation of level sets and exporting

%Iso-surface
[f,v] = isosurface(U,V,W,funx_grid,ksi);
c=zeros(size(f,1),1);

%Isocaps
if isocap==1
    %Compute isocaps
    [fc,vc] = isocaps(U,V,W,funx_grid,ksi);

    nc=patchNormal(fc,vc);
    cc=zeros(size(fc,1),1);
    cc(nc(:,1)<-0.5)=1;
    cc(nc(:,1)>0.5)=2;
    cc(nc(:,2)<-0.5)=3;
    cc(nc(:,2)>0.5)=4;
    cc(nc(:,3)<-0.5)=5;
    cc(nc(:,3)>0.5)=6;    
    
    %Join sets
    [f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc});
    
end

%% Merge nodes and clean-up mesh 

%Merge nodes
[f,v]=mergeVertices(f,v); 

%Check for unique faces
[~,indUni,~]=unique(sort(f,2),'rows');
f=f(indUni,:); %Keep unique faces
c=c(indUni);

%Remove collapsed faces
[f,logicKeep]=patchRemoveCollapsed(f); 
c=c(logicKeep);

%Remove unused points
[f,v]=patchCleanUnused(f,v); 

%Invert faces
f=fliplr(f); 

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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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


