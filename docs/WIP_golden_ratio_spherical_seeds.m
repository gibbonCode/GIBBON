clear; close all; clc; 

%% Plot settings

cMap=(kvw(250));

%%

N=100; %Number of seeds
phi=(1+sqrt(5))/2; %Ratio for seed depositing (e.g. golden ratio) 

nRefine=5; %Visualization sphere refinement steps from icosahedron(too dense=slow)

%% Visualization components

[F,V]=geoSphere(nRefine,1); %Sphere patch data for visualization
[~,~,Nv]=patchNormal(F,V); %Sphere vertex normal vectors for morphing

[vn,c]=sphereCoord(phi,N); %Get coordinates
[Vn,dn]=distMorph(V,vn,Nv); %Morph sphere 

%% Initialize figure

hf=cFigure; hold on; 
hp=gpatch(F,Vn,dn,'none',1);
hp.FaceColor='interp';
% hs1=scatterV(v,500,c,'filled');
hs1=plotV(vn,'k.','MarkerSize',50,'Visible',1);
ht=gtitle(['n=',num2str(N)]);
ht.FontSize=15; 
axisGeom; axis off; 
camlight headlight; lighting('gouraud');
colormap(cMap); % colorbar; 
axis manual;
drawnow;

%%
%Populate the animation structure

N_range=2:1:N; 
nSteps=numel(N_range);
animStruct.Time=N_range; %Create the time vector
for q=1:1:nSteps  
    [vn,cn]=sphereCoord(phi,N_range(q)); %Get coordinates    
    [Vn,dn]=distMorph(V,vn,Nv); %Morph sphere
    
    %Set entries in animation structure
    animStruct.Handles{q}=[hs1 hs1 hs1 hp hp ht]; %Handles of objects to animate
    animStruct.Props{q}={'XData','YData','ZData','Vertices','CData','String'}; %Properties of objects to animate
    animStruct.Set{q}={vn(:,1),vn(:,2),vn(:,3),Vn,dn,['n=',num2str(N_range(q))]}; %Property values for to set in order to animate
end

%% Animate 
anim8(hf,animStruct);

%%

function [v,c]=sphereCoord(f,N)

n=1:N-1;
z = 2*n./N - 1;
r = sqrt(1-z.^2);
theta = (2*pi)*n/f;
x = r.*sin(theta);
y = r.*cos(theta);

v=[x(:) y(:) z(:)];
c=n(:);

end

function [Vn,dn]=distMorph(V,vn,Nv)
    dn=minDist(V,vn);
    t=dn; 
    t=t-min(t(:)); 
    t=t./max(t(:));
    vs=Nv.*0.2.*t.^2;
    Vn=V-vs;
end
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
