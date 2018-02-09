%% edgeListToCurve
% Below is a basic demonstration of the features of the |edgeListToCurve| function.

%% 

clear; close all; clc;

%% CREATING A REGION MESH

% Creating boundary curves 

%Boundary 1
ns=150;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
r=6+2.*sin(5*t);
[x,y] = pol2cart(t,r);
V1=[x(:) y(:)];

%Boundary 2
[x,y] = pol2cart(t,ones(size(t)));
V2=[x(:) y(:)+4];

%Boundary 3
[x,y] = pol2cart(t,2*ones(size(t)));
V3=[x(:) y(:)-0.5];

%Defining a region
regionCell={V1,V2,V3}; %A region between V1 and V2 (V2 forms a hole inside V1)

plotOn=1; %This turns on/off plotting

%Desired point spacing
pointSpacing=0.5; 

[F,V]=regionTriMesh2D(regionCell,pointSpacing,1,0);

%%
cFigure; hold on; 
gpatch(F,V,'g');

axisGeom; view(2);
drawnow;

%%
E=patchBoundary(F,V);
% E=E(all(ismember(E,find(V(:,1)<4.15)),2),:);


% [indList]=edgeListToCurve(E);
[G,G_iter]=tesgroup(E);

for q=1:1:size(G,2)
        
    E_now=E(G(:,q),:);

    plotV(V(E_now,:),'b.','markersize',25);


    [indListNow]=edgeListToCurve(E_now);
    
    plotV(V(indListNow,:),'b.-','markersize',25,'lineWidt',5);
    for w=1:1:numel(indListNow)
        plotV(V(indListNow(w),:),'r.','markersize',25);
        drawnow;
    end
    
end

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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
