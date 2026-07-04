%% hex8_hex20
% Below is a demonstration of the features of the |hex8_hex20| function

%% Syntax
% |[E_HEX20,V_HEX20,V_HEX20_cell,Fb_HEX20,Fb_HEX20_QUAD8]=hex8_hex20(E_HEX8,V_HEX8,V_HEX8_cell,Fb_HEX8);|

%% Description
% The |hex8_hex20| converts 4-node tetrahedral elements to 10-node
% tetrahedral elements. 

%% Examples

%%
clear; close all; clc;

% Plot settings
fontSize=15;
faceAlpha=0.5;
edgeColor='k';
edgeWidth1=2;
edgeWidth2=1;
markerSize1=75;
markerSize2=10;

%% CONVERSION FROM HEX8 TO HEX20, EXAMPLE FOR A SINGLE HEXAHEDRON
% Creating a single 4-node HEXAHEDRON
[V8,~]=platonic_solid(2,1); %q indicates solid type, r is the radius
HEX8=[1:8];
[F8,~]=element2patch(HEX8,[],'hex8');

%%
% Converting to a single 10-node HEXAHEDRON
[HEX20,V20,~]=hex8_hex20(HEX8,V8,{});
[F20,~]=element2patch(HEX20,[],'hex20');

%%
% Plotting elements
cFigure;
subplot(1,2,1); hold on;
title('A linear HEXAHEDRON','FontSize',fontSize);

hp=gpatch(F8,V8,'gw','k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize1;

patchNormPlot(F8,V8,0.75); %Plotting face normals

for q=1:1:size(HEX8,2)
    text(V8(HEX8(1,q),1),V8(HEX8(1,q),2),V8(HEX8(1,q),3),['  ',num2str(q)],'FontSize',fontSize*2,'color','b');
end

axisGeom(gca,fontSize);
axis off;
camlight('headlight'); lighting flat;

subplot(1,2,2); hold on;
title('A quadratic HEXAHEDRON','FontSize',fontSize);

hp=gpatch(F20,V20,'rw','k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize1;
patchNormPlot(F20,V20,0.75); %Plotting face normals

for q=1:1:size(HEX20,2)
    text(V20(HEX20(1,q),1),V20(HEX20(1,q),2),V20(HEX20(1,q),3),[' ',num2str(q)],'FontSize',fontSize*2,'color','b');
end

axisGeom(gca,fontSize);
axis off;
camlight('headlight'); lighting flat;

drawnow; 

%% CONVERSION FROM HEX8 TO HEX20, EXAMPLE FOR A HEXAHEDRON MESH WITH NODAL PARAMETERS
n=3;
for q=1:1:n
    [HEX8,V8]=subHex(HEX8,V8,1);
end

[F8,~]=element2patch(HEX8,[],'hex8');

%Create nodal result e.g. displacement and color
V4d=V8;
V4d(:,1)=V4d(:,1)+0.5.*sin(pi*V4d(:,3));
V4d(:,2)=V4d(:,2)+0.5.*cos(pi*V4d(:,3));
D4=V4d-V8;

C4=V8(:,1); %Color towards X

%%
% Converting to a single 10-node HEXAHEDRON
HEX8_cell={D4,C4};
[HEX20,V20,HEX20_cell]=hex8_hex20(HEX8,V8,HEX8_cell);
[F20,~]=element2patch(HEX20,[],'hex20');
D10=HEX20_cell{1};
C10=HEX20_cell{2};
V10d=V20+D10;

%%
% Plotting elements
hf=cFigure; % Open figure for plotting
subplot(2,2,1); hold on;
title('A linear HEXAHEDRON mesh','FontSize',fontSize);

hp=gpatch(F8,V8,'gw','k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize2;

axisGeom(gca,fontSize); axis off;
view([-50,12])
camlight('headlight'); lighting flat;

subplot(2,2,3); hold on;
title('Data on linear HEXAHEDRON mesh','FontSize',fontSize);

hp=gpatch(F8,V4d,C4,'k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize2;

axisGeom(gca,fontSize); axis off;
view([-50,12])
camlight('headlight'); lighting flat;

subplot(2,2,2); hold on;
title('A quadratic HEXAHEDRON mesh','FontSize',fontSize);

hp=gpatch(F20,V20,'rw','k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize2;

axisGeom(gca,fontSize); axis off;
view([-50,12])
camlight('headlight'); lighting flat;

subplot(2,2,4); hold on;
title('Mapped data on quadratic HEXAHEDRON mesh','FontSize',fontSize);

hp=gpatch(F20,V10d,C10,'k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize2;

axisGeom(gca,fontSize); axis off;
view([-50,12])
camlight('headlight'); lighting flat;

drawnow; 

%% CONVERSION FROM HEX8 TO HEX20, EXAMPLE FOR KEEPING TRACK OF BOUNDARY FACES

[F8,~]=element2patch(HEX8,[],'hex8');
[indBoundary]=tesBoundary(F8);
Fb4=F8(indBoundary,:);

[HEX20,V20,HEX20_cell,Fb10]=hex8_hex20(HEX8,V8,{},Fb4);
[F20,~]=element2patch(HEX20,[],'hex20');

% [indBoundary]=tesBoundary(F20,V20);
% Fb20=F20(indBoundary,:);

%%
% 

hf=cFigure; % Open figure for plotting
subplot(1,2,1); hold on;
title('A linear HEXAHEDRON mesh','FontSize',fontSize);

gpatch(F8,V8,'gw','k',faceAlpha);

gpatch(Fb4,V8,'rw','k',faceAlpha);
axisGeom(gca,fontSize); axis off;
camlight('headlight'); lighting flat;

subplot(1,2,2); hold on;
title('A quadratic HEXAHEDRON mesh','FontSize',fontSize);

hp=gpatch(Fb10,V20,'rw','k',faceAlpha);
hp.Marker='.';
hp.MarkerSize=markerSize2;

axisGeom(gca,fontSize); axis off;
camlight('headlight'); lighting flat;

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
