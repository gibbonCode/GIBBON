%% parHipImplant
% Below is a demonstration of the features of the |parHipImplant| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,C,curveSet]=parHipImplant(hipParStruct);|

%% Description
% This function creates a hip implant object based on an input structure
% containing control parameters. 
% The hip implant does not current represent any particular implant
% geometry in use in clinical practise. Instead it is currently merely an
% excercise in the parameterisation of this type of objects and can be used
% to develop approximate models investigating him implant loading.

%% Examples

%%
% Plot settings

fontSize=20;
lineWidth=5;

%% 

%Define the parameter set
hipParStruct.ballRadius=20;
hipParStruct.stickRadius=7;
hipParStruct.stickLength=21;
hipParStruct.stickLengthStraight=hipParStruct.stickLength-6;
hipParStruct.neckRadius=15;
hipParStruct.neckEllipseScale=2;
hipParStruct.collarThickness=3; 
hipParStruct.loftOffset=20;
hipParStruct.loftLenght=40;
hipParStruct.stemRadius=8;
hipParStruct.stemLength=50;
hipParStruct.stemAngle=0.25*pi;
hipParStruct.pointSpacing=3;

[F,V,C,curveSet]=parHipImplant(hipParStruct);

%%
% Visualizing model

plotColors=gjet(numel(curveSet));

cFigure; 
gtitle('The parameterized hip implant');
subplot(1,2,1); hold on;
hp=gpatch(F,V,C,'k',1,0.5);
legend(hp,{'Surface'},'Location','NorthWest');
axisGeom(gca,fontSize); 
camlight headlight; view(2);
colormap gjet; icolorbar; 

subplot(1,2,2); hold on;
h=gobjects(numel(curveSet),1);
for q=1:1:numel(curveSet)
    h(q)=plotV(V(curveSet{q},:),'k-','LineWidth',lineWidth);
    h(q).Color=plotColors(q,:);
end
h(end+1)=gpatch(F,V,'w','none',0.25);
legend(h,{'Curve 1','Curve 2','Curve 3','Curve 4','Surface'},'Location','NorthWest');
axisGeom(gca,fontSize); view(2);
camlight headlight;

gdrawnow; 

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
