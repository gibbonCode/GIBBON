%% patchVolumepatchVolume
% Below is a demonstration of the features of the |patchVolume| function

%%
clear; close all; clc;

%% Syntax
% |[volEst]=patchVolume(F,V);|

%% Description
% The patchVolume function computes the volume enclosed by the patch
% elements defined by the faces F and vertices V. 

%% Examples
%
makeAbsolute=0;
for testCase=1:9

switch testCase    
    case 1 %Trianglulated sphere
        r=2;
        ns=5;
        [F,V]=geoSphere(ns,r);         
        volTotalTrue=4/3*pi*r^3; %True theoretical volume
    case 2 %Trianglulated sphere
        r=2;
        ns=6;
        [F,V]=quadSphere(ns,r);         
        volTotalTrue=4/3*pi*r^3; %True theoretical volume

        %Shift randomly to create non-zero centred patch data
        V=V+10.*randn(1,3);
    case 3 %
        %Torus parameters
        r=1; %Sphere radius
        R=2; %Central radius
        nr=76;
        nc=150;
        [F,V]=patchTorus(r,nr,R,nc,'tri');
        volTotalTrue=2*pi^2*r*R; %True theoretical volume
    case 4 %
        %Torus parameters
        r=1; %Sphere radius
        R=2; %Central radius
        nr=76;
        nc=150;
        [F,V]=patchTorus(r,nr,R,nc,'quad');
        volTotalTrue=2*pi^2*r*R; %True theoretical volume
    case 5 %
        %Torus parameters
        r=1; %Sphere radius
        R=2; %Central radius        
        nr=76;
        nc=150;
        [F,V]=patchTorus(r,nr,R,nc,'honey');
        volTotalTrue=2*pi^2*r*R; %True theoretical volume
    case 6 % triangulated box
        boxDim=[5 2 1];
        pointSpacing=1;
        [F,V]=triBox(boxDim,pointSpacing); 
        volTotalTrue=prod(boxDim); %True theoretical volume
    case 7 % quadrangulated box
        boxDim=[5 2 1];
        [F,V]=quadBox(boxDim,[5 2 1]);
        volTotalTrue=prod(boxDim); %True theoretical volume
    case 8 % A mixed mesh consisting of pentahedra and hexahedra
        r=2;
        ns=5;
        [Ft,Vt]=geoSphere(ns,r);
        [V,F]=patch_dual(Vt,Ft);        
        volTotalTrue=4/3*pi*r^3; %True theoretical volume
    case 9 %Trianglulated sphere
        r=2;
        ns=5;
        [F,V]=geoSphere(ns,r);         
        F=fliplr(F);
        volTotalTrue=-4/3*pi*r^3; %True theoretical volume
end

%%

%Compute volume for the patch surface 
volEst=patchVolume(F,V,makeAbsolute)

volTotalTrue

%%
% Visualization

cFigure;
gpatch(F,V,'w');
axisGeom; camlight headlight;
gdrawnow;

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
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
