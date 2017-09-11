function [Ep,Vp,Cp]=pillowHex(E,V,C,shrinkFactor)

%Defive core element
[E_core,V_core,C_core]=scalePatch(E,V,C,shrinkFactor);

Vp=[V;V_core]; %Collect nodes 
E_core=E_core+size(V,1); %Fix node indices in E_core

%% FORMAT FOR FACES
% F =[E(:,[4 3 2 1]);... %top
%     E(:,[5 6 7 8]);... %bottom
%     E(:,[1 2 6 5]);... %side 1
%     E(:,[3 4 8 7]);... %side 2
%     E(:,[2 3 7 6]);... %front
%     E(:,[5 8 4 1]);]; %back
%C=repmat(C,6,1);
        
indTop=1:1:4;
indBottom=5:1:8;
indFront=[2 3 7 6];
indBack=[5 8 4 1];
indSide1=[1 2 6 5];
indSide2=[3 4 8 7];

%% THE TOP ELEMENT
%   top-of-core top-of-outer
E1=[E_core(:,indTop) E(:,indTop)];

%% THE BOTTOM ELEMENT
%   bottom-of-core bottom-of-outer
E2=[E_core(:,indBottom) E(:,indBottom)];

%% THE FRONT ELEMENT
%   bottom=front of E2 top=front of E1
E3=[E2(:,indFront) E1(:,indFront)];

%% THE BACK ELEMENT
%   bottom=back of E2 top=back of E1
E4=[E2(:,indBack) E1(:,indBack)];

%% THE SIDE1 ELEMENT
%   bottom=side1 of E2 top=side1 of E1
E5=[E2(:,indSide1) E1(:,indSide1)];

%% THE SIDE2 ELEMENT
%   bottom=side2 of E2 top=side2 of E1
E6=[E2(:,indSide2) E1(:,indSide2)];

%% FIX FACES ORIENTATION 

%% COLLECT ALL ELEMENTS
Ep=[E_core; E1; E2; E3; E4; E5; E6];
Cp=repmat(C,[7,1]);
% Ep=[E_core;E1;E2];

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
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
