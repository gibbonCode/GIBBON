function [FT,VT,CT]=patchClosedCylinder(inputStruct)



% 2017/04/18

%%
if isfield(inputStruct,'numDigitKeep')
    numDigitKeep=inputStruct.numDigitKeep;
else
    numDigitKeep=[];     
end

%%
cylRadius=inputStruct.cylRadius;
pointSpacing=inputStruct.pointSpacing;
numRadial=round((2*pi*cylRadius)/pointSpacing);
cylHeight=inputStruct.cylHeight;
numHeight=round(inputStruct.cylHeight/((pointSpacing/2)*sqrt(3)));
numHeight=numHeight+iseven(numHeight); %Force uneven for 'tri' method

%%

t=linspace(0,2*pi,numRadial+1);
t=t(1:end-1);
x=cylRadius*cos(t);
y=cylRadius*sin(t);
Vc=[x(:) y(:)];
Vc(:,3)=0; 

cPar.numSteps=numHeight;
cPar.depth=cylHeight; 
cPar.patchType='tri'; 
cPar.dir=0;
cPar.closeLoopOpt=1; 

[F,V]=polyExtrude(Vc,cPar);
C=1*ones(size(F,1),1); %Color for side faces

%%
% Indices for the top and bottom points can be obtained as follows
indTop=numHeight:numHeight:size(V,1);
indBottom=1:numHeight:size(V,1);

%%
% The top and bottom can be meshed using |regionTriMesh2D|

[Ft,Vt]=regionTriMesh2D({V(indTop,[1 2])},[],0);
Vt(:,3)=mean(V(indTop,3));
Ct=2*ones(size(Ft,1),1); %Color for top faces

[Fb,Vb]=regionTriMesh2D({V(indBottom,[1 2])},[],0);
Vb(:,3)=mean(V(indBottom,3));
Fb=fliplr(Fb);
Cb=3*ones(size(Fb,1),1); %Color for bottom faces 

%%

%Join sets 
[FT,VT,CT]=joinElementSets({F,Ft,Fb},{V,Vt,Vb},{C,Ct,Cb});

%Merges nodes
[FT,VT]=mergeVertices(FT,VT,numDigitKeep);

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
