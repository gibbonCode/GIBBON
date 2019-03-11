function [F,V,faceBoundaryMarker]=quadBox(boxDim,boxEl)


%%
F=[]; V=[]; faceBoundaryMarker=[];

[Xtb,Ytb]=meshgrid(0:boxEl(1),0:boxEl(2));


[Ylr,Zlr]=meshgrid(0:boxEl(2),0:boxEl(3));
[f,v]=surf2patch(0*ones(size(Ylr)),Ylr,Zlr); %Left
f=fliplr(f);
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 1*ones(size(f,1),1)];

[f,v]=surf2patch(boxEl(1)*ones(size(Ylr)),Ylr,Zlr); %Right
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 2*ones(size(f,1),1)];

[Yfb,Zfb]=meshgrid(0:boxEl(1),0:boxEl(3));
[f,v]=surf2patch(Yfb,0*ones(size(Yfb)),Zfb); %Front
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 3*ones(size(f,1),1)];

[f,v]=surf2patch(Yfb,boxEl(2)*ones(size(Yfb)),Zfb); %Back
f=fliplr(f);
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 4*ones(size(f,1),1)];

[f,v]=surf2patch(Xtb,Ytb,0*ones(size(Xtb))); %Bottom
f=fliplr(f);
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 5*ones(size(f,1),1)];

[f,v]=surf2patch(Xtb,Ytb,boxEl(3)*ones(size(Xtb))); %Top
F=[F;f+size(V,1)]; V=[V;v;]; faceBoundaryMarker=[faceBoundaryMarker; 6*ones(size(f,1),1)];


%%
% Remove double nodes by merging
[F,V]=mergeVertices(F,V);

%%
% Scale coordinates
maxV=max(V,[],1);
V=V./maxV(ones(size(V,1),1),:);
V=V.*boxDim(ones(size(V,1),1),:);
V=V-boxDim(ones(size(V,1),1),:)/2;
 
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
