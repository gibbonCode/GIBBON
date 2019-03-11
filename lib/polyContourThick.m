function [C]=polyContourThick(varargin)

switch nargin
    case 1
        V=varargin{1};
        k=[];
        v=[];
        pointSpacing=[];
        resampleMethod='linear';
    case 2
        V=varargin{1};
        k=varargin{2};
        v=[];
        pointSpacing=[];
        resampleMethod='linear';
    case 3
        V=varargin{1};
        k=varargin{2};
        v=varargin{3};
        pointSpacing=[];
        resampleMethod='linear';
    case 4
        V=varargin{1};
        k=varargin{2};
        v=varargin{3};
        pointSpacing=varargin{4};
        resampleMethod='linear';
    case 5        
        V=varargin{1};
        k=varargin{2};
        v=varargin{3};
        pointSpacing=varargin{4};
        resampleMethod=varargin{5};
end

%%

if size(V,2)==3
   warning('polyContourThick is for 2D polygons. 3D detected, ignoring 3rd dimension');  
   V=V(:,[1 2]);
end

D=sqrt(sum(diff(V,1,1).^2,2));
if isempty(k)
    k=mean(D)*2;
end

if isempty(v)
    v=mean(D)/3;
end

if isempty(pointSpacing)
    pointSpacing=mean(D);
end

%%
%Create coordinate matrices
minV=min(V,[],1)-2*k;
maxV=max(V,[],1)+2*k;
range_V=maxV-minV;
n=round(range_V./v);
[X,Y]=meshgrid(linspace(minV(1),maxV(1),n(1)),linspace(minV(2),maxV(2),n(2)));
Vg=[X(:) Y(:)];

%Derive distance based level-set image
D=minDist(Vg,V); %Nearest point distance between image grid points and polygon
M=reshape(D,size(X)); %Level-set image

%Compute contour
[C]=gcontour(X,Y,M,k,pointSpacing,resampleMethod);

 
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
