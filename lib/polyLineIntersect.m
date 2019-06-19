function [Vn]=polyLineIntersect(varargin)

% function [Vn]=polyLineIntersect(V,n_cut,cutLevel,isClosed)

switch nargin
    case 2
        V=varargin{1};
        n_cut=varargin{2};
        cutLevel=0;
        isClosed=0;
    case 3
        V=varargin{1};
        n_cut=varargin{2};
        cutLevel=varargin{3};
        isClosed=0;
    case 4
        V=varargin{1};
        n_cut=varargin{2};
        cutLevel=varargin{3};
        isClosed=varargin{4};
end

%%

if isClosed==1
    E=[(1:size(V,1))' [2:size(V,1) 1]'];
else
    E=[(1:size(V,1)-1)' (2:size(V,1))'];
end

VV=V(E(:,2),:)-V(E(:,1),:);
N=vecnormalize(VV);
P=V(E(:,1),:);

N_cut=n_cut(ones(size(P,1),1),:);

d=cutLevel-dot(P,N_cut,2);

s=dot(N,N_cut,2);
f=d./s;

a=(f.*N);
b=dot(N,a,2);
c=dot(N,VV,2);

f(b>c)=NaN;
f(f<0)=NaN;

Vn=P+(f.*N);

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
