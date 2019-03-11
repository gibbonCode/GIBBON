function [Q]=kabschRotationMatrix(V1,V2)

%2017/01/19 Fixed bug in relation to forcing right handed coordinate
%system. This will avoid inverting as well. 

%%
%Force input to 3D
if size(V1,2)==2
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

%% 

%Force input to 3D
if size(V1,2)==2
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

%Centering on the mean
V1_m=mean(V1,1);
V1=V1-V1_m(ones(size(V1,1),1),:);

V2_m=mean(V2,1);
V2=V2-V2_m(ones(size(V2,1),1),:);

%% Use Kabsch algorithm

A=V1'*V2;

[U,~,V] = svd(A);

d=sign(det(U'*V));
D=eye(3,3);
D(3,3)=d; 
Q=V*D*U';

 
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
