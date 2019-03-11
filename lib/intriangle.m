function [L]=intriangle(TRI,N,V,P,d)

% clear all; close all; clc;
% TRI=[1 2 3];
% V=eye(3,3);
% V(3,3)=0;
% [N,Vn]=trinorm(TRI,V);
% %In, on edge, along edge
% % P=[0.25 0.25 0;];
% P=[0.9 0 0;];
% % P=[ 1.1 0 0;];
% d=5;
% [L]=intriangle(TRI,N,V,P,d)
% 
% figure; plot3(V(:,1),V(:,2),V(:,3),'b-');view(2); grid on; hold on; 
% plot3(P(1),P(2),P(3),'r.');
% A-----B
%  \    /
%   \  /
%    C

A=V(TRI(:,1),:); B=V(TRI(:,2),:); C=V(TRI(:,3),:);
% [N,Vn]=trinorm(TRI,V);

AB=vec_normalize(B-A);
AP=vec_normalize(P-A);

BC=vec_normalize(C-B);
BP=vec_normalize(P-B);

CA=vec_normalize(A-C);
CP=vec_normalize(P-C);

c1=vec_normalize(cross(AB,AP,2));
c2=vec_normalize(cross(BC,BP,2));
c3=vec_normalize(cross(CA,CP,2));

C=[sum(c1,2) sum(c2,2) sum(c3,2)]; 
C(isnan(C))=0;
C=pround(C,d); 
Lc=(any(C==0,2)& all(C>=0,2));

dot_1=pround(dot(c1,N,2),d); 
dot_2=pround(dot(c2,N,2),d); 
dot_3=pround(dot(c3,N,2),d); 
Ld=(dot_1==1) & (dot_2==1) & (dot_3==1);

L=Lc | Ld;
end
 
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
