function [E,Ev,Fe,Ve]=patch2edgeIND(F,V)

%% DERIVE EDGE MATRIX

% E=[]; %Can be improved through memory allocation and indexing in loop
% for i=1:1:size(F,2)-1;
%   E=[E; F(:,i) F(:,i+1)]; 
% end
% E=[E; F(:,end) F(:,1)]; 

%Format of column index in F
EColumnInd=[(1:size(F,2)); (1:size(F,2))];
EColumnInd=[EColumnInd(2:end) EColumnInd(1)];

%Derive edges matrix
E=F(:,EColumnInd)'; %Use index into F to create edges matrix
E=reshape(E,2,numel(E)/2)'; 

%%
E=sort(E,2); %Sort edge order
[E,~,ind2] = unique(E,'rows'); %Removing double edges, i.e. [1  4] = [4  1]

%%

Fe=reshape(1:numel(F),size(F,2),size(F,1))';
Fe=ind2(Fe);

ind_E=(1:size(E,1))'*ones(1,2);

Ve=sort((sparse(E(:),ind_E(:),ind_E(:),size(V,1),size(E,1))),2);
[~,J]=find(Ve);
Ve=full(Ve(:,min(J):end));

A=E(Ve(Ve>0),:);
B=zeros(size(Ve));
B(Ve>0)=A(:,1);
C=zeros(size(Ve));
C(Ve>0)=A(:,2);
D=(1:size(Ve,1))'*ones(1,size(Ve,2));
B(B==D)=0;
C(C==D)=0;
Ev=B+C;


 
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
