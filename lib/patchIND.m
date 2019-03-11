function [IND_F,IND_V]=patchIND(F,V,formOpt)

% function [IND_F,IND_V]=patchIND(F,V,formOpt)
% ------------------------------------------------------------------------
%
% This function computes the neighbouring faces (read elements in the case
% of >3D tesselations) and vertices for the input tesselation specified by
% F and V. The function generates matrices whereby row entries correspond
% to vertices and column entries correspond to neighbouring (according to
% the connectivity in the tesselation specified) vertices in the case IND_V
% and faces in the case of IND_F. 
% If formOpt==2 the output arrays IND_F and IND_V are sparse arrays of size
% [size(V,1),size(V,1)] and [size(V,1),size(F,1)] respectively. However if
% formOpt==1 (DEFAULT) or not specified then the output array is instead a
% full array whose size in the column direction is smaller and depends on
% the maximum number of vertex and face neighbours encountered. E.g. in a
% triangulated mesh where each vertex is connected to 6 vertices the array
% IND_V is thus [size(V,1),6] in size. If some locations have less than
% this maximum number of vertices zeros are used to fill up the array. The
% same line of arguments holds for the IND_F array where instead entries
% reflect neighbouring faces. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/01/26 %added option of outputting as sparse arrays
%------------------------------------------------------------------------

%%

%Check if formOpt is provided
if nargin==2 %DEFAULT
    formOpt=2;  %Output is cropped and full  
end

switch formOpt
    case 1
        sparseOpt=1;
    case 2
        sparseOpt=0;
end

[IND_F,IND_V]=tesIND(F,V,sparseOpt);

% %Check if formOpt is provided
% if nargin==2 %DEFAULT
%     formOpt=1;  %Output is cropped and full  
% end
%
% JF=(1:1:size(F,1))'*ones(1,size(F,2)); %Face indices copied for each vertex entry
% JF=JF(:); %Face indices as column
% IF=F(:); %The vertex indices in the face matrix as column
% %The above generates a column of face numbers indices (JF) which can be
% %compared to a column of corresponding vertex indices (IF). 
% 
% %Creating a sparse array where on the rows IF and columns JF we place the
% %number JF. Rows indicate vertex index and the colum the face index
% VF_IND_sp=sparse(IF,JF,JF,size(V,1),size(F,1));
% 
% %Creating a sparse array where on the rows IF and columns JF we place the
% %number IF. Rows indicate vertex index and the colum the face index
% FV_IND_sp=sparse(IF,JF,IF,size(V,1),size(F,1));
% 
% %Finding number of times vertices are used 
% V_count=full(sum(FV_IND_sp>0,2)); %Vertex use count
% V_count_max=max(V_count(:)); %Max. vertex use count
% 
% %Preparing index matrix with a number of rows equal to the number of points
% %and V_count_max columns. But transposed.
% IND=(ones(size(V,1),1)*(1:1:V_count_max))'; %Allocate memory for index array using ones
% L=(IND<=(V_count*ones(1,V_count_max))'); %Logic for 1<=V_count, i.e. V_count>=1 i.e. points used for faces
% 
% 
% VF_IND_sp=VF_IND_sp'; %tranpose of sparse array
% IND_F=zeros(V_count_max,size(V,1)); %Face index matrix allocation as zeros
% [I,~,~] = find(VF_IND_sp); %Row index of non-zero VF_IND_sp elements
% IND_F(L)=I; %Where vertices are used IND_F is set equal to the row index of nonzero points in VF_IND_sp
% IND_F=IND_F'; %Transpose
% 
% IND_V=[]; %Allocated as empty
% for q=1:1:size(F,2) %For each column in F
%     A=IND_F; 
%     A(A~=0)=F(A(A~=0),q);
%     IND_V(:,end+1:end+size(A,2))=A;
% end
% IND_V=sparse(IND_V);
% 
% %Creating sparse array
% [I,~,v] = find(IND_V);
% Iv=[I v];
% [Iv_uni, ~, ~] = unique(Iv,'rows');
% I=Iv_uni(:,1); v=Iv_uni(:,2);
% IND_V=sparse(I,v,v,size(IND_V,1),size(IND_V,1));
% 
% switch formOpt    
%     case 1
%         %Sorting and cropping sparse array
%         IND_V=sort(IND_V,2);
%         [~,J,~] = find(IND_V);
%         IND_V=full(IND_V(:,min(J):end));
%         IND_L=(1:1:size(IND_V,1))'*ones(1,size(IND_V,2));
%         IND_V(IND_V==IND_L)=0;
%         IND_V=sort(IND_V,2);
%         IND_V=IND_V(:,2:end);
%     case 2
%         IND_F=VF_IND_sp'; %Output sparse array instead
% end
 
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
