function [IND_F,IND_V]=patchIND_old(F,V,formOpt)
%Check if formOpt is provided
if nargin==2 %DEFAULT
    formOpt=1;  %Output is cropped and full
end

JF=(1:1:size(F,1))'*ones(1,size(F,2)); %Face indices copied for each vertex entry
JF=JF(:); %Face indices as column
IF=F(:); %The vertex indices in the face matrix as column
%The above generates a column of face numbers indices (JF) which can be
%compared to a column of corresponding vertex indices (IF).

%Creating a sparse array where on the rows IF and columns JF we place the
%number JF. Rows indicate vertex index and the colum the face index
VF_IND_sp=sparse(IF,JF,JF,size(V,1),size(F,1));

%Creating a sparse array where on the rows IF and columns JF we place the
%number IF. Rows indicate vertex index and the colum the face index
FV_IND_sp=sparse(IF,JF,IF,size(V,1),size(F,1));

%Finding number of times vertices are used
V_count=full(sum(FV_IND_sp>0,2)); %Vertex use count
V_count_max=max(V_count(:)); %Max. vertex use count

%Preparing index matrix with a number of rows equal to the number of points
%and V_count_max columns. But transposed.
IND=(ones(size(V,1),1)*(1:1:V_count_max))'; %Allocate memory for index array using ones
L=(IND<=(V_count*ones(1,V_count_max))'); %Logic for 1<=V_count, i.e. V_count>=1 i.e. points used for faces


VF_IND_sp=VF_IND_sp'; %tranpose of sparse array
IND_F=zeros(V_count_max,size(V,1)); %Face index matrix allocation as zeros
[I,~,~] = find(VF_IND_sp); %Row index of non-zero VF_IND_sp elements
IND_F(L)=I; %Where vertices are used IND_F is set equal to the row index of nonzero points in VF_IND_sp
IND_F=IND_F'; %Transpose

IND_V=[]; %Allocated as empty
for q=1:1:size(F,2) %For each column in F
    A=IND_F;
    A(A~=0)=F(A(A~=0),q);
    IND_V(:,end+1:end+size(A,2))=A;
end
IND_V=sparse(IND_V);

%Creating sparse array
[I,~,v] = find(IND_V);
Iv=[I v];
[Iv_uni, ~, ~] = unique(Iv,'rows');
I=Iv_uni(:,1); v=Iv_uni(:,2);
IND_V=sparse(I,v,v,size(IND_V,1),size(IND_V,1));

switch formOpt
    case 1
        %Sorting and cropping sparse array
        IND_V=sort(IND_V,2);
        [~,J,~] = find(IND_V);
        IND_V=full(IND_V(:,min(J):end));
        IND_L=(1:1:size(IND_V,1))'*ones(1,size(IND_V,2));
        IND_V(IND_V==IND_L)=0;
        IND_V=sort(IND_V,2);
        IND_V=IND_V(:,2:end);
    case 2
        IND_F=VF_IND_sp'; %Output sparse array instead
end
