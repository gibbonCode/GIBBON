function [F,V]=patchNanFix(F,V)

numPoints=size(V,1); %Get original number of entries in vertex list

logicValid=~any(isnan(V),2); %Logic for valid vertices
logicValid_F=all(logicValid(F),2); %Logic vor valid faces

F=F(logicValid_F,:); %Get valid faces
V=V(logicValid,:); %Get valid vertices

%Fix indices in faces matrix
indFix1=1:nnz(logicValid);
indFix2=zeros(numPoints,1);
indFix2(logicValid)=indFix1;
F=indFix2(F); %Fix indices