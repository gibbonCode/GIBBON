function [Fc,Vc,indFix2]=patchCleanUnused(F,V)

numPoints=size(V,1);
indUni=unique(F(:)); %Indices of used vertices
Vc=V(indUni,:); %Select relevant points

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(numPoints,1);
indFix2(indUni)=indFix1;
Fc=indFix2(F);


