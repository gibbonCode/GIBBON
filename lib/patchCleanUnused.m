function [Fc,Vc,indFix2]=patchCleanUnused(F,V)

logicValid = (F>0);% & (F>size(V,1)); %Treat 0,NaN,inf, extreme entries

numPoints=size(V,1);
indUni=unique(F(logicValid)); %Indices of used vertices
Vc=V(indUni,:); %Select relevant points

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(numPoints,1);
indFix2(indUni)=indFix1;
Fc=F; 
Fc(logicValid)=indFix2(F(logicValid));


