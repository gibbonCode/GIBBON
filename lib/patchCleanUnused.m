function [varargout]=patchCleanUnused(F,V)

logicValid =F>0;% %Treat 0,NaN,inf

numPoints=size(V,1);
indUni=unique(F(logicValid)); %Unique indices of used vertices
Vc=V(indUni,:); %Select relevant points

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(numPoints,1);
indFix2(indUni)=indFix1;
Fc=F; 
Fc(logicValid)=indFix2(F(logicValid));

%Output
varargout{1}=Fc;
varargout{2}=Vc;
varargout{3}=indFix2;