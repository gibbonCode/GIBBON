function [Et,Vt]=patchExtrude(F,V,zRange)

numSteps=numel(zRange);

%Deriving coordinates
Vt=repmat(V,numSteps,1);
Z_add=ones(size(V,1),1)*zRange; 
Z_add=Z_add(:);
Vt(:,3)=Vt(:,3)+Z_add;

%Replicated faces matrix
F_rep=repmat(F,numSteps-1,1);

%Fix indices since points are copied also
indFix=0:(numSteps-2); 
indFix=indFix(ones(1,size(F,1)),:);
indFix=indFix(:);
indFix=indFix(:,ones(1,size(F_rep,2)));

%Create element matrix
Et=[F_rep+(size(V,1)*indFix) F_rep+(size(V,1)*(indFix+1))];