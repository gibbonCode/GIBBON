function [Vf,indFix]=curvePathOrderFix(V)

%%

%Deriving distance matrix
try
    D=dist(V,V');
catch
    D=distND(V,V);
end

D(eye(size(D))==1)=nan;

%Start loop
indFix=ones(1,size(V,1));
% hw = waitbar(0,'curvePathOrderFix:');
numPoints=size(V,1);
for qIter=2:1:numPoints
%     waitbar(qIter/numPoints);    
    [~,indMin]=nanmin(D,[],2); 
    indFix(qIter)=indMin(indFix(qIter-1));    
    D(:,indFix)=NaN;                
end
Vf=V(indFix,:);
% close(hw);