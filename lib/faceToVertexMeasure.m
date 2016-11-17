function [CV]=faceToVertexMeasure(F,V,CF)

[IND_F]=tesIND(F,V,0);

L=IND_F>0;

CV=ones(size(V,1),size(CF,2));
for q=1:1:size(CF,2)
    cf=CF(:,q);        
    cv=nan(size(IND_F));     
    cv(L)=cf(IND_F(L));
    cv=nanmean(cv,2);
    CV(:,q)=cv; 
end
