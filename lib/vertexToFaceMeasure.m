function [CF]=vertexToFaceMeasure(F,CV)

CF=ones(size(F,1),size(CV,2));
for q=1:1:size(CV,2);
    cf=CV(:,q);            
    cv=cf(F);
    cv=mean(cv,2);
    CF(:,q)=cv; 
end
