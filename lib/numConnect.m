function [N]=numConnect(F,V)

%Get patch face/vertex connectivity matrices
[~,IND_V]=tesIND(F,V,1);

%Count point connectivity
N=sum(IND_V>0,2); %Number of vertex neighbours