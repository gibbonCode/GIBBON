function [eulerVal]=eulerChar(F,V)

E=patchEdges(F,1);
eulerVal=size(V,1)-size(E,1)+size(F,1);
