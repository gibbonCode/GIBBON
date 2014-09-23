function [Vv]=shiftV(V,v)
Vv=V+v(ones(size(V,1),1),:);

