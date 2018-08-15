function [L]=ind2logic(ind,siz)

%Assume vector if size has one entry
if numel(siz)==1
    siz=[siz 1];
end
L=false(siz);
L(ind)=1;