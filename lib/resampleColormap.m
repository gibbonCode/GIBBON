function [cmap_i]=resampleColormap(cmap,n)

ind=(1:1:size(cmap,1))';
ind_i=linspace(1,size(cmap,1),n)';
cmap_i=zeros(n,size(cmap,2));
for q=1:1:size(cmap,2);
    cmap_i(:,q)=interp1(ind,cmap(:,q),ind_i,'linear');
end
