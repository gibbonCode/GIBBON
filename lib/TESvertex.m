function [TES_dist]=TESvertex(TES,x,y,z)

x_vert=x(TES);
y_vert=y(TES);
z_vert=z(TES);

c_from=1:1:size(TES,2);
c_upto=[c_from(end) c_from(1:end-1)];
TES_dist=zeros(size(TES,1),numel(c_from));
for c=c_from;
    dc=sqrt((x_vert(:,c_from(c))-x_vert(:,c_upto(c))).^2+(y_vert(:,c_from(c))-y_vert(:,c_upto(c))).^2+(z_vert(:,c_from(c))-z_vert(:,c_upto(c))).^2);
    TES_dist(:,c)=dc;
end





