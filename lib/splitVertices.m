function [Fn,Vn]=splitVertices(F,V)

%%

%If vertices are shared they need to be disconnected by adding vertices

Vn=zeros(numel(F),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q); %The coordinate set
    XF=X(F); %The coordinates for each face vertex
    XFt=XF'; %transpose to prepare for column array
    Vn(:,q)=XFt(:); %Add column as part of new vertex matrix
end

%Fix indices in face matrix
indVn=1:numel(F);
Fn=reshape(indVn(:),size(F,2),size(F,1))';

