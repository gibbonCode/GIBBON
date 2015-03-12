function [E]=gramSchmidtOrtho(E)

k=size(E,1);
for i=1:1:k
    E(i,:)=E(i,:)./norm(E(i,:));
    for j=i+1:1:k
        Ei=E(i,:);
        Ej=E(j,:);
        proj_ei_ej=dot(Ei,Ej).*Ei./norm(Ei);
        E(j,:)=E(j,:)-proj_ei_ej;
    end
end