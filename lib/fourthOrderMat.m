function [CM]=fourthOrderMat(C)

ind_C_all=1:numel(C);
[I,J,K,L]=ind2sub(size(C),ind_C_all);

%Create the 9x9 matrix indices
p=3*(I-1)+K;
q=3*(J-1)+L;

%Treat posible symbolic class
switch class(C)
    case 'sym'
        CM=sym(zeros(9,9));
    otherwise
        CM=zeros(9,9);
end

%Set values
[ind_pq]=sub2ind(size(CM),p,q);
CM(ind_pq(:))=C(:);