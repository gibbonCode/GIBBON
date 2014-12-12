function [A_norm]=vecnormalize(A)

if isvector(A)
    H=sqrt(sum(A.^2));
else
    H=sqrt(sum(A.^2,2));    
end
H(H<eps)=eps; %Avoids introducing NAN due to devision by zero for zero length vectors

A_norm=A./H(:,ones(size(A,2),1));

end