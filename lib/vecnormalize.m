function [A_norm]=vecnormalize(A)

if isvector(A)
    H=sqrt(sum(A.^2));
else
    H=sqrt(sum(A.^2,2));
    H(H==0)=eps;
    H=H*ones(1,size(A,2));
end

A_norm=A./H;

end