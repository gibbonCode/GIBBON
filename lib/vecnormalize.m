function [V_norm]=vecnormalize(V)

if isvector(V)
    H=sqrt(sum(V.^2));
else
    H=sqrt(sum(V.^2,2));    
end
logicInvalid=H==0;

V_norm=V./H(:,ones(size(V,2),1));
V_norm(logicInvalid,:)=0; %Set invalid lengths to 0
end