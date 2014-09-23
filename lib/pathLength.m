function D=pathLength(V)

% function [F,V,C]=ind2patch(IND,M,ptype)
% ------------------------------------------------------------------------
%
% This function calculates the "current" curve patch length for each of the
% points of the curve defined by the input argument V. The curve may be
% multidimensional. The output D is a vector of size [size(V,1) 1] whereby
% each entry is defined as the sum of the point-to-point (Euclidean)
% distances (i.e. the curve patch length) leading up to that point. Hence
% it is clear that D(1) is 0 and D(end) is the total curve length. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 08/05/2013
%------------------------------------------------------------------------

%Compute distance metric
D=zeros(size(V,1),1); %Initialise with zeros (first values stays zero)
D(2:end)=cumsum(sqrt(sum(diff(V,1,1).^2,2)),1);


