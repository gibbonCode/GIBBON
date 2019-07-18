function [E]=curveToEdgeList(N)
% function [E]=curveToEdgeList(N)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%%

if numel(N)==1 %the size of the list is specified
    indList=(1:1:N)';
elseif ismatrix(N) %ordered vertices are provided
    indList=(1:1:size(N,1))';
else %The indices are provided
    indList=N(:);
end

%%

E=[indList(1:end-1) indList(2:end)];


