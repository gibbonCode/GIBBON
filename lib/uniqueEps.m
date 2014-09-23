function [C,ia,ic]=uniqueEps(A,opt,epsOrderDiff)

% function [C,ia,ic]=uniqueEps(A,opt,epsOrderDiff)
% ------------------------------------------------------------------------
% Creates a unique set of entries or rows (if opt='rows') whereby entries
% are deemed the same if they vary by less than is implied by the
% epsOrderDiff. 
%
% See also: eps
%
% EXAMPLE
% [pi_uni,~,~]=uniqueEps([pi pi+eps(pi)],[],5)
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 08/04/2013
%------------------------------------------------------------------------

epsMax=max(eps(A),[],2);
epsFac=10.^(-(round(log10(epsMax))+epsOrderDiff));
if isempty(opt)
    [~,ia,ic] = unique(round(A.*(epsFac*ones(1,size(A,2)))));
    C=A(ia);
else %'rows' option assumed
    [~,ia,ic] = unique(round(A.*(epsFac*ones(1,size(A,2)))),opt);
    C=A(ia,:);
end




