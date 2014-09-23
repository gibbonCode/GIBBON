function [X_mean]=wmean(X,W)

% function [X_mean]=wmean(X,W)
% ------------------------------------------------------------------------
% This function calculates the weighted mean of X using the weights defined
% in the vector W. 
%
% Uses the following formula: X_mean=sum(X.*W)./sum(W);
%
% 12/09/2008
% ------------------------------------------------------------------------

%%

X=X(:);
W=W(:);
X_mean=sum(X.*W)./sum(W);

%% END