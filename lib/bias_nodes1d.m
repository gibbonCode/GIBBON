function [x]=bias_nodes1d(x,f)

% function [x]=bias_nodes1d(x,f)
% ------------------------------------------------------------------------
% This function biases the spacing for the entries in x using the factor f
% such that the output array goes from x(1) to x(end) but biased using x^f
% (i.e. the spacing decreasing depending on the magnitude of x). 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

min_x=min(x,[],2)*ones(1,size(x,2)); 
x=x-min_x; 
max_x=max(x,[],2)*ones(1,size(x,2)); 
x=max_x.*((x.^f)./(max_x.^f)); 
x=x+min_x;
