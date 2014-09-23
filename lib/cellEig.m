function [V,D]=cellEig(C)

% function [V,D]=cellEig(C)
% ------------------------------------------------------------------------
% Computes eigenvalues and eigenvectors for each matrix contained in the
% cell array C, i.e. [v,d]=eig(c) is executed for each cell entry. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/11/2013
%------------------------------------------------------------------------

[V,D]=cellfun(@eig,C,'UniformOutput',0);