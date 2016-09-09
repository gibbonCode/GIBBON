function [V,D]=cellEig(C)

% function [V,D]=cellEig(C)
% ------------------------------------------------------------------------
% Computes eigenvalues and eigenvectors for each matrix contained in the
% cell array C, i.e. [v,d]=eig(c) is executed for each cell entry. The
% output is two cell arrays, i.e. the cell V containing the eigenvectors
% and the cell D containing the eigenvalues.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 28/11/2013
% 2016/09/09 Updated documentation
%------------------------------------------------------------------------

[V,D]=cellfun(@eig,C,'UniformOutput',0);