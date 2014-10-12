function [B]=cellTrace(A)

% function [B]=cellTrace(A)
% ------------------------------------------------------------------------
% Computes the trace for each matrix contained in the cell array A, i.e.
% b=trace(a) is executed for each cell entry.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 28/11/2013
%------------------------------------------------------------------------

[B]=cellfun(@trace,A,'UniformOutput',0);