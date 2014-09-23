function [M_pad]=padrep(M,repSize)

% function [M_pad]=padrep(M,repSize)
% ------------------------------------------------------------------------
% This function replicates the outer elements of a matrix 'M' outwards
% according to repSize to create the new matrix M_new
%
%
% 28/07/2008
% ------------------------------------------------------------------------


%%

numDims = numel(repSize);
idex_no   = cell(1,numDims);
for k = 1:numDims
    onesVector = ones(1,repSize(k));
    idex_no{k} = [onesVector 1:size(M,k) size(M,k)*onesVector];
end
M_pad = M(idex_no{:});

%% END