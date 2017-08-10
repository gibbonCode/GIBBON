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
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
