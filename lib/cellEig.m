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
 
%% 
% ********** _license boilerplate_ **********
% 
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
