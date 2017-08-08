function [W]=mask_design(Z)

% function [Im Jm Km Wm]=mask_design(I, J, K, Z)
% ------------------------------------------------------------------------
% This function designs the weights for a mask. It assumes the vector 'Z'
% contains the normalised intensities found in a region of interest. The
% mask weight vector 'W' is calculated such that the value sum(W.*Z) is
% minimum when the mask is centered on the region of interest and maximum
% when it is not.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 13/08/2008
% ------------------------------------------------------------------------


%%

W=Z;
W(Z==1)=-(10./(numel(Z(Z==1))*(Z(Z==1))));
W((Z>0 & Z<1))=(10./(numel(Z(Z>0 & Z<1))*(Z((Z>0 & Z<1)))));
W(Z==0)=999;

%% END
 
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
