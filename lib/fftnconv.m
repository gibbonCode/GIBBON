function [MF]=fftnconv(M,F)

%K-space image
M_fft=fftn(M);

%K-space filter
F_fft=fftn(F);

%Filtering
MF_fft=M_fft.*F_fft;

%Reconstructing image
MF=ifftshift(ifftn(MF_fft));








 
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
