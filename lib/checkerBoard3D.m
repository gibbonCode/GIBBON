function M=checkerBoard3D(siz)

% function M=checkerBoard3D(siz)
% ------------------------------------------------------------------------
% This function creates a checkboard image of the size siz whereby elements
% are either black (0) or white (1). The first element is white.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% ------------------------------------------------------------------------

%%

%Coping with 1D or 2D input
if numel(siz)==2
    siz(3)=1; 
elseif numel(siz)==1
    siz(2:3)=[1 1];
end

%%
[I,J,K]=ndgrid(1:1:siz(1),1:1:siz(2),1:1:siz(3));
logic_ij=((iseven(I)| iseven(J)) & ((iseven(I)~=iseven(J))));
M=false(siz);
M(iseven(K))=logic_ij(iseven(K));
M(~iseven(K))=~logic_ij(~iseven(K));
 
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
