function [CM]=fourthOrderMat(C)

ind_C_all=1:numel(C);
[I,J,K,L]=ind2sub(size(C),ind_C_all);

%Create the 9x9 matrix indices
p=3*(I-1)+K;
q=3*(J-1)+L;

%Treat posible symbolic class
switch class(C)
    case 'sym'
        CM=sym(zeros(9,9));
    otherwise
        CM=zeros(9,9);
end

%Set values
[ind_pq]=sub2ind(size(CM),p,q);
CM(ind_pq(:))=C(:);
 
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
