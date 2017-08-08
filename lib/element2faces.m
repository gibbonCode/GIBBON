function [F,C]=element2faces(E,C)

%Converts elements to faces enabling PATCH based visualisation colordata is
%copied for each face. Double faces for bordering elements may occur. 
%
%
%19/01/2012, Kevin Mattheus Moerman

switch size(E,2)
    case 4 %4 node tetrahedral elements
        F=[E(:,[1 2 3]);... 
           E(:,[1 2 4]);... 
           E(:,[2 3 4]);... 
           E(:,[3 1 4])]; 
       C=repmat(C(:),4,1);
    case 8 %8 node hexahedral elements
        F=[E(:,[1 2 3 4]);... %top
           E(:,[5 6 7 8]);... %bottom
           E(:,[1 2 6 5]);... %side 1
           E(:,[3 4 8 7]);... %side 2
           E(:,[2 3 7 6]);... %front
           E(:,[1 4 8 5]);]; %back       
       C=repmat(C(:),6,1);
    otherwise
        error('MATLAB:ELEMENT2FACES: size(E,1) not consistent with tetrahedral or hexahedral element');
end

end
 
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
