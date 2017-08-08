function [subIndOut,L_valid]=snapSubInd(subIndIn,siz)

%Initialize outputs
L_valid=false(size(subIndIn));
subIndOut=subIndIn;

%Fixing out of range indices
for q=1:1:size(subIndIn,2)
    subInd=subIndIn(:,q); %The current subscript index set
    L_valid(:,q)=(subInd>0) & (subInd<=siz(q)); %Logic for valid indices not out of range
    
    subIndToFix=subInd(~L_valid(:,q)); %Indices to fix
    
    %Snapping out of range indices to boundary
    subIndToFix=(subIndToFix.*(subIndToFix>1))+(subIndToFix<1); %Fix smaller than 1
    subIndToFix=(subIndToFix.*(subIndToFix<=siz(q)))+siz(q).*(subIndToFix>siz(q)); %Fix larger than siz(q)

    %Storing fixed indices in output
    subIndOut_current=subInd;
    subIndOut_current(~L_valid(:,q))=subIndToFix;
    subIndOut(:,q)=subIndOut_current;
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
