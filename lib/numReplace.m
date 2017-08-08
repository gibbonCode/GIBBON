function A=numReplace(A,a,b)

% function A=numReplace(A,a,b)
% ------------------------------------------------------------------------
% Replacing the entries in a occuring in A with the entries in b. 
% 
% Example: 
% %Define input array
% A=[0,-6,3,0,0;...
%     -1,-4,-7,-9,-9;...
%     -4,-7,11,-5,12;...
%     10,-7,5,-7,13;...
%     0,11,-2,-6,2;...
%     12,4,2,-5,NaN];
% a=[2 -5 nan 0]; %Entries to replace
% b=[991 992 993 994]; %Entries to take their place
% B=numReplace(A,a,b); % Replacing the numbers using |numReplace|
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/01
%------------------------------------------------------------------------

if numel(a)~=numel(b)
    error('numel(a)~=numel(b)')
end

if numel(a)~=numel(unique(a))
    error('numel(a)~=numel(unique(a))')
end

%Dealing with NaN entries
L_nan=isnan(A); %Logic for the nan entries
nanRep=max(A(~L_nan))+1; %Value to replace nan's with (arbitrary number not in A)
A(L_nan)=nanRep; 
a(isnan(a))=nanRep;

%Sorting input vectors;
[a,I]=sort(a); %Values to replace
b=b(I); %Numbers to appear in output in place of the entries in a

%Replacing numbers in a with those in b
L_mem=ismember(A,a); %Finding members common to A and a
valRep=A(L_mem); %These are the (non-unique) entries to be replaces
[~,~,indReplace]=unique(valRep); %Use to unique operation to get indices into b
A(L_mem)=b(indReplace); %Now replace the entries
 
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
