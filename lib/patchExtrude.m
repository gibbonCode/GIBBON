function [Et,Vt]=patchExtrude(F,V,zRange)

numSteps=numel(zRange);

%Deriving coordinates
Vt=repmat(V,numSteps,1);
Z_add=ones(size(V,1),1)*zRange; 
Z_add=Z_add(:);
Vt(:,3)=Vt(:,3)+Z_add;

%Replicated faces matrix
F_rep=repmat(F,numSteps-1,1);

%Fix indices since points are copied also
indFix=0:(numSteps-2); 
indFix=indFix(ones(1,size(F,1)),:);
indFix=indFix(:);
indFix=indFix(:,ones(1,size(F_rep,2)));

%Create element matrix
Et=[F_rep+(size(V,1)*indFix) F_rep+(size(V,1)*(indFix+1))];
 
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
