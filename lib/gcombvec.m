function [B]=gcombvec(varargin)

numElements=cellfun(@numel,varargin);
numVecs=numel(varargin);
siz=prod(numElements);
B=zeros(numVecs,siz);

for q=1:1:numVecs
     vecNow=varargin{q};
     if q==1
         nRep1=1;
     else
         nRep1=prod(numElements(1:q-1));
     end     
     if q==numVecs
         nRep2=1;
     else
         nRep2=prod(numElements(q+1:end));         
     end
     w=repmat(vecNow,[nRep1 1]);    
     w=repmat(w(:),[nRep2 1]);          
     B(q,:)=w;
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
