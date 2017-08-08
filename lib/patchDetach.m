function [Fc,Vc]=patchDetach(F,V,shrinkFactor)

Vc=zeros(size(F,1)*size(F,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    FX_mean=mean(FX,2);
    FX=((FX-FX_mean)*shrinkFactor)+FX_mean;
    Vc(:,q)=FX(:);
end
    
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));
 
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
