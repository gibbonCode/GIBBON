function [Zm]=SVD_filter(Z,P,T)

%Computer SVD
[U,S,V] = svd(Z);

%Normalising singular values
Ss=diag(S)-min(diag(S));
Ss=Ss./max(Ss);

%Creating smoothening parameters
p_max=P(1); %1=No blurring
p_min=P(2); %0=straight line fit
p=(Ss.*(p_max-p_min))+p_min; %Scale towards singular values

Zm=nan(size(Z));
for i=1:1:size(U,2);
    v=V(:,i); u=U(:,i); s=S(i,i); %components    
    if Ss(i)<T %Filter after threshold
                us = csaps(1:numel(u),u,p(i),1:numel(u))'; %Smooth u
                vs = csaps(1:numel(v),v,p(i),1:numel(v))'; %Smooth v
    else
        vs=v; us=u; %Keep unsmoothened
    end
    z=us*s*vs'; %sub-data
    Zm(:,:,i)=z;
end
Zm=sum(Zm,3);
 
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
