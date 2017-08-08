function p=csapsPar(varargin)

switch nargin
    case 1
        V=varargin{1};
        pw=1/1.1;
    case 2
        V=varargin{1};
        pw=varargin{2};
end

%%

% The interesting range for p is close to 1./(1+((h.^3)/6)). The following
% form is used introducing the factor f: p=1./(1+(((h.^3)/6)*f)). By using
% f=10 we obtain p=1./(1+((h.^3)/60)) which should result in a close
% following of the data. If instead f=0.1 is used, leading to
% p=1./(1+((h.^3)/0.6)), a smoother result is obtained.

%%
if pw<0
    pw=0;
end

if pw>1
    pw=1;
end

%%
f=(1/pw)-1;

%Calculate point spacings
hVec=sqrt(sum(diff(V,1,1).^2,2));
h=mean(hVec(:)); %Average point spacing

%Estimate smoothening parameter based on f and point spacing
p=1./(1+(((h.^3)/6)*f));

 
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
