function c=voigtUnMap(cVoigt)


if isvector(cVoigt) %assume that c is a 4th order tensor
    siz_c=[3 3];
    cVoigt(4:end)=(1/2).*cVoigt(4:end); %Undo doubling
else
    siz_c=[3 3 3 3];
end

[linearIndexVoigt,linearIndexFourthOrder]=tensor2voigtMap(zeros(siz_c));
    
switch class(cVoigt)
    case 'double'
        c=zeros(siz_c);
    case 'sym'
        c=sym(zeros(siz_c));
end

c(linearIndexFourthOrder)=cVoigt(linearIndexVoigt);

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
