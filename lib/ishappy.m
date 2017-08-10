function L=ishappy(n)

N=n;
L=0;
NHistory=n;
while L==0
    N=sum(str2num(num2str(N).').^2); %Sum digits 
    if N==1 %n is a happy number
        L=1;
    elseif ismember(N,NHistory); %Kill we are in a loop
        break
    end
    NHistory=[NHistory N]; 
end
NHistory=[NHistory N]; 
 
%% <-- GIBBON footer text --> 
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
