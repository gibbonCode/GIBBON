function [MU,k]=vonMisesStat(T)

MU=angle(mean(exp(1i.*T)));
Rsq=mean(cos(T)).^2+mean(sin(T)).^2;
R=sqrt(Rsq);

p=2; 
ki=R.*(p-Rsq)./(1-Rsq);
diffTol=ki/1000;
qIter=1;
while 1;
    kn=ki;
    A=besseli(p/2,ki)./ besseli((p/2)-1,ki);
    ki=ki-((A-R)./(1-A^2-((p-1)/ki)*A));                
    if abs(kn-ki)<=diffTol
        break
    end
    qIter=qIter+1
end
k=ki;

 
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
