function d=divisors(x,arg)

% Determine divisors of input x (which needs to be natural number)
% By default the output provides a vector containing all positive divisors
% but this can be altered depending on arg. 


%Default if no option is given
if exist('arg')==0
    arg='pos';
end

x=abs(x);
X=0:1:x;

d=X(ismember(X,abs(x)./X));
switch arg
    case 'pos'
    case 'neg'
        d=-fliplr(d);
    case 'all'
        d=[-fliplr(d) d];
    otherwise
        error('False input!');
end

end
 
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
