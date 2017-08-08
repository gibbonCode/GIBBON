function [P]=tesSmooth(TES,V,IND_V,cPar)

% function [P]=tesSmooth(TES,V,IND_V,cPar)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/06/02
%------------------------------------------------------------------------

%% 

%Get/set method
if isfield(cPar,'Method')
    smoothMethod=cPar.Method;
else
    smoothMethod='LAP'; %DEFAULT
end

%Smooth
switch smoothMethod
    case 'LAP'
        [P]=tesSmooth_LAP(TES,V,IND_V,cPar);
    case 'HC'
        [P]=tesSmooth_HC(TES,V,IND_V,cPar);
    otherwise
        error('Invalid smooth method specified');
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
