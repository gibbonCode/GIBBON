function Fmin=fmin_symbolic_form(Pn,OPT_struct,F)

% %Getting function variables 
% for i=1:1:numel(P)    
%     eval(['P',num2str(i),'=P(i);']);
% end

X=OPT_struct.FX;

%Evaluate function
P=sym('P',size(Pn)); %Create P as symbolic
Ff=subs(subs(F,P,Pn)); %Substitution of P and other

%Derive optimisation function value
switch OPT_struct.method
    case 'fminsearch'
        Fmin=sum((Ff(:)-OPT_struct.FY(:)).^2);
    case 'lsqnonlin' 
        Fmin=(Ff-OPT_struct.FY);
end
Fmin=double(Fmin);

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
