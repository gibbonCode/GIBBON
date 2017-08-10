function GB=gaborFilter(Vr,S,f,tagShift,d,siz,p)

S=S(:)';
Sm=ones(size(Vr,1),1)*S;

%Gaussian component
G=exp( -0.5.*(sum((Vr./Sm).^2,2)) );
G=G./sum(G(:)); %Normalising

%Harmonic component
H=real(exp(1i.*(2*pi*f*(Vr(:,d)-tagShift))));

%Defining Gabor filter form
GB=G.*H;

%Normalising
% normFactor=abs(2.*sqrt(2).*exp(-2.*f.^2.*S(d).^2.*pi.^(3/2)).*prod(S).*cos(phaseOffset));
% GB=GB./normFactor;
% GB=GB./abs(sum(GB(:))); 

% Defining Gabor filter
GB=reshape(GB,siz); 

%Flip sign if required
GB=p.*GB;




 
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
