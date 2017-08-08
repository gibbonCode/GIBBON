function GB=gabor_kernel(Vr,S,f,d,siz,p)

%Gaussian component
G=exp( -0.5.*((Vr(:,1)./S(1)).^2+(Vr(:,2)./S(2)).^2+(Vr(:,3)./S(3)).^2) );

%Harmonic component
H=p.*cos(2*pi*f*Vr(:,d));

% Defining Gabor filter
GB=reshape(G.*H,siz); 


% % Gaussian component
% G=exp( -0.5.*((Vr(:,1)./S(1)).^2+(Vr(:,2)./S(2)).^2+(Vr(:,3)./S(3)).^2) );
% G=G./sum(G(:));  %Normalise to sum to 1
% 
% %Harmonic component
% H=p.*cos(2*pi*f*Vr(:,d));
% 
% % Defining Gabor filter
% GB=reshape(G.*H,siz); 
% GB=GB./sum(abs(GB(:))); %Normalise such that sum(abs(GB))==1
% % GB=GB./abs(sum(GB(:))); %Normalise such that sum(GB)==1



 
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
