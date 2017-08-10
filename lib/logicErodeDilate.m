function [Lk]=logicErodeDilate(L,k,erodeDilateOption)

%Normalise structural kernel
k=k./sum(k(:));

%Convolve logic with kernel
M=double(L);
p=size(k)-1;
M_rep=zeros(size(M)+p);
M_rep(p(1)-1+(1:size(M,1)),p(2)-1+(1:size(M,2)),p(3)-1+(1:size(M,3)))=M;
Mk = convn(M_rep,k,'valid');
Mk=Mk./max(Mk(:)); %Scale max to 1

epsMax=max(eps(Mk(:)));

%Erode or dilate logic
switch erodeDilateOption
    case 'erode'
        Lk=(Mk>=(1-epsMax)); %Stayed nearly 1
    case 'dilate'    
        Lk=Mk>=(0+epsMax); %Became higher than 0
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
