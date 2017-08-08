function Is=unsharp(I,k,S,L,m)

if  ndims(I)==2
    [x,y] = meshgrid(linspace(-((k-1)/2),((k-1)/2),k));
    hg=exp(-(x.^2 + y.^2)/(2*S^2));
    hc=ones(k,k);

elseif  ndims(I)==3
    [x,y,z] = meshgrid(linspace(-((k-1)/2),((k-1)/2),k));
    hg=exp(-(x.^2 + y.^2 + z.^2)/(2*S^2));
    hc=ones(k,k,k);
end

hg(hg<eps*max(hg(:))) = 0;
sumh = sum(hg(:));
if sumh ~= 0;
    hg=hg/sumh;
end

Ig=convn(I,hg,'same'); %Gaussian blurred image

if m==1
    Is=I+(L.*(I-Ig));
elseif m==2
    hc=(1/numel(hc)).*hc;
    Im=convn(I,hc,'same');
    Id=(I-Im).^2; %A squared deviation image
    Iv=convn(Id,hc,'same'); %A local variance image
    Lv=L.*sqrt(1+Iv);
    Is=I+(Lv.*(I-Ig)); %Sharpened image
end

%% Alternative method
% %See also the matlab image processing toolbox command fspecial, h = fspecial('unsharp',0.2);
% a=0.2;
% h1=a/(a+1);
% h2=(1-a)/(a+1);
% h=[-h1 -h2 -h1;-h2 1+(4/(a+1)) -h2; -h1 -h2 -h1];
% Is=convn(I,h,'same');
% Is(Is<0)=0;
 
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
