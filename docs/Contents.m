function [dY,fft_dY]=FFT_derivative(Y,dt,dimDir,derOrder)

N=size(Y,dimDir); 

%Calculate frequencies
if iseven(N)    
    nx = ((-N/2:N/2-1)/N);
else    
    nx = ((-(N-1)/2:(N-1)/2)/N);
end
kx = ifftshift((2*pi/dt).*nx); 

%Reshape for bsxfun
if dimDir == 1
    kx = reshape(kx, N, 1);
else
    kx = reshape(kx, [ones(1,dimDir-1), N]);
end

fft_dY=bsxfun(@times,(1i*kx).^derOrder,fft(Y,[],dimDir));

dY=ifft(fft_dY,[],dimDir,'symmetric');
 
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
