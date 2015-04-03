function [dY,fft_dY]=FFT_derivative(Y,t,dimDir,derOrder)

N=size(Y,dimDir); 

dt=diff(t,1,dimDir);
dt=dt(1);

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


