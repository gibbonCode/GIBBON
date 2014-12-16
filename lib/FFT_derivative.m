function [dY,fft_dY]=FFT_derivative(Y,X,d,pCutOff)

siz=size(Y);

N=siz(d); 
sampleStepSize = (max(X(:))-min(X(:)))/(N-1); %Sampling step size
samplingFrequency = 1/sampleStepSize; %Sampling frequency

dimsY=1:ndims(Y);
dimsYd=dimsY(dimsY~=d);

sizLinspacen=[siz(dimsYd) 1];

k=linspacen(-samplingFrequency/2.*ones(sizLinspacen),samplingFrequency/2.*ones(sizLinspacen),N); %Frequencies space
k=permute(k,[dimsYd d]);

fft_Y=fftshift(fft(Y,[],d)); %Fourier transform of input
fft_dY=1i.*2*pi*k.*fft_Y; %Fourier transform of derivative
if ~isempty(pCutOff) %filter option
    fft_dY(abs(k)>(pCutOff.*(samplingFrequency/2)))=0; %Remove frequencies above kc
end
dY=ifft(ifftshift(fft_dY),[],d); %Obtain derivative through inverse Fourier transform
