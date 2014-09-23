function [MF]=fftnconv(M,F)

%K-space image
M_fft=fftn(M);

%K-space filter
F_fft=fftn(F);

%Filtering
MF_fft=M_fft.*F_fft;

%Reconstructing image
MF=ifftshift(ifftn(MF_fft));








