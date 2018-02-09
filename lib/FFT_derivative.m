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
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
