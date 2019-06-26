function [varargout]=FFT_derivative(varargin)

% function [dY,fft_dY]=FFT_derivative(dt,Y,dimDer,derOrder,numPad)
% ------------------------------------------------------------------------
%
%
%
%
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        dt=varargin{1};
        Y=varargin{2};
        dimDer=1;
        derOrder=1;
        numPad=0;
    case 3
        dt=varargin{1};
        Y=varargin{2};
        dimDer=varargin{3};
        derOrder=1; 
        numPad=0;
    case 4
        dt=varargin{1};
        Y=varargin{2};
        dimDer=varargin{3};
        derOrder=varargin{4};
        numPad=0;
    case 5
        dt=varargin{1};
        Y=varargin{2};
        dimDer=varargin{3};
        derOrder=varargin{4};
        numPad=varargin{5};
end

%% Pad linearly

if numPad>0    
    siz=size(Y);
    [Y,indOriginal]=padLinDim(Y,numPad,dimDer,'both');
end

%%

% Get N
N=size(Y,dimDer); 

%Calculate frequencies
if iseven(N)    
    nx = ((-N/2:N/2-1)/N);
else    
    nx = ((-(N-1)/2:(N-1)/2)/N);
end
kx = ifftshift((2*pi./dt).*nx); 

%Reshape for bsxfun
if dimDer == 1
    kx = reshape(kx, N, 1);
else
    kx = reshape(kx, [ones(1,dimDer-1), N]);
end

fft_dY=(1i*kx).^derOrder.*fft(Y,[],dimDer);

dY=ifft(fft_dY,[],dimDer,'symmetric');

%% Crop padded back
if numPad>0  
    dY=reshape(dY(indOriginal),siz);
end

%% Collect output

varargout{1}=dY; 
varargout{2}=fft_dY; 
varargout{3}=kx; 
if numPad>0
    varargout{4}=indOriginal;
end

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
