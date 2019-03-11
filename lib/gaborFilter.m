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
