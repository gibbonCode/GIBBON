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
