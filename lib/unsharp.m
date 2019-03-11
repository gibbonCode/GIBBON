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
