function p=csapsPar(varargin)

switch nargin
    case 1
        V=varargin{1};
        pw=1/1.1;
    case 2
        V=varargin{1};
        pw=varargin{2};
end

%%

% The interesting range for p is close to 1./(1+((h.^3)/6)). The following
% form is used introducing the factor f: p=1./(1+(((h.^3)/6)*f)). By using
% f=10 we obtain p=1./(1+((h.^3)/60)) which should result in a close
% following of the data. If instead f=0.1 is used, leading to
% p=1./(1+((h.^3)/0.6)), a smoother result is obtained.

%%
if pw<0
    pw=0;
end

if pw>1
    pw=1;
end

%%
f=(1/pw)-1;

%Calculate point spacings
hVec=sqrt(sum(diff(V,1,1).^2,2));
h=mean(hVec(:)); %Average point spacing

%Estimate smoothening parameter based on f and point spacing
p=1./(1+(((h.^3)/6)*f));

 
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
