function S=triplyPeriodicMinimal(varargin)

%Parse input
switch nargin
    case 2
        P=varargin{1};
        X=P(:,1); Y=P(:,2); Z=P(:,3); %Get coordinates
        typeStr=varargin{2};
    case 4
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        typeStr=varargin{4};
    otherwise
        error('False number of input arguments.');
end

%Evaluate metric on coordinates
switch typeStr
    case 'p' %Schwarz P
        S=(cos(X)+cos(Y)+cos(Z));
    case 'd' %Schwarz D
        S=(sin(X).*sin(Y).*sin(Z))+(sin(X).*cos(Y).*cos(Z))+(cos(X).*sin(Y).*cos(Z))+(cos(X).*cos(Y).*sin(Z));
    case 'g' %Schoen Gyroid
        S=(sin(X).*cos(Y))+(sin(Y).*cos(Z))+(cos(X).*sin(Z)); %Schoen Gyroid
    case 'n' %Neovius
        S=3*(cos(X)+ cos(Y)+ cos(Z))+ (4*cos(X).*cos(Y).*cos(Z));
    case 'w'
        S=cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X)-cos(X).*cos(Y).*cos(Z);
    case 'pw'
        S=(4.*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-3.*cos(X).*cos(Y).*cos(Z))+2.4;
    otherwise
        error('unknown surface type requested')
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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
