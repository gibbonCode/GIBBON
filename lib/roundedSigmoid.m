function [Vc]=roundedSigmoid(optionStruct)

%% Parse input 

defaultOptionStruct.w = 2;
defaultOptionStruct.h = 2;
defaultOptionStruct.p = 1; 
defaultOptionStruct.n = 100; 
defaultOptionStruct.np = 100; 
defaultOptionStruct.eps_level = 1e-9;

%Complete input structure with defaults if missing
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

% Access parameters
w = optionStruct.w; 
h = optionStruct.h; 
p = optionStruct.p; 
n = optionStruct.n; 
np = optionStruct.np; 
eps_level = optionStruct.eps_level;

%% Compute sigmoid

if p > eps_level 
    r = (h^2*p + h^2*p^2 + p*w^2 + sqrt(h^2*p^2*w^2 - h^2*p^4*w^2 + p^2*w^4))/(2*(h + 2*h*p + h*p^2));
    a = acos(1-((p*h)/(2*r)));

    t = linspace(1.5*pi,1.5*pi+a,np)';
    Vc1 = [r.*cos(t)-w/2 r.*sin(t)+r-h/2];

    t = linspace(pi/2+a,pi/2,np)';
    Vc2 = [r.*cos(t)+w/2 r.*sin(t)-r+h/2];

    if (1-p)<eps_level                
        Vc = [Vc1;Vc2(2:end,:)];
    else
        Vc = [Vc1;Vc2];
    end    
else 
    Vc = [-w/2 -h/2; w/2 h/2];
end

Vc = evenlySampleCurve(Vc,n);

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
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
