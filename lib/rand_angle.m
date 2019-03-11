function a=rand_angle(siz)

% function a=rand_angle(siz)
% ------------------------------------------------------------------------
% This function generates a matrix or array of uniformly distributed random
% angles (of the size siz) in the range [0 2*pi]. The function operates by
% first creating angles in the range 0-pi/2. The angles are next converted
% to unit vectors in the unit circle. The X and Y components of these
% vectors are next independantly and randomly negated. This negating or
% flipping operation causes the vectors to uniformly but randomly span the
% entire cirle. Nexts the vectors are converted to an angle in the range [0
% 2*pi]. 
%
% 2018/12/04 Created and added to GIBBON
% ------------------------------------------------------------------------

%%

%Initial angles in range [0 pi/2]
a_temp=rand(siz)*(pi/2); %Random uniformly distributed angles

%Use angles to generate X component of unit vector. The Y component is not
%needed as the angles themselves are negated later to achieve the same
%effect. 
nX=cos(a_temp); %X components of the random unit vectors

%Flip vectors by randomly negating components
logicFlip_X=randi([0 1],siz)==1; %Logic for flipping X
logicFlip_Y=randi([0 1],siz)==1; %Logic for flipping Y
nX(logicFlip_X)=-nX(logicFlip_X); %Negate X components

%Derive angles for all vectors (now in range -pi pi degrees)
a=acos(nX); %Derive angles in range [0 pi]
a(logicFlip_Y)=-a(logicFlip_Y); %Negate angles to achieve Y component flipping

%Offset so all angles are positive in range 0 2*pi degrees
logicNegative=a<0;
a(logicNegative)=(2*pi)+a(logicNegative);

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
