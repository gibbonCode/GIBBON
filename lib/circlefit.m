function [Vc,R]=circlefit(V)

% function [Vc,R]=circlefit(V)
% ------------------------------------------------------------------------
% This function fits a circle to 2 or more points. The input is an array
% containing the points (vertices) V in 2D or 3D. Although 3D data is
% supported the Z-coordinates should all be equal such that the fitting
% problem remains 2D.
% The output consists of the circle centre Vc and radius R. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2019/11/21 Added to GIBBON
% ------------------------------------------------------------------------

%% Parse input

%Check if empty or not numeric
if isempty(V)
    error('Empty input')
elseif ~isnumeric(V)
    error('Input is not numeric');
end

%Force real
if ~isreal(V)
    V=real(V);
    warning('Complex data encountered. Complex parts ignored');
end

%Force double
if ~isa(V,'double')
    V=double(V); %Attempt conversion to double
end

%Check size
sizV=size(V);
if sizV(1)==1
    error('Only 1 point provided')
end

%Force 3D
if sizV(2)==2
    V(:,3)=0;
end

%% Fit cirlce
Vm=mean(V,1); %Get mean of point set
VD=V-Vm(ones(sizV(1),1),:); % shift points on mean
D=sqrt(sum(VD.^2,2)); %Compute distance to centre

if sizV(1)==2 %If only two points provided
    Vc=Vm; %Cirlc centre set to mean
    R=mean(D); %Radius set to mean distance
else %Fit circle to more than 2 points    
    % Scale data
    scaleFactor=max(max(D,eps));
    VD=VD./scaleFactor; D=D./scaleFactor;
    
    % Solve sytem
    VF=[2*VD(:,1) 2*VD(:,2) ones(sizV(1),1)]\D.^2;
    
    % Unscale, collect output
    Q=[Vm(1) Vm(2) 0]+[VF(1) VF(2) sqrt(VF(3)+hypot(VF(1),VF(2)).^2)]*scaleFactor;
    
    Vc=[Q(1) Q(2) Vm(3)]; %Cirle centre
    R=Q(3); %Radius
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
