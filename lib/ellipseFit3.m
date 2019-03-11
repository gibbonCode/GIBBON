function [e] = ellipseFit3(varargin)

% function [e] = ellipseFit3(V,optMethod,numSample)


%% Parse input

switch nargin
    case 1
        V=varargin{1};
        optMethod=2;
        numSample=[];
    case 2
        V=varargin{1};
        optMethod=varargin{2};
        numSample=[];
    case 3
        V=varargin{1};
        optMethod=varargin{2};
        numSample=varargin{3};
end

%%

mean_V=mean(V,1);
Vk=V-mean_V(ones(size(V,1),1),:);

[Q,~,~]=svd(Vk',0);
Vf=(Q'*Vk')';

%Fit ellipse
[A] = ellipseFit(Vf,optMethod,numSample);

[R,~]=euler2DCM([0 0 -A(5)]);

R1=eye(4,4);
R1(1:3,1:3)=R';

T1=eye(4,4);
T1(1:3,4)=[A(1) A(2) 0];

R2=eye(4,4);
R2(1:3,1:3)=Q;

T2=eye(4,4);
T2(1:3,4)=mean_V;

M=T2*R2*T1*R1';

%%

e.centre=M(1:3,4)'; 
e.radii=A(3:4);
e.axes=M(1:3,1:3);
e.tform=M;
 
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
