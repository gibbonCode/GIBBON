function [Vg]=sweepCurveSmooth(varargin)

switch nargin
    case 5
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numSteps=varargin{5};
        p=[];
        f=0.1;        
    case 6    
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numSteps=varargin{5};
        p=varargin{6};
        f=0.1;
    case 7
        V1=varargin{1};
        V2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numSteps=varargin{5};
        p=varargin{6};
        f=varargin{7};
end

%%

d=sqrt(sum((V1-V2).^2))*f; %Distance between points

V=[V1-d*n1; V1; V1+d*n1; V2-d*n2; V2; V2+d*n2];

%Compute distance metric used for parametric representation
D=pathLength(V);

%Redefine distance metric for evenly spaced points
Dg=linspace(D(2),D(end-1),numSteps)';

p=csapsPar(V,p); %Scale invariant version of smoothening parameter
W=[1 1e9 1 1 1e9 1]; %Weights vector
Vg=zeros(numSteps,size(V,2)); %Initialize curve
for q=1:1:size(V,2) %Loop across dimensions
    v=V(:,q);
    pp = csaps(D,v,p,[],W)'; %Smoothened ppform
    Vg(:,q)=ppval(pp,Dg); %Evaluate ppform
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
