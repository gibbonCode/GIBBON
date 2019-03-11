function [varargout]=polythick(varargin)

% function [Xu Yu Xl Yl]=polythick(x,y,t,nSteps,resampleFactor)
% ------------------------------------------------------------------------
% This function thickens a polygon defined by x, y. It generates the upper
% (Xu, Yu) and lower (Xl, Yl) coordinates depending on the thickness t. 
%
% The slope is calculated at each coordinate and points are copied upwards
% and downwards (thickening) orthogonal to the local slope.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 27/08/2008 Created
% 2017/03/16 Updated for GIBBON
% 2017/03/16 Added varargout option

%------------------------------------------------------------------------

%% 

switch nargin
    case 2
        x=varargin{1};
        y=varargin{2};
        t=[];
        nSteps=[];
        resampleFactor=3;
    case 3
        x=varargin{1};
        y=varargin{2};
        t=varargin{3};
        nSteps=[];
        resampleFactor=3;
    case 4
        x=varargin{1};
        y=varargin{2};
        t=varargin{3};
        nSteps=varargin{4};
        resampleFactor=3;
    case 5
        x=varargin{1};
        y=varargin{2};
        t=varargin{3};
        nSteps=varargin{4};
        resampleFactor=varargin{5};
end

%%
V_curve=[x(:) y(:)];
V_curve(:,3)=0;

D=sqrt(sum(diff(V_curve,1,1).^2,2));
if isempty(t)
    t=mean(D);
end

if isempty(nSteps)
    nSteps=ceil(4*(t./mean(D)));
end

%%
V_curve1=V_curve;
V_curve2=V_curve;

for q=1:1:nSteps

    [~,~,N]=polyNormal(V_curve1);
    N(:,3)=0;
    V_curve1=V_curve1-(t/2/nSteps)*N;    
    
    [~,~,N]=polyNormal(V_curve2);
    N(:,3)=0;
    V_curve2=V_curve2+(t/2/nSteps)*N;    
    
    [V_curve1] = evenlySampleCurve(V_curve1,size(V_curve,1)*resampleFactor,'linear',0);
    [V_curve2] = evenlySampleCurve(V_curve2,size(V_curve,1)*resampleFactor,'linear',0);
end

[V_curve1] = evenlySampleCurve(V_curve1,size(V_curve,1),'linear',0);
[V_curve2] = evenlySampleCurve(V_curve2,size(V_curve,1),'linear',0);

X1=V_curve1(:,1);
Y1=V_curve1(:,2);

X2=V_curve2(:,1);
Y2=V_curve2(:,2);

%%

switch nargout
    case 2
        varargout{1}=V_curve1;
        varargout{2}=V_curve2;
    case 4
        varargout{1}=X1;
        varargout{2}=Y1;
        varargout{3}=X2;
        varargout{4}=Y2;
    otherwise
        error('Wrong number of output arguments. Use 2 or 4 output arguments');
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
