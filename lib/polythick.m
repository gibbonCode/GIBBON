function [X1,Y1,X2,Y2]=polythick(varargin)

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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
