function [Vs]=subCurve(varargin)

% function [Vs]=subCurve(V,n,closeLoopOpt)
%------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/12/08 Updated to included closed loop option
%------------------------------------------------------------------------

%% Parse Input

switch nargin
    case 1
        V=varargin{1};
        n=1;
        closeLoopOpt=0;
    case 2
        V=varargin{1};
        n=varargin{2};
        closeLoopOpt=0;
    case 3
        V=varargin{1};
        n=varargin{2};
        closeLoopOpt=varargin{3};
    otherwise
        error('Wrong number of input arguments')
end
%%

if closeLoopOpt
    V(end+1,:)=V(1,:); %Add first point to end to close loop
end
    
if n==0 %No subdevision
    Vs=V;     
elseif n>1 %Subdevision of segments
    Vs=zeros(size(V,1)+(size(V,1)-1)*n,size(V,2));
    for q=1:1:size(V,2);
        X=V(:,q);
        XX=linspacen(X(1:end-1),X(2:end),n+2);
        XX=XX(:,1:end-1)';
        Vs(1:end-1,q)=XX(:);
        Vs(end,q)=V(end,q);
    end
elseif n<1 %Throw error
    error('n should be >=0')
end

if closeLoopOpt
    Vs=Vs(1:end-1,:); %Take away last point
end
 
%% <-- GIBBON footer text --> 
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
