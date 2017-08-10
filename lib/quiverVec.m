function [varargout]=quiverVec(varargin)


%% Parse input

switch nargin
    case 2
        P=varargin{1};
        V=varargin{2};
        vecSize=[];
        colorSpec=[];
        edgeColorOpt='none';
        quiverStyleOpt=1;
        alphaLevel=1;
    case 3
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=[];
        edgeColorOpt='none';
        quiverStyleOpt=1;
        alphaLevel=1;
    case 4
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};        
        edgeColorOpt='none';
        quiverStyleOpt=1;
        alphaLevel=1;
    case 5
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        edgeColorOpt=varargin{5};
        quiverStyleOpt=1;
        alphaLevel=1;
    case 6
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        edgeColorOpt=varargin{5};
        quiverStyleOpt=varargin{6};
        alphaLevel=1;
    case 7
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        edgeColorOpt=varargin{5};
        quiverStyleOpt=varargin{6};
        alphaLevel=varargin{7};
end
   
if isempty(edgeColorOpt)
    edgeColorOpt='none';
end

switch quiverStyleOpt
    case 1 %Depart from origin
        %Keep as is
    case 2 %Arrive at origin
        P=P-V;
    case 3 %Pass through origin
        P=P-(V/2);
    case 4 %Two-sided
        P=[P;P];
        V=[V;-V];
        if ~ischar(colorSpec) && size(colorSpec,1)>1
            colorSpec=[colorSpec;colorSpec];
        end
end

if size(P,2)==2
    P(:,3)=0;
end

if size(V,2)==2
    V(:,3)=0;
end

if numel(vecSize)==1
    vecSize=vecSize*ones(1,2);
end

if ischar(colorSpec)    
    [F,P,~]=quiver3Dpatch(P(:,1),P(:,2),P(:,3),V(:,1),V(:,2),V(:,3),[],vecSize);    
    C=colorSpec;
else
    if size(colorSpec,1)==1 %If only 1 color is provided
        colorSpec=colorSpec(ones(size(P,1),1),:); %copy for all vectors
    end
        
    [F,P,C]=quiver3Dpatch(P(:,1),P(:,2),P(:,3),V(:,1),V(:,2),V(:,3),colorSpec,vecSize);    
end

varargout{1}=gpatch(F,P,C,edgeColorOpt,alphaLevel);
 
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
