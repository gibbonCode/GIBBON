function [F,V]=patchcylinder(varargin)

%------------------------------------------------------------------------
% function [F,V]=patchcylinder(varargin)


% 2017/18/04
% 2017/18/04 Added varargin style with defaults for missing parameters
%------------------------------------------------------------------------

%%

switch nargin    
    case 1
        inputStruct=varargin{1};
    case 5
        inputStruct.cylRadius=varargin{1};
        inputStruct.numRadial=varargin{2};
        inputStruct.cylHeight=varargin{3};
        inputStruct.numHeight=varargin{4};
        inputStruct.meshType=varargin{5};        
    otherwise
        error('Wrong numer of input arguments');
end

if isfield(inputStruct,'cylRadius')
    cylRadius=inputStruct.cylRadius;
else
    cylRadius=1;
end
if isempty(cylRadius)
    cylRadius=1;
end

if isfield(inputStruct,'numRadial')
    numRadial=inputStruct.numRadial;
else
    numRadial=[];
end
if isempty(numRadial)
    numRadial=10;
end

if isfield(inputStruct,'cylHeight')
    cylHeight=inputStruct.cylHeight;
else
    cylHeight=2*cylRadius;
end
if isempty(cylHeight)
    cylHeight=2*cylRadius;
end

if isfield(inputStruct,'numHeight')
    numHeight=inputStruct.numHeight;
else
    numHeight=numRadial;
end
if isempty(numHeight)
    numHeight=numRadial;
end

if isfield(inputStruct,'meshType')
    meshType=inputStruct.meshType;
else
    meshType='quad';
end
if isempty(meshType)
    meshType='quad';
end

%%

t=linspace(0,2*pi,numRadial+1);
t=t(1:end-1);
x=cylRadius*cos(t);
y=cylRadius*sin(t);
Vc=[x(:) y(:)];
Vc(:,3)=0; 

cPar.numSteps=numHeight;
cPar.depth=cylHeight; 
cPar.patchType=meshType; 
cPar.dir=0;
cPar.closeLoopOpt=1; 

[F,V]=polyExtrude(Vc,cPar);

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
