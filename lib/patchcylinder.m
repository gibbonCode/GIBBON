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
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
