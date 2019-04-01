function [varargout]=patchcylinder(varargin)

%------------------------------------------------------------------------
% function [F,V,C]=patchcylinder(optionStruct)

% 2017/04/18
% 2017/04/18 Added varargin style with defaults for missing parameters
% 2018/03/15 Updated input parsing based on default structure 
% 2018/03/15 Added option to output a closed cylinder
% 2018/03/15 Added varargout and color output
%------------------------------------------------------------------------

%% Parse input

%Create default structure
defaultOptionStruct.cylRadius=1;
defaultOptionStruct.numRadial=10;
defaultOptionStruct.cylHeight=2*defaultOptionStruct.cylRadius;
defaultOptionStruct.numHeight=[];
defaultOptionStruct.meshType='quad';
defaultOptionStruct.closeOpt=0;
switch nargin    
    case 1
        optionStruct=varargin{1};
    case 5
        optionStruct.cylRadius=varargin{1};
        optionStruct.numRadial=varargin{2};
        optionStruct.cylHeight=varargin{3};
        optionStruct.numHeight=varargin{4};
        optionStruct.meshType=varargin{5};        
    otherwise
        error('Wrong numer of input arguments');
end

%Complement/fix input with default
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

if isempty(optionStruct.numHeight)
    optionStruct.numHeight=round(optionStruct.cylHeight/((2*pi*optionStruct.cylRadius)/optionStruct.numRadial));
end

%Access parameters
cylRadius=optionStruct.cylRadius;
numRadial=optionStruct.numRadial;
cylHeight=optionStruct.cylHeight;
numHeight=optionStruct.numHeight;
meshType=optionStruct.meshType;
closeOpt=optionStruct.closeOpt;

%% Create cylinder

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

%% 

indTop=numHeight:numHeight:size(V,1);
indBottom=1:numHeight:size(V,1);

%% Cap ends if requested

if closeOpt==1
    [Ft,Vt]=regionTriMesh2D({V(indTop,[1 2])},[],0);
    Vt(:,3)=mean(V(indTop,3));
    
    [Fb,Vb]=regionTriMesh2D({V(indBottom,[1 2])},[],0);
    Vb(:,3)=mean(V(indBottom,3));
    Fb=fliplr(Fb);
    
    [F,V,C]=joinElementSets({F,Ft,Fb},{V,Vt,Vb});
    
    [F,V,~,ind2]=mergeVertices(F,V);
    indTop=ind2(indTop);
    indBottom=ind2(indBottom);
else
    C=ones(size(F,1),1);
end

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
varargout{4}=indTop;
varargout{5}=indBottom;

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
