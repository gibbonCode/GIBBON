function [varargout]=meshCleave(varargin)

% function [logicAt,logicAbove,logicBelow]=meshCleave(E,V,P,n,inclusiveSwitch)
% ------------------------------------------------------------------------
% This function creates logic arrays for the mesh components (e.g. elements
% or faces) which are at, above, or below a plane defined by the point P,
% and the normal direction n. 
% The optional inclusiveSwitch is a 2-component vector (default [0 0]) and
% sets how "inclusive", the below/above logic is, i.e. they set wether <
% and > is used ([0 0]), or <= and >= are used ([1 1]). A combination may
% also be used e.g. [1 0] results in below checks which features <= and
% above checks using >. 
% 
% Kevin Mattheus Moerman
% 2020/05/18: Created
%
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};        
        P=mean(V,1);
        n=[0 0 1]; 
        inclusiveSwitch=[0 0];
        toleranceLevel=max(eps(V(:)));
    case 3
        E=varargin{1};
        V=varargin{2};        
        P=varargin{3};
        n=[0 0 1];        
        inclusiveSwitch=[0 0];
        toleranceLevel=max(eps(V(:)));
    case 4
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};        
        inclusiveSwitch=[0 0];
        toleranceLevel=max(eps(V(:)));
    case 5
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};
        inclusiveSwitch=varargin{5};
        toleranceLevel=max(eps(V(:)));
    case 6
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};
        inclusiveSwitch=varargin{5};
        toleranceLevel=varargin{6};
end

%Check P
if isempty(P)
    P=mean(V,1); %Use default as mean of coordinates
end

%Check normal direction
if isempty(n)
    n=[0 0 1]; %Default z-direction
end
n=vecnormalize(n); %Normalize n vector

%% Construct rotation matrix
nz=[0 0 1]; %The z-vector
R=vecPair2Rot(n,nz);

%%

%Rotate coordinates
Vr=V-P(ones(size(V,1),1),:);
Vr=Vr*R';
Z=Vr(:,3);

%%

if inclusiveSwitch(1)==1
    logicVerticesBelow=Z<=0; %Logic for nodes below slice
else
    logicVerticesBelow=Z<-toleranceLevel; %Logic for nodes below slice
end

if inclusiveSwitch(2)==1
    logicVerticesAbove=Z>=0; %Logic for nodes above slice
else
    logicVerticesAbove=Z>toleranceLevel; %Logic for nodes above slice
end

logicAt=any(logicVerticesBelow(E),2) & any(logicVerticesAbove(E),2); %Logic for slice elements
logicBelow=all(logicVerticesBelow(E),2); %Logic for elements with nodes all below
logicAbove=all(logicVerticesAbove(E),2); %Logic for elements with nodes all above

%%

varargout{1}=logicAt;
varargout{2}=logicAbove;
varargout{3}=logicBelow;

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
