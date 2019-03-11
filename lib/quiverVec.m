function [varargout]=quiverVec(varargin)

% function [h]=quiverVec(P,V,vecSize,colorSpec,edgeColorOpt,quiverStyleOpt,alphaLevel)


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

if isempty(quiverStyleOpt)
    quiverStyleOpt=1;
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

h=gpatch(F,P,C,edgeColorOpt,alphaLevel);

if nargout>0
    varargout{1}=h;
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
