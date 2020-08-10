function [varargout]=lineVec(varargin)

% function [h]=lineVec(P,V,vecSize,colorSpec,quiverStyleOpt,lineWidth,alphaLevel)
% ------------------------------------------------------------------------
% 
% 
% 
% 
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        P=varargin{1};
        V=varargin{2};
        vecSize=[];
        colorSpec=[];        
        quiverStyleOpt=1;
        lineWidth=1;
        alphaLevel=1;
    case 3
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=[];        
        quiverStyleOpt=1;
        lineWidth=1;
        alphaLevel=1;
    case 4
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};        
        quiverStyleOpt=1;
        lineWidth=1;
        alphaLevel=1;
    case 5
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};        
        quiverStyleOpt=varargin{5};
        lineWidth=1;
        alphaLevel=1;
    case 6
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};        
        quiverStyleOpt=varargin{5};
        lineWidth=varargin{6};
        alphaLevel=1;
    case 7
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        quiverStyleOpt=varargin{5};
        lineWidth=varargin{6};
        alphaLevel=varargin{7};
end

if isempty(vecSize)
    vecSize=sqrt(sum(V.^2,2));
end

if isempty(colorSpec)
   colorSpec='k'; 
end

if ~isempty(P) || ~isempty(V)
        
    %Force data to be 3D
    if size(P,2)==2 
        P(:,3)=0;
    end    
    if size(V,2)==2
        V(:,3)=0;
    end
    
    %Expand vecsize to be 2x1 if needed
    if numel(vecSize)==1
        vecSize=vecSize*ones(1,2);
    end
    
    if ~isempty(vecSize)
        %Normalize vectors
        N=vecnormalize(V);        
        W=sqrt(sum(V.^2,2));
        W=W-min(W(:)); W=W./max(W(:)); 
        W=W*abs(vecSize(2)-vecSize(1))+vecSize(1);
        V=N.*W; %Scaled lengths
    end
        
    switch quiverStyleOpt
        case 1 %Depart from origin
            %Keep as is
        case 2 %Arrive at origin
            P=P-V;
        case 3 %Pass through origin
            P=P-V/2;
        case 4 %Two-sided
            P=[P;P];
            V=[V;-V];
            if ~ischar(colorSpec) && size(colorSpec,1)>1
                colorSpec=[colorSpec;colorSpec];
            end
    end
    
    %Prepare vertices and "edges" for plotting    
    
    if ~ischar(colorSpec) && size(colorSpec,1)>1
        colorSpec=[colorSpec;colorSpec]; %Double so vertex coloring is used (not face coloring)
    end
    
    VE=[P; P+V]; %Vertices for edges
    ind=(1:1:size(P,1))'; %Indices for points used to construct edges
    E=[ind ind+size(P,1)]; %Edges
    
    if ischar(colorSpec)
        h=patch('Faces',E,'Vertices',VE,'FaceColor','none','EdgeColor',colorSpec,'LineWidth',lineWidth,'EdgeAlpha',alphaLevel);
    else
        if size(colorSpec,1)==1 %If only 1 color is provided
            colorSpec=colorSpec(ones(size(VE,1),1),:); %copy for all points
        end        
        h=patch('Faces',E,'Vertices',VE,'FaceColor','flat','EdgeColor','flat','FaceVertexCData',colorSpec,'LineWidth',lineWidth,'EdgeAlpha',alphaLevel);
    end
    % h.EdgeLighting='flat';
    
else 
    h=[];    
end

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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
