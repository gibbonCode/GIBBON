function [varargout]=quiverLine(varargin)

% function [h]=quiverLine(P,V,vecSize,colorSpec,lineWidth,quiverStyleOpt,alphaLevel)
% ------------------------------------------------------------------------
% Plots the vectors V with origins P as patch data arrows. The vectors size
% can be specified as vecSize allong with the colors in colorSpec. Other
% options include the edge color, quiver style, and transparency level.
%
%
% Kevin Mattheus Moerman
%
% Change log:
% 2020/07/14 Added single origin support
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        P=varargin{1};
        V=varargin{2};
        vecSize=[];
        colorSpec=[];
        lineWidth=2;
        quiverStyleOpt=1;
        alphaLevel=1;
    case 3
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=[];
        lineWidth=2;
        quiverStyleOpt=1;
        alphaLevel=1;
    case 4
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        lineWidth=2;
        quiverStyleOpt=1;
        alphaLevel=1;
    case 5
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        lineWidth=varargin{5};
        quiverStyleOpt=1;
        alphaLevel=1;
    case 6
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        lineWidth=varargin{5};
        quiverStyleOpt=varargin{6};
        alphaLevel=1;
    case 7
        P=varargin{1};
        V=varargin{2};
        vecSize=varargin{3};
        colorSpec=varargin{4};
        lineWidth=varargin{5};
        quiverStyleOpt=varargin{6};
        alphaLevel=varargin{7};
end

if ~isempty(P) && ~isempty(V)

    if size(P,1)==1 %Expand to match V if single origin is given
        P=P(ones(size(V,1),1),:);
    end
    
    R_vec=sqrt(sum(V.^2,2)); %Magnitude

    if ~isempty(vecSize)
        if numel(vecSize)==1
            vecSize=vecSize*ones(1,2);
        end

        %Line lengths
        if  min(R_vec(:))==max(R_vec(:)) %If all radii are equal, or if just 1 vector is used
            lineLengths=vecSize(2)*ones(size(R_vec)); %All lengths become 2nd entry
        else
            %Scale lengths between 1st and 2nd entry
            lineLengths=R_vec-min(R_vec(:));
            lineLengths=vecSize(1)+((lineLengths./max(lineLengths(:))).*(vecSize(2)-vecSize(1)));
        end
        V=vecnormalize(V).*lineLengths; %Set vector lenghts
    end

    if size(P,2)==2
        P(:,3)=0;
    end

    if size(V,2)==2
        V(:,3)=0;
    end
    
    if isempty(quiverStyleOpt)
        quiverStyleOpt=1;
    end

    %Create "edges"
    E=[(1:size(V,1))' (1:size(V,1))'+size(V,1)];

    switch quiverStyleOpt
        case 1 %Depart from origin
            PE=[P; P+V];
        case 2 %Arrive at origin
            PE=[P-V; P];
        case 3 %Pass through origin
            PE=[P-V/2; P+V/2];
        case 4 %Two-sided
            PE=[P-V; P+V];
    end

    h=gedge(E,PE,colorSpec,lineWidth,alphaLevel);
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
