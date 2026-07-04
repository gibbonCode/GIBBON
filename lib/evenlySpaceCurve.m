function [Vg]=evenlySpaceCurve(varargin)

% function [Vg]=evenlySpaceCurve(V,pointSpacing,interpMethod,closeLoopOpt,indMust)
% ------------------------------------------------------------------------
% This function evenly sample the input curve V using the point spacing
% pointSpacing, and uses interpMethod as interpolation methods. The open
% closeLoopOpt will consider the curve closed. The indices indMust define
% points which must be included in the sampling. 
%
% Kevin Mattheus Moerman
% 
% 2020/05/06: Created
% 2021/04/12: Fixed bug for resampling of non-closed curve with must-points
% (end point removed from must points as it is already a must point for
% non-closed curves). 
%
% ------------------------------------------------------------------------

%% Parse input

V=varargin{1};
switch nargin
    case 2
        pointSpacing=varargin{2};
        interpMethod='pchip';
        closeLoopOpt=0;
        indMust=[];
    case 3
        pointSpacing=varargin{2};
        interpMethod=varargin{3};
        closeLoopOpt=0;
        indMust=[];
    case 4
        pointSpacing=varargin{2};
        interpMethod=varargin{3};
        closeLoopOpt=varargin{4};
        indMust=[];
    case 5
        pointSpacing=varargin{2};
        interpMethod=varargin{3};
        closeLoopOpt=varargin{4};
        indMust=varargin{5};
end

if isempty(closeLoopOpt)
    closeLoopOpt=0;
end

if isempty(interpMethod)
    interpMethod='pchip';    
end

if isempty(pointSpacing)
    if closeLoopOpt==1
        d=max(pathLength([V;V(1,:)]));
    else
        d=max(pathLength(V));
    end
    pointSpacing=d/size(V,1);
end

%%

%Remove first, if present, since it is already a must point
if ~isempty(indMust)
    indMust=indMust(indMust~=1); %Remove first if member
end

%Remove last for a non-closed curve since it is already a must point
if ~isempty(indMust) && closeLoopOpt==0
    indMust=indMust(indMust~=size(V,1));
end

if isempty(indMust) %Normal resample
    [Vg]=evenlySampleCurve(V,pointSpacing,interpMethod,closeLoopOpt,1);
else %Use must points

    %Append last point to curve if the loop is closed
    if closeLoopOpt==1
        V=[V;V(1,:)];
    end
    
    indMust=sort(indMust); %Sort must points
    
    Vg=V(1,:); %Initialize as first point
    indStart=1; %Initialize start index for curve segment as first point
    for q=1:1:numel(indMust)+1
        %Get end index for curve segment
        if q~=numel(indMust)+1
            indEnd=indMust(q); %Use the next must point
        else
            indEnd=size(V,1); %Use the end of the curve
        end
        %Get current curve segment
        V_now=V(indStart:indEnd,:);
        
        %Resample current segment
        [Vg_now]=evenlySampleCurve(V_now,pointSpacing,interpMethod,0,1);
        
        %Append new point set
        Vg=[Vg; Vg_now(2:end,:)]; 
        
        indStart=indEnd; %Update start for next segment
    end
    
    if closeLoopOpt==1        
        Vg=Vg(1:end-1,:);
    end
end

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
