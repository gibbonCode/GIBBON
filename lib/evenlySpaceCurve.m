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

