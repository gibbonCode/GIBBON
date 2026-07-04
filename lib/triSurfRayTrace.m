function [varargout]=triSurfRayTrace(varargin)

% function [P,indIntersect,d,TUV]=triSurfRayTrace(p_origin,nr,F,V,optionStruct)
%-------------------------------------------------------------------------
% This function performs ray-tracings using the rays or lines defined by Nr
% with origin P_origin. The rays are traced to the surface defined by the
% faces F and the vertices V. An option structure can be added to control
% the tolerance level to use to consider intersections on ray/line and on
% triangle.
% 
% DEFAULTS:
% defaultOptionStruct.tolEps = 1e-5; %Tolerance level
% defaultOptionStruct.triSide = 0; %Triangle sides to consider
% defaultOptionStruct.rayType = 'ray'; %Use 'ray' type rather than 'line'
% defaultOptionStruct.exclusionType = 'inclusive'; %Include within tolerance
% defaultOptionStruct.paired=0; 
%
% The output consists of: 
% P            : An nx3 array for the intersection point coordinates
% indIntersect : An nx2 array with indices for rays and faces or the first/second columns respectively. 
% d            : An nx1 array with distances from the ray origin to the intersection point
% TUV          : An nx3 array with the t-parameter, and the barycentric coordinates u and v as columns
%
% Change log: 
% 2021/07/22: KMM Created as replacement for triangleRayIntersection
%-------------------------------------------------------------------------

%% Parse input

%Create default option structure
defaultOptionStruct.tolEps = 1e-5;
defaultOptionStruct.triSide = 0;
defaultOptionStruct.rayType = 'ray'; % or line
defaultOptionStruct.exclusionType = 'inclusive'; % or exclusice
defaultOptionStruct.paired=0; % 1 for paired, 0 for non-paired

switch nargin
    case 4
        La=varargin{1};
        Lab=varargin{2};
        F=varargin{3};
        V=varargin{4};
        optionStruct=defaultOptionStruct;
    case 5
        La=varargin{1};
        Lab=varargin{2};
        F=varargin{3};
        V=varargin{4};
        optionStruct=varargin{5};
end

%Complete input structure with defaults if missing
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

%%

if size(La,1)~=size(Lab,1)
    error('Number of ray origins does not match number of ray vectors');
end

if optionStruct.paired==1
    %Check number of rays against number of faces
    if size(La,1)~=size(F,1)
        error('For paired analysis the number of rays should match the number of faces.');
    end
end

%% Access options

tolEps=optionStruct.tolEps;
triSide=optionStruct.triSide;
rayType=optionStruct.rayType;
exclusionType=optionStruct.exclusionType;

switch exclusionType
    case 'normal'
        tolInc=0;
    case 'inclusive'
        tolInc=tolEps;
    case 'exclusive'
        tolInc=-tolEps;
end

%%

%Get face corner points
P0=V(F(:,1),:);
P1=V(F(:,2),:);
P2=V(F(:,3),:);

P01=P1-P0;
P02=P2-P0;

P01_P02_cross=cross(P01,P02,2);

if optionStruct.paired==1    
    [P,indIntersect,t,u,v]=evalRayIntersect(La,Lab,P0,P01,P02,P01_P02_cross,tolInc,tolEps,triSide,rayType);
    indRay=indIntersect;
    indFace=indIntersect;
    TUV=[t u v];
else
    P=[];
    indRay=[];
    indFace=[];
    TUV=[];    
    for q=1:1:size(La,1)
        LA=La(q*ones(size(F,1),1),:);
        LAB=Lab(q*ones(size(F,1),1),:);
        [P_now,indIntersect,t,u,v]=evalRayIntersect(LA,LAB,P0,P01,P02,P01_P02_cross,tolInc,tolEps,triSide,rayType);
        if ~isempty(P_now)
            P=[P;P_now];
            indRay=[indRay; q*ones(size(P_now,1),1);];
            indFace=[indFace; indIntersect];
            TUV=[TUV; [t u v]];
        end
    end
end

%% Collect output

varargout{1}=P; %The intersection points
varargout{2}=[indRay indFace]; %Ray and face indices

if nargout>2 %Also output distances to ray/line origins    
    if ~isempty(P)
        m=sqrt(sum(Lab(indRay,:).^2,2)); %Line lengths/magnitudes
        d=TUV(:,1).*m; %Compute distances        
    else
        d=[];
    end
    varargout{3}=d; %Add distances to output
end

if nargout==4 %Also output t,u,v
    varargout{4}=TUV;  %Add TUV to output
end

end

function [P,indIntersect,t,u,v]=evalRayIntersect(La,Lab,P0,P01,P02,P01_P02_cross,tolInc,tolEps,triSide,rayType)

    La_P0=La-P0;
    det_vec=dot(-Lab,P01_P02_cross,2); %det([-Lab' P01' P02'])
    
    switch triSide
        case 1
            logicDetSelect=det_vec>tolEps;
        case 0
            logicDetSelect=abs(det_vec)>tolEps;
        case -1
            logicDetSelect=det_vec<tolEps;
    end
    
    t=dot(P01_P02_cross(logicDetSelect,:),La_P0(logicDetSelect,:),2)./det_vec(logicDetSelect,:);    
    u=dot(cross(P02(logicDetSelect,:),-Lab(logicDetSelect,:),2),La_P0(logicDetSelect,:),2)./det_vec(logicDetSelect,:);    
    v=dot(cross(-Lab(logicDetSelect,:),P01(logicDetSelect,:),2),La_P0(logicDetSelect,:),2)./det_vec(logicDetSelect,:);
       
    %%
    %Create logic to check if intersections are on the triangle
    logicOnTriangle = u>=-tolInc & u<=1+tolInc & v>=-tolInc & v<=1+tolInc & u+v<=1+tolInc;
    
    switch rayType
        case 'ray'
            logicIntersect=logicOnTriangle;
        case 'line'
            logicOnline = t>=-tolInc & t<=1-tolInc;
            logicIntersect=logicOnline & logicOnTriangle;
    end
    P0_detSelect=P0(logicDetSelect,:);
    P01_detSelect=P01(logicDetSelect,:);
    P02_detSelect=P02(logicDetSelect,:);
    
    u=u(logicIntersect); 
    v=v(logicIntersect); 
    t=t(logicIntersect); 
    
    P=P0_detSelect(logicIntersect,:)+u.*P01_detSelect(logicIntersect,:)+v.*P02_detSelect(logicIntersect,:);
    
    indDetSelect=find(logicDetSelect);
    indIntersect=indDetSelect(logicIntersect);
    
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
