function [varargout]=interp_spherical(varargin)

% function [Ri,Ci]=interp_spherical(T,P,R,Ti,Pi,interpMethod,numberInterpSteps)
% ------------------------------------------------------------------------
% The function |interp_spherical| interpolates in a spherical coordinate
% sytem using standard interp2 type interpolation methods or those based on 
% Delaunay tesselations in the angular space such as natural neighbour
% interpolation method. Standard spherical interpolation of this type
% creates artifacts at the poles. Hence |interp_spherical| splits the
% interpolation up into a number of steps (set by numberInterpSteps). The
% function aims to interpolate at the "equator" such that polar artifacts
% can be minimized. For each interpolation step the interpolation problem is
% rotated such that the currect "equatorial band" is centered at the
% equator. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2017/09/13 Removed demo mode. Improved input parsing
%
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 5        
        T=varargin{1};
        P=varargin{2};
        R=varargin{3};
        Ti=varargin{4};
        Pi=varargin{5};
        interpMethod='natural';
        numberInterpSteps=10;
    case 6        
        T=varargin{1};
        P=varargin{2};
        R=varargin{3};
        Ti=varargin{4};
        Pi=varargin{5};
        interpMethod=varargin{6};
        numberInterpSteps=10;
    case 7
        T=varargin{1};
        P=varargin{2};
        R=varargin{3};
        Ti=varargin{4};
        Pi=varargin{5};
        interpMethod=varargin{6};
        numberInterpSteps=varargin{7};
    otherwise
        error('Wrong number of input arguments');        
end
    
%%
if numberInterpSteps==1
    [Ri]=interp_spherical_part(T,P,R,Ti,Pi,interpMethod);
    switch nargout
        case 1
            varargout{1}=Ri;
        case 2
            varargout{1}=Ri;
            varargout{2}=Ri;
    end
else
    %Getting source vertices
    V=zeros(size(R,1),3);
    [V(:,1),V(:,2),V(:,3)] = sph2cart(T,P,R);
    
    %Getting target vertices
    Vi=zeros(size(Pi,1),3);
    [Vi(:,1),Vi(:,2),Vi(:,3)] = sph2cart(Ti,Pi,ones(size(Pi)));
    
    %Determine the rotation settings
    rotationRange=linspace(0,pi,numberInterpSteps+1);
    rotationRange=rotationRange(1:end-1);
    phiThreshold=diff(rotationRange(1:2))/2;
    
    %Step wise spherical interpolation around equator
    indDone=[];
    Ri=ones(size(Vi,1),1);        
    Ci=ones(size(Vi,1),1);
    for q=1:1:numberInterpSteps
        %Rotate so that points of interest are at equator
        [DCM,~]=euler2DCM([rotationRange(q) 0 0]); %Direction cosines matrix
        Vr=V*DCM; %Rotated source vertices
        Vri=Vi*DCM; %Rotated target vertices
        
        %Convert to spherical coordinates
        [Tr,Pr,Rr]=cart2sph(Vr(:,1),Vr(:,2),Vr(:,3)); %Source
        [Ti,Pi,~]=cart2sph(Vri(:,1),Vri(:,2),Vri(:,3)); %Target
        
        %Get indices for the current target vertices
        indDoNow=find(Pi>=-phiThreshold & Pi<=phiThreshold);
        indDoNow=indDoNow(~ismember(indDoNow,indDone)); %remove ones done already
        
        %Interpolate region
        [Ri_step]=interp_spherical_part(Tr,Pr,Rr,Ti(indDoNow),Pi(indDoNow),interpMethod);        
        Ri(indDoNow)=Ri_step; %Setting new radii        
        Ci(indDoNow)=q; %Store iteration index for current points
        
        indDone=unique([indDone; indDoNow]); %Adding vertix indices to done list        
    end
    
    varargout{1}=Ri;
    varargout{2}=Ci;
    
end

end

%%

function [Ri]=interp_spherical_part(T,P,R,Ti,Pi,interpMethod)
%Tesselate data above, below, left and right to aid interpolation (and
%avoid some polar- and edge artifacts)
T=[T+2*pi; T+2*pi; T+2*pi; T; T; T; T-2*pi; T-2*pi; T-2*pi];
P=[P-pi; P; P+pi; P-pi; P; P+pi; P-pi; P; P+pi];
R=repmat(R,[9,1]);

%Removing double points
fRound=1e5; %Rounding factor for unique test
[~,indUni,~]=unique(round([T P R]*fRound)/fRound,'rows');
P=P(indUni);
T=T(indUni);
R=R(indUni);

if strcmp(interpMethod,'natural') || strcmp(interpMethod,'linear') || strcmp(interpMethod,'nearest') %TriScatterdInterp function
    F_delaunay=scatteredInterpolant([T P],R,interpMethod); %interpolator
    Ri=F_delaunay([Ti Pi]);
elseif strcmp(interpMethod,'cubic') %Griddata function
    Ri = griddata(T,P,R,Ti,Pi,interpMethod);
else
    error('Invalid interpolation method. The following methods are supported: linear, nearest, and cubic')
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
