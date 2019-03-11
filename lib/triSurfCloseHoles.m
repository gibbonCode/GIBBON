function [varargout]=triSurfCloseHoles(varargin)

% function [FT,VT,CT]=triSurfCloseHoles(F,V,E,G)
%-------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------

%% Parse input

switch nargin    
    case 2
        Fc=varargin{1};
        Vc=varargin{2};
        pointSpacing=[];
        Eb=[];
        G=[];        
    case 3
        Fc=varargin{1};
        Vc=varargin{2};
        pointSpacing=varargin{3};
        Eb=[];
        G=[];
    case 4
        Fc=varargin{1};
        Vc=varargin{2};
        pointSpacing=varargin{3};
        Eb=varargin{4};
        G=[];
    case 5
        Fc=varargin{1};
        Vc=varargin{2};
        pointSpacing=varargin{3};
        Eb=varargin{4};
        G=varargin{5};
end

%% Create edge lists and groupings if missing

if isempty(Eb)
    [Eb]=patchBoundary(Fc,Vc);
end

if isempty(G)
    [G]=tesgroup(Eb);
end

[~,~,Nv]=patchNormal(Fc,Vc);

%%

resampleCurveOpt=0;
interpMethod='natural'; %or 'natural'

FT={Fc};
VT={Vc};
for q=1:1:size(G,2)
    logicNow=G(:,q);
    indBoundaryNow=edgeListToCurve(Eb(logicNow,:));
    indBoundaryNow=indBoundaryNow(1:end-1);
    
    Vt=Vc(indBoundaryNow,:);
    
%     plotV(Vt,'y-','LineWidth',5);
    
    if ~isempty(pointSpacing)
        regionCell={Vt};
        [Fh,Vh]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);
        
        [Ebh,~,~]=patchBoundary(Fh,Vh);
        indBoundary_h=unique(Ebh(:));
        [D,indMin]=minDist(Vh(indBoundary_h,:),Vt);
        Vh(indBoundary_h,:)=Vt(indMin,:); %overwrite boundary points with originals
    else
        Vt_mean=mean(Vt,1);
        Vtc=Vt-Vt_mean(ones(size(Vt,1),1),:);
        [R_fit]=pointSetPrincipalDir(Vtc); %Get orientation
        Vtcr=(R_fit'*Vtc')';
        
        ns= size(Vtcr,1);
        C=[(1:ns)' [2:ns 1]'];
        
        %Initial Delaunay triangulation
        DT = delaunayTriangulation(Vtcr(:,1),Vtcr(:,2),C);
        Vh=Vt;%DT.Points;
        Fh=DT.ConnectivityList;
        
        %Remove faces not inside region
        L = isInterior(DT);
        Fh=Fh(L,:);
        
    end
            
    %Check normal directions
    N_now=mean(Nv(indBoundaryNow,:),1);
    [N_now_fit,~,~]=patchNormal(Fh,Vh);
    N_now_fit=mean(N_now_fit,1);
    if dot(N_now_fit,N_now)<0
        Fh=fliplr(Fh);
    end
     
%     gpatch(Fh,Vh,'g','k',1);
%     patchNormPlot(Fh,Vh);
%     drawnow;
    
    FT{end+1}=Fh;
    VT{end+1}=Vh;
end

%%
%Join element sets
[FT,VT,CT]=joinElementSets(FT,VT); 
 
%Merge nodes
[FT,VT]=mergeVertices(FT,VT,5);

%%

varargout{1}=FT;
varargout{2}=VT;
if nargout==3
    varargout{3}=CT;
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
