function [varargout]=splitCurveSetMesh(varargin)

% function [F,V,curveIndices,faceMarker]=splitCurveSetMesh(V_cell,ns,patchType,smoothPar,splitMethod,w)
% ------------------------------------------------------------------------
%
% 2012
% 2018/11/08 Added varargout
% 2018/11/08 
% ------------------------------------------------------------------------

%% Parse input

switch nargin   
    case 5
        V_cell=varargin{1};
        numSteps=varargin{2};
        patchType=varargin{3};
        smoothPar=varargin{4};
        splitMethod=varargin{5};
        w=0;
    case 6
        V_cell=varargin{1};
        numSteps=varargin{2};
        patchType=varargin{3};
        smoothPar=varargin{4};
        splitMethod=varargin{5};
        w=varargin{6};
end

%%

V1=V_cell{1};
V2=V_cell{2};
V3=V_cell{3};

n1=size(V1,1);

%Compute distance metric used for parametric representation
D=pathLength(V1);
pointSpacing=max(D)./n1;

%%
%Derive contour centres
meanV1=mean(V1,1);
meanV2=mean(V2,1);
meanV3=mean(V3,1);

%Splice mother contour
switch splitMethod
    case 'ortho' 
        
        %Splice mother contour
        P1=meanV2-meanV3;
        N1=vecnormalize(P1);
        N2=cross(N1,[0 0 1]);
        [tN2,~]=cart2pol(N2(1),N2(2));       
        
        V1c=V1-meanV1(ones(size(V1,1),1),:); %Centered curve
        [tV1c,~]=cart2pol(V1c(:,1),V1c(:,2));        
       
        [~,indMin]=min(abs(tN2-tV1c));
        
        %Find opposite side
        indMax=indMin+round(size(V1,1)/2);
        if indMax>size(V1,1)
            indMax=indMax-size(V1,1);
        end
        
    case 'nearMid'
        
        [D,indMinD_all]=minDist(V2,V3);
        [~,indMinD_sub]=min(D);
        indMin_V2=indMinD_sub;
        indMin_V3=indMinD_all(indMinD_sub);
        
        V_test=0.5*(V2(indMin_V2,:)+V3(indMin_V3,:));
        D=sqrt(sum((V1-V_test(ones(size(V1,1),1),:)).^2,2));
        [~,indMin]=min(D);
        
        %Find opposite side
        indMax=indMin+round(size(V1,1)/2);
        indMaxRange=(indMax-round(size(V1,1)/4)):(indMax+round(size(V1,1)/4));
        indMaxRange(indMaxRange>size(V1,1))=indMaxRange(indMaxRange>size(V1,1))-size(V1,1);
          
        [~,indMinMaxRange]=min(D(indMaxRange));
        indMax=indMaxRange(indMinMaxRange);

end

indPoints=[indMin,indMax];
indMin=min(indPoints);
indMax=max(indPoints);

%Determine length of crossing patch
distCrossing=sqrt(sum((V1(indMin,:)-V1(indMax,:)).^2));

ind1=min(indPoints):max(indPoints);
ind2=[max(indPoints):size(V1,1) 1:min(indPoints)];

Vc1=V1(ind1,:);
Vc2=V1(ind2,:);

nc3=ceil(distCrossing/pointSpacing)+1;

[d,indDistMin]=minDist(V2,V3);
[~,indIndDistMin]=min(d); 
indClosest1=indIndDistMin; 
indClosest2=indDistMin(indIndDistMin); 

V_midSaddle=w*((V2(indClosest1,:)+V3(indClosest2,:))/2)+(1-w)*((V1(indMin,:)+V1(indMax,:))/2);

Vc3=[V1(indMin,:); V_midSaddle; V1(indMax,:);];
[Vc3]=evenlySampleCurve(Vc3,nc3+2,'spline',0);

V1s1=[flipud(Vc3); Vc1(2:end-1,:);];
V1s2=[(Vc3); Vc2(2:end-1,:); ];

%Find curves for split sides. Normally pairing V2 with V1s1 however
%switching may be required
meanV1s1=mean(V1s1,1);
d12=sqrt(sum((meanV1s1-meanV2).^2));
d13=sqrt(sum((meanV1s1-meanV3).^2));
if d13<d12 %Swith curve sets
    switchDone=1;
    V_temp=V1s2;
    V1s2=V1s1;
    V1s1=V_temp;
else
    switchDone=0;
end


%Derive "grid"
[V2s]=evenlySampleCurve(V2,size(V1s1,1),'pchip',1); %Resampling curve
[V2f]=minPolyTwist(V1s1,V2s); %Fix curve order
X=linspacen(V1s1(:,1),V2f(:,1),numSteps)';
Y=linspacen(V1s1(:,2),V2f(:,2),numSteps)';
Z=linspacen(V1s1(:,3),V2f(:,3),numSteps)';
[F1s,V1s,~]=patchCylSurfClose(X,Y,Z,[]);

%Derive "grid"
[V3s]=evenlySampleCurve(V3,size(V1s2,1),'pchip',1); %Resampling curve
[V3f]=minPolyTwist(V1s2,V3s); %Fix curve order
X=linspacen(V1s2(:,1),V3f(:,1),numSteps)';
Y=linspacen(V1s2(:,2),V3f(:,2),numSteps)';
Z=linspacen(V1s2(:,3),V3f(:,3),numSteps)';
[F2s,V2s,~]=patchCylSurfClose(X,Y,Z,[]);

%%

V=[V1s;V2s];
F=[F1s;F2s+size(V1s,1)];
faceMarker=[2*ones(size(F1s,1),1); 3*ones(size(F2s,1),1);];

[F,V]=mergeVertices(F,V);

%% FIND BOUNDARIES


boundEdges = patchBoundary(F,V);

epsMax=pointSpacing/2;%max(eps(V(:)))*10;
edgeGroups=tesgroup(boundEdges);
boundaryMarker=nan(size(V,1),1);
curveIndices=cell(1,3);
for q=1:1:size(edgeGroups,2)
    logicNow=edgeGroups(:,q)>0;
    edgesNow=boundEdges(logicNow,:);
    
    [indList]=edgeListToCurve(edgesNow);
    indList=indList(1:end-1);
    V_now=V(indList,:);
        
    mean_V_now=mean(V_now,1);
    D1=sqrt(sum((mean_V_now-meanV1).^2));
    D2=sqrt(sum((mean_V_now-meanV2).^2));
    D3=sqrt(sum((mean_V_now-meanV3).^2));
    if D1<epsMax
        curveIndices{1}=indList;
        boundaryMarker(indList)=1;
    elseif D2<epsMax
        curveIndices{2}=indList;
        boundaryMarker(indList)=2;
    elseif D3<epsMax
        curveIndices{3}=indList;
        boundaryMarker(indList)=3;
    end
end

%%

switch patchType
    case 'tri'
        F=Ft;
        faceMarker=faceMarker_t;
    case 'quad'
        
end

%% CONSTRAINED SMOOTHENING OF THE MESH

if ~isempty(smoothPar)
    %Adding point constaints
    boundaryInd=unique(boundEdges(:));
    smoothPar.RigidConstraints=boundaryInd; %Points to hold on to    
    V=patchSmooth(F,V,[],smoothPar);
end

%%

varargout{1}=F; 
varargout{2}=V; 
varargout{3}=curveIndices;
varargout{4}=faceMarker;
 
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
