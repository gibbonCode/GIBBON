function [varargout]=coneTriMesh(varargin)

% function [F,V,C]=coneTriMesh(coneRadius,coneHeight,pointSpacing,nTopMin,closeBaseOpt)
% ------------------------------------------------------------------------
% The |coneTriMesh| function creates a triangulated mesh for a cone with a
% base radius coneRadius and a height coneHeight. The desired node spacing
% is set by pointSpacing. The parameter nTopMin sets the number of repeated
% "pie-segments" used to build the cone, and hence also sets the minimum
% number of points the tip node is connected to. 
%
% 2023/06/10 KMM: Created
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 0
        coneRadius=1;
        coneHeight=1;
        pointSpacing=coneRadius/10;
        nTopMin=6;
        closeBaseOpt=0;
    case 1
        coneRadius=varargin{1};
        coneHeight=1;
        pointSpacing=coneRadius/10;
        nTopMin=6;
        closeBaseOpt=0;
    case 2
        coneRadius=varargin{1};
        coneHeight=varargin{2};
        pointSpacing=coneRadius/10;
        nTopMin=6;
        closeBaseOpt=0;
    case 3
        coneRadius=varargin{1};
        coneHeight=varargin{2};
        pointSpacing=varargin{3};
        nTopMin=6;
        closeBaseOpt=0;
    case 4
        coneRadius=varargin{1};
        coneHeight=varargin{2};
        pointSpacing=varargin{3};
        nTopMin=varargin{4};
        closeBaseOpt=0;
    case 5
        coneRadius=varargin{1};
        coneHeight=varargin{2};
        pointSpacing=varargin{3};
        nTopMin=varargin{4};
        closeBaseOpt=varargin{5};
end

if isempty(pointSpacing)
    pointSpacing=coneRadius/10;
end

%%


l = sqrt(coneRadius^2+coneHeight^2); 

thetaMax=2*pi/(l/coneRadius);
thetaMaxSplit=thetaMax/nTopMin; 

n=ceil((thetaMaxSplit*l)/pointSpacing);

t = linspace(0,thetaMaxSplit,n)'; 

Vc = l.*[cos(t) sin(t)];

nr=ceil(l/pointSpacing);
V1 = [linspace(0,l,nr)' zeros(nr,1)];
V2 = [linspace(Vc(end,1),0,nr)' linspace(Vc(end,2),0,nr)'];
Vm = [0 0];

Vt=[V1(1:end-1,:); Vc; V2(2:end-1,:)];

inputStructure.regionCell={Vt};
inputStructure.pointSpacing=pointSpacing;
inputStructure.resampleCurveOpt=0;
inputStructure.plotOn=0;
% inputStructure.mustPointsInner=Vm;
[F,V]=regionTriMesh2D(inputStructure);
V(:,3)=0;

Fr=cell(nTopMin,1);
Vr=cell(nTopMin,1);
for q=0:1:nTopMin-1
    R=euler2DCM([0 0 -q*thetaMaxSplit]);    
    Vr{q+1}=V*R;
    Fr{q+1}=F;
end
[F,V]=joinElementSets(Fr,Vr);

V=V(:,[1 2]);

[F,V]=mergeVertices(F,V); 
Eb=patchBoundary(F);

optionStructSmooth.n=25;
optionStructSmooth.Method='LAP';
optionStructSmooth.RigidConstraints=unique(Eb(:));
optionStructSmooth.Tolerance=pointSpacing/1000;
[V]=patchSmooth(F,V,[],optionStructSmooth);

[thV,rV]=cart2pol(V(:,1),V(:,2));
rV=coneRadius.*(rV./l);
thV(thV<0)=2*pi-abs(thV(thV<0));
thV=2*pi*(thV./max(thV));

V2=V;
[V2(:,1),V2(:,2)]=pol2cart(thV,rV);

V2(:,3) = coneHeight.* (coneRadius-sqrt(sum(V2.^2,2)))./coneRadius;

F2=F;
[F2,V2,~,ind2]=mergeVertices(F2,V2); 

%%

Eb1=ind2(Eb);
Eb2=patchBoundary(F2);
[~,indMax]=max(V2(:,3));

indAll=1:1:size(V2,1);
indBoundaryTop=[unique(Eb2(:)); indMax];
indZipNotBoundaryTop=Eb1(~ismember(Eb1,indBoundaryTop));
indFix=[find(~ismember(indAll,indZipNotBoundaryTop)) indBoundaryTop'];

clear optionStructSmooth
optionStructSmooth.n=1;
optionStructSmooth.Method='HC';
optionStructSmooth.RigidConstraints=indFix;
optionStructSmooth.Tolerance=pointSpacing/1000;

for q=1:1:5
    [V2]=patchSmooth(F2,V2,[],optionStructSmooth);
    [thV,~,zV]=cart2pol(V2(:,1),V2(:,2),V2(:,3));
    rV=coneRadius-((coneRadius/coneHeight).*zV);
    V2(:,1)=rV.*cos(thV);
    V2(:,2)=rV.*sin(thV);
    V2(:,3)=zV;
    [V2(:,1),V2(:,2),V2(:,3)]=pol2cart(thV,rV,zV);
end

%%

if closeBaseOpt
    Eb=patchBoundary(F2);
    indBoundaryList=edgeListToCurve(Eb); indBoundaryList=indBoundaryList(1:end-1);
    [Fb,Vb]=regionTriMesh2D({V2(indBoundaryList,[1 2])},pointSpacing,0,0);
    Fb=fliplr(Fb);
    Vb(:,3)=0;
    [F2,V2,C2]=joinElementSets({F2,Fb},{V2,Vb});
    [F2,V2]=mergeVertices(F2,V2);
    
else
    C2=ones(size(F2,1),1); 
end

%% 

varargout{1}=F2;
varargout{2}=V2;
varargout{3}=C2;

end

