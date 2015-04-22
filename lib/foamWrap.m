function [FT,VT,CT]=foamWrap(varargin)

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        C=(1:1:size(F,1))';        
        cPar=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};                
        cPar=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        cPar=varargin{4};
end

if isempty(cPar); %Use defaults
    cParSmooth.Method='HC';
    cParSmooth.n=50;
    cPar.n=4;
    cPar.Smooth=cParSmooth;
    cPar.dirFlip=1;
    [D]=patchEdgeLengths(F,V); %Edge lengths
    cPar.foamThickness=mean(D)/3; %1/3 of the mean edge length
end

if ~isfield(cPar,'n');
    cPar.n=4;
end

if ~isfield(cPar,'dirFlip');
    cPar.dirFlip=1;
end

if ~isfield(cPar,'Smooth');
    cParSmooth.Method='HC';
    cParSmooth.n=50;    
    cPar.Smooth=cParSmooth;
end

if ~isfield(cPar,'foamThickness');
    [D]=patchEdgeLengths(F,V); %Edge lengths
    cPar.foamThickness=mean(D)/3; %1/3 of the mean edge length
end

%% Get settings

n=cPar.n;
foamThickness=cPar.foamThickness; 
cParSmooth=cPar.Smooth;
dirFlip=cPar.dirFlip; 

%%

%Check if n can be achieved through splitting
nSplitIterations=log2(n+1); %Check for integer solution
logicInteger=abs(round(nSplitIterations)-nSplitIterations)<eps(nSplitIterations);

if logicInteger
    subMethod='split';
else
    subMethod='seed';
end

[Fs,Vs]=subtri(F,V,n,1);

%%
Ci=(1:size(F,1))';
switch  subMethod
    case 'split'
        Cs=repmat(Ci,[1,size(Fs,1)./size(F,1)]);
    case 'seed'
        Cs=Ci;
        Cs=Cs(:,ones(1,size(Fs,1)./size(F,1)));
        Cs=Cs';
        Cs=Cs(:);
end
Cs=C(Cs);

%%

[IND_Fs]=tesIND(Fs,Vs,0);
L=IND_Fs>0;     
Cv=nan(size(IND_Fs));
Cv(L)=Cs(IND_Fs(L));
Cv_max=nanmax(Cv,[],2);
Cv_min=nanmin(Cv,[],2);
logicMulti=~(Cv_max==Cv_min);

%%

%%
logicMultiFaces=any(logicMulti(Fs),2);

Fse=Fs(logicMultiFaces,:); 
Fsi=Fs(~logicMultiFaces,:); 
Fss=[Fse;Fsi];
Css=[ones(size(Fse,1),1);2*ones(size(Fsi,1),1)];


indVe=unique(Fse);
logicIndVe=false(size(Vs,1),1);
logicIndVe(indVe)=1; 

indAll=(1:1:size(Vs,1))';
indFix1=[indAll(logicIndVe);indAll(~logicIndVe)];
[~,indSort]=sort(indFix1);
Vss=Vs(indFix1,:);
Fss=indSort(Fss);

indStripVertices=1:numel(indVe);
indFaces=find(Css==1);

[IND_Fss]=tesIND(Fss,Vss,0);
L=IND_Fss>0;     
Cv=nan(size(IND_Fss));
Cv(L)=Css(IND_Fss(L));
Cv_max=nanmax(Cv,[],2);
Cv_min=nanmin(Cv,[],2);
logicMulti=~(Cv_max==Cv_min);
indBoundary=find(logicMulti);

%%
%Get vertex normals
[~,~,Nv]=patchNormal(Fss,Vss);

%Flip if desired
Nv=Nv.*dirFlip;

%Thicken inwards
VT2=Vss+foamThickness.*Nv;


E=patchEdges(Fss(indFaces,:),1);
logicBoundary=all(ismember(E,indBoundary),2);

Vn=[Vss; VT2(indStripVertices,:)];
Fn=[Fss; Fss(indFaces,:)+size(Vss,1)];
Cn=[Css; 1*ones(numel(indFaces),1)];
indBoundaryVertices1=indBoundary;
indBoundaryVertices2=indBoundary+size(Vss,1);
Eb1=E(logicBoundary,:);
Eb2=Eb1+size(Vss,1);

Fq=[Eb1 fliplr(Eb2)];

Ftn=[Fq(:,[1 2 3]); Fq(:,[3 4 1]);];

VT=Vn;
FT=[Fn;Ftn];
CT=[Cn;ones(size(Ftn,1),1);];

%%

indBlack=FT(CT==1,:);
indBlack=unique(indBlack(:));
L=true(size(VT,1),1);
L(indBlack)=0; 

indRigid=find(L);
indRigid=[indRigid;indBoundaryVertices1; indStripVertices(:)];
indBlack=unique(indRigid);

cParSmooth.RigidConstraints=indRigid;
[VT]=patchSmooth(FT,VT,[],cParSmooth);

%%
