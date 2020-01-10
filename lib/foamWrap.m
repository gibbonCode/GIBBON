function [varargout]=foamWrap(varargin)

% function [FT,VT,CT,CT_c]=foamWrap(F,V,C,cPar)
% ------------------------------------------------------------------------
%
% Adds a foam like structural offset from input mesh. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015
% 2016/04/19: Updated with surface color book keeping
%------------------------------------------------------------------------
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

if ~isfield(cPar,'n')
    cPar.n=4;
elseif isempty(cPar.n);
    cPar.n=4;
end

if ~isfield(cPar,'dirFlip');
    cPar.dirFlip=1;
elseif isempty(cPar.dirFlip);    
    cPar.dirFlip=1;
end

if ~isfield(cPar,'Smooth');
    cParSmooth.Method='HC';
    cParSmooth.n=50;    
    cPar.Smooth=cParSmooth;
elseif  isempty(cPar.cParSmooth);
end

if ~isfield(cPar,'foamThickness');
    [D]=patchEdgeLengths(F,V); %Edge lengths
    cPar.foamThickness=mean(D)/2; %1/2 of the mean edge length
elseif  isempty(cPar.foamThickness);
    [D]=patchEdgeLengths(F,V); %Edge lengths
    cPar.foamThickness=mean(D)/2; %1/2 of the mean edge length
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
Cv_max=gnanmax(Cv,[],2);
Cv_min=gnanmin(Cv,[],2);
logicMulti=~(Cv_max==Cv_min);

%%

% XX=Vs(:,3);
% XX_F=mean(XX(Fs),2);
% logicLowFaces=XX_F<30; 
% lowColors=unique(Cs(logicLowFaces));
% logicLowFaceCells=ismember(Cs,lowColors);
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
Cv_max=gnanmax(Cv,[],2);
Cv_min=gnanmin(Cv,[],2);
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
Cn_c=[Css; 4*ones(numel(indFaces),1)];
indBoundaryVertices1=indBoundary;
indBoundaryVertices2=indBoundary+size(Vss,1);
Eb1=E(logicBoundary,:);
Eb2=Eb1+size(Vss,1);

Fq=[Eb1 fliplr(Eb2)];

Ftn=[Fq(:,[1 2 3]); Fq(:,[3 4 1]);];

VT=Vn;
FT=[Fn;Ftn];
CT=[Cn;ones(size(Ftn,1),1);];
CT_c=[Cn_c; 3*ones(size(Ftn,1),1);];
CT_cc=CT_c;
CT_c(CT_cc==1)=2; 
CT_c(CT_cc==2)=1; 

%%

indBlack=FT(CT==1,:);
indBlack=unique(indBlack(:));
L=true(size(VT,1),1);
L(indBlack)=0; 

indRigid=find(L);
indRigid=[indRigid;indBoundaryVertices1; indStripVertices(:)];
indBlack=unique(indRigid);

cParSmooth.RigidConstraints=indBlack;
[VT]=patchSmooth(FT,VT,[],cParSmooth);

%%

varargout{1}=FT; %Faces 
varargout{2}=VT; %Vertices
varargout{3}=CT; %Two colors
varargout{4}=CT_c; %Four colors 1=input, 2=base for offset, 3=sides for offset, 4=top for offset

 
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
