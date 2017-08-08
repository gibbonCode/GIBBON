function [Es,Vs,Cs]=element2HexLattice(varargin)

% function [Es,Vs]=element2HexLattice(E,V,cPar)

%% Parse input

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        cPar=[];
    case 3
        E=varargin{1};
        V=varargin{2};
        cPar=varargin{3};
end

%The default control parameters
cParDefault.growSteps=1;
cParDefault.latticeSide=1;
cParDefault.indBoundary=[];

if isempty(cPar)
    cPar=cParDefault;
else
    if ~isfield(cPar,'growSteps')
        cPar.growSteps=cParDefault.growSteps;
    end
    if ~isfield(cPar,'latticeSide')
        cPar.latticeSide=cParDefault.latticeSide;
    end
end

%%

switch size(E,2)
    case 4
        elementType='tet4';
    case 8
        elementType='hex8';
    otherwise
        error('Element type not supported');
end

%%

switch elementType
    case 'tet4'
        [Eh,Vh,CVh]=tet2hex(E,V,1);
        ind=find(CVh==2 | CVh==3);
    case 'hex8'
        [Eh,Vh,~,CVh]=subHex(E,V,1,1);        
        ind=find(CVh==2 | CVh==3);
end
[Eh,Vh]=subHex(Eh,Vh,1,1);

LE=~any(ismember(Eh,ind),2);

%%

if cPar.growSteps<0
    LE=~LE;
end

for q=1:1:abs(cPar.growSteps)
    [Eh,Vh]=subHex(Eh,Vh,1,1);
    LE=repmat(LE,[8,1]);
    
    ind_LE=Eh(LE,:);
    ind_LE=unique(ind_LE(:));
    LE=any(ismember(Eh,ind_LE),2);
end

if cPar.growSteps<0
    LE=~LE;
end
    
%%


if isempty(cPar.latticeSide)
    Es=[Eh(LE==0,:);Eh(LE==1,:)];
    Cs=[ones(nnz(LE==0),1); 2*ones(nnz(LE==1),1)];
elseif cPar.latticeSide==1
    Es=Eh(LE==1,:);
    Cs=ones(size(Es,1),1);
elseif cPar.latticeSide==2
    Es=Eh(LE==0,:);
    Cs=ones(size(Es,1),1);
end

Vs=Vh;
[Es,Vs]=patchCleanUnused(Es,Vs);




 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
