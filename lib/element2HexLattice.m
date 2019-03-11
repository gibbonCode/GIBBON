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
