
function [outputStruct]=import_obj(varargin)

% function [outputStruct]=import_obj(fileName,objOptionStruct)
% -----------------------------------------------------------------------
% This function imports wavefront OBJ files defined by fileName. A single
% body is assumed.
%
% Change log:
% 2021/04/30 KMM: Fixed bug relating to / symbols for faces
% 2021/04/30 KMM: Speeded up by using single cellfun based loop
% 2023/06/01 KMM: Improved speed, added handling of mixed face types
% -----------------------------------------------------------------------

%% Parse input

%Default option structure
defaultOptionStruct.fullMode=1;
defaultOptionStruct.colorMappingOpt='face-center';

switch nargin
    case 1
        fileName=varargin{1};
        optionStruct=defaultOptionStruct;
    case 2
        fileName=varargin{1};
        optionStruct=varargin{2};
end

%Complete input structure if needed
optionStruct=structComplete(optionStruct,defaultOptionStruct,1);

%Access input structure variables
fullMode=optionStruct.fullMode;
colorMappingOpt=optionStruct.colorMappingOpt;

%% Import and parse OBJ
% Import OBJ lines
S_OBJ = readlines(fileName);

%Process OBJ text to detect features
Lc=contains(S_OBJ,'# ');
L_mtl=contains(S_OBJ,'mtllib');
L_mat=contains(S_OBJ,'usemtl');

if fullMode %Also check for MTL and texture image
    fileName_mtl=sscanf(S_OBJ(L_mtl),'mtllib %s');
    %name_mat=sscanf(S_OBJ(L_mat),'usemtl %s');

    % Load texture image
    [filePath,~,~]=fileparts(fileName);
    S_MTL = readlines(fullfile(filePath,fileName_mtl));
    L_Kd = contains(S_MTL,'map_Kd');
    imageName_Kd=sscanf(S_MTL(L_Kd),'map_Kd %s');
    m=imread(fullfile(filePath,imageName_Kd));
end

L_exc = Lc | L_mtl | L_mat; %Logic of lines to exclude (all none v, vn, vt, or f lines).
Lv=contains(S_OBJ,'v ') & ~L_exc;
Lvn=contains(S_OBJ,'vn ') & ~L_exc;
Lvt=contains(S_OBJ,'vt ') & ~L_exc;
Lf=contains(S_OBJ,'f ') & ~L_exc;

%% Get vertices

nVertices = nnz(Lv);
V=zeros(nVertices,3);
s=S_OBJ(Lv);
for q=1:1:nVertices
    V(q,:)=sscanf(s(q),'v %f %f %f')';
end

%% Get normals/texture data

if fullMode
    nVertNorm = nnz(Lvn);
    nVertTexture = nnz(Lvt);
    N=zeros(nVertNorm,3);
    s=S_OBJ(Lvn);
    for q=1:1:nVertNorm
        N(q,:)=sscanf(s(q),'vn %f %f %f')';
    end

    T=zeros(nVertTexture,2);
    s=S_OBJ(Lvt);
    for q=1:1:nVertTexture
        T(q,:)=sscanf(s(q),'vt %f %f ')';
    end

    %Convert UV to image coordinates (subscript indices into array)
    ij_M=T; %Copy T
    ij_M=flip(ij_M,2); %Swap columns so they refer to rows then columns
    ij_M(:,1)=1-ij_M(:,1); %Invert the row coordinates
    ij_M(:,1)=(ij_M(:,1).*(size(m,1)-1))+1; %Scale to pixel coordinates
    ij_M(:,2)=(ij_M(:,2).*(size(m,2)-1))+1; %Scale to pixel coordinates
end

%% Get faces, and associated texture/normal indices if requested

% Allocate cell array for triangles (with vertex/normal data) for the moment (multiple face types may occur)
nFaces = nnz(Lf);
if fullMode==1
    FC=repmat({zeros(1,3)},nFaces,1);
    FD=repmat({zeros(1,3*3)},nFaces,1);
    FtC=FC;

    nf=zeros(nFaces,1);
    s=S_OBJ(Lf);
    sf=regexprep(s,'/+[0-9]+|f', ' '); %Remove f and get just face vertex indices
    sd=regexprep(s,'/|f', ' '); %Remove f and /
    for q=1:1:nFaces
        FC{q}=sscanf(sf(q),'%d')';
        FD{q}=sscanf(sd(q),'%d')';
        nf(q)=numel(FC{q});
        if size(FD{q},2)==2*nf(q)
            i=(1:nf(q))*2-1;
        elseif size(FD{q},2)==3*nf(q)
            i=(1:nf(q))*3-1;
        else
            disp('No way!')
        end
        FtC{q}=FD{q}(i);
    end

    faceTypeSet=unique(nf);
    numFaceTypes=numel(faceTypeSet);
    if numFaceTypes>1
        F=cell(numFaceTypes,1);
        F_uv=cell(numFaceTypes,1);
        for q=1:1:numFaceTypes
            logicNow = nf==faceTypeSet(q);
            F{q}=cell2mat(FC(logicNow));
            F_uv{q}=cell2mat(FtC(logicNow));
        end
    else
        logicNow = nf==faceTypeSet;
        F=cell2mat(FC(logicNow));
        F_uv=cell2mat(FtC(logicNow));
    end
    C=textureCoord2FaceColor(F_uv,ij_M,m,colorMappingOpt);
else
    FC=repmat({zeros(1,3)},nFaces,1);
    nf=zeros(nFaces,1);
    s=S_OBJ(Lf);
    sf=regexprep(s,'/+[0-9]+|f', ' '); %Remove f and get just face vertex indices
    for q=1:1:nFaces
        FC{q}=sscanf(sf(q),'%d')';
        nf(q)=numel(FC{q});
    end

    faceTypeSet=unique(nf);
    numFaceTypes=numel(faceTypeSet);
    if numFaceTypes>1
        F=cell(numFaceTypes,1);
        for q=1:1:numFaceTypes
            logicNow = nf==faceTypeSet(q);
            F{q}=cell2mat(FC(logicNow));
        end
    else % Just one face type
        logicNow = nf==faceTypeSet;
        F=cell2mat(FC(logicNow));
    end
end

%% Collect output

if fullMode==1
    outputStruct.F=F;
    outputStruct.V=V;
    outputStruct.C=C;
    outputStruct.F_uv=F_uv;
    outputStruct.ij_M=ij_M;
    outputStruct.m=m;
else
    outputStruct.F=F;
    outputStruct.V=V;
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
