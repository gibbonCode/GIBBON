%% import_obj
% Below is a demonstration of the features of the |import_obj| function

%%
clear; close all; clc;

%% Syntax
% |[F,V,C]=import_obj(fileName);|

%% Description 
% This function imports the geometry (faces and nodes0 contained in an OBJ
% file. All texture/material data is ignored. 

%% Examples 
% 

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
loadPath=fullfile(defaultFolder,'data','OBJ');
fileName=fullfile(loadPath,'gibbon.obj');

%%

[filePath,~,~]=fileparts(fileName);

S_OBJ = readlines(fileName);
Lc=contains(S_OBJ,'# ');
L_mtl=contains(S_OBJ,'mtllib');
fileName_mtl=sscanf(S_OBJ(L_mtl),'mtllib %s');
L_mat=contains(S_OBJ,'usemtl');
name_mat=sscanf(S_OBJ(L_mat),'usemtl %s');

S_MTL = readlines(fullfile(filePath,fileName_mtl));
L_Kd = contains(S_MTL,'map_Kd');
imageName_Kd=sscanf(S_MTL(L_Kd),'map_Kd %s');

m=imread(fullfile(filePath,imageName_Kd)); 
m=permute(m,[2 1 3]);

cFigure; hold on;  
title('Texture image')
image(m);

M=double(m)./255;

L_exc = Lc | L_mtl | L_mat; 

Lv=contains(S_OBJ,'v ') & ~L_exc;
Lvn=contains(S_OBJ,'vn ') & ~L_exc;
Lvt=contains(S_OBJ,'vt ') & ~L_exc;
Lf=contains(S_OBJ,'f ') & ~L_exc;

nVertices = nnz(Lv);
nVertNorm = nnz(Lvn);
nVertTexture = nnz(Lvt);
nFaces = nnz(Lf);

V=zeros(nVertices,3);
s=S_OBJ(Lv);
for q=1:1:nVertices
    V(q,:)=sscanf(s(q),'v %f %f %f')';
end

N=zeros(nVertices,3);
s=S_OBJ(Lvn);
for q=1:1:nVertices
    N(q,:)=sscanf(s(q),'vn %f %f %f')';
end

T=zeros(nVertTexture,2);
s=S_OBJ(Lvt);
for q=1:1:nVertTexture
    T(q,:)=sscanf(s(q),'vt %f %f ')';
end
ij_M=T; 
ij_M(:,1)=round((ij_M(:,1).*(size(M,1)-1))+1); 
ij_M(:,2)=round((ij_M(:,2).*(size(M,2)-1))+1);

%Allocate cell array for triangles (with vertex/normal data) for the moment (multiple face types may occur)
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
F=cell(numFaceTypes,1);
C=cell(numFaceTypes,1);
for q=1:1:numFaceTypes
    logicNow = nf==faceTypeSet(q);
    F{q}=cell2mat(FC(logicNow));

    ft=cell2mat(FtC(logicNow));

    RGB_now=zeros(size(ft,1),3);
    for qc=1:1:size(ft,2)
        ij_M_now=ij_M(ft(:,qc),:);
        ind_M_now_R=sub2ind(size(M),ij_M(ft(:,qc),1),ij_M(ft(:,qc),2),1*ones(size(ft(:,qc))));
        ind_M_now_G=sub2ind(size(M),ij_M(ft(:,qc),1),ij_M(ft(:,qc),2),2*ones(size(ft(:,qc))));
        ind_M_now_B=sub2ind(size(M),ij_M(ft(:,qc),1),ij_M(ft(:,qc),2),3*ones(size(ft(:,qc))));        
        RGB_now=RGB_now+[M(ind_M_now_R) M(ind_M_now_G) M(ind_M_now_B)]/3;
    end
    C{q}=RGB_now;
end

%%

warning('Texture is still wrong')
cFigure; 
gpatch(F,V,C,'k')
% patchNormPlot(F,V);
axisGeom; camlight headlight; 
gdrawnow; 

% F
% nf



%S(Lf)

%%
% [F,V]=import_obj_geom(fileName); 

% %%
% % Visualize
% 
% cFigure; 
% gpatch(F,V,'w','k');
% axisGeom; camlight headlight; 
% gdrawnow;

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
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
