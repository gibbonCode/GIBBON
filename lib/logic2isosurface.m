function [varargout]=logic2isosurface(L,controlPar)

% function [F,V,C]=logic2isosurface(L,controlPar)
% ------------------------------------------------------------------------
% 
%
% Kevin Moerman
%
% Change log: 
% 2019/05/01, Kevin Moerman: Added varargout support
% 2019/05/01, Kevin Moerman: Added colordata export option so caps can be retrieved
% 2019/05/01, Kevin Moerman: Minor updates to comments and code organisation
% ------------------------------------------------------------------------

%% Parse input

if isfield(controlPar,'kernelSize')    
    kernelSize=controlPar.kernelSize;
else %DEFAULT
    kernelSize=3; 
end

if isfield(controlPar,'contourLevel')    
    contourLevel=controlPar.contourLevel;
else %DEFAULT
    contourLevel=0.5; 
end

if isfield(controlPar,'voxelSize')    
    voxelSize=controlPar.voxelSize;   
else %DEFAULT
    voxelSize=[1 1 1];
end

if isfield(controlPar,'capOpt')    
    capOpt=controlPar.capOpt;
else %DEFAULT
    capOpt=1;
end

if isfield(controlPar,'nSub')    
    nSub=controlPar.nSub;
else %DEFAULT
    nSub=[1 1 1];
end

if isfield(controlPar,'kernelType')    
    kernelType=controlPar.kernelType;
else %DEFAULT
    kernelType=1;
end

siz=size(L);

%% Smoothen logic

if ~isnan(kernelSize) && ~isnan(kernelType)
    %Define smoothening kernel
    switch kernelType
        case 1 %Simple averaging kernel
            hg=ones(kernelSize,kernelSize,kernelSize);
            hg=hg./sum(hg(:));
        case 2
            hg=gauss_kernel(kernelSize,ndims(L),2,'width');
    end
    
    %Convolve logic
    L=convn(double(L),hg,'same');
else
    L=double(L);
end

%% CREATE ISO-SURFACE

%Pad with zeros if desired
if capOpt==1
    Lpad=zeros(size(L)+2*kernelSize);
    Lpad(1+kernelSize:siz(1)+kernelSize,1+kernelSize:siz(2)+kernelSize,1+kernelSize:siz(3)+kernelSize)=L;
    L=Lpad;
    siz=size(L);
end

% Get image coordinates
[J,I,K]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3));
if capOpt==1
    I=I-kernelSize; J=J-kernelSize; K=K-kernelSize; %Correct for padding
end

%Convert to Cartesian coordinates using voxels size if provided
[X,Y,Z]=im2cart(I,J,K,voxelSize);

%Resample mesh and image
nNew=round(siz./nSub);
useRange_I=unique(round(linspace(1,siz(1),nNew(1))));
useRange_J=unique(round(linspace(1,siz(2),nNew(2))));
useRange_K=unique(round(linspace(1,siz(3),nNew(3))));

X_iso=X(useRange_I,useRange_J,useRange_K);
Y_iso=Y(useRange_I,useRange_J,useRange_K);
Z_iso=Z(useRange_I,useRange_J,useRange_K);
L_iso=L(useRange_I,useRange_J,useRange_K);

% X_iso=X(1:nSub(1):end,1:nSub(2):end,1:nSub(3):end);
% Y_iso=Y(1:nSub(1):end,1:nSub(2):end,1:nSub(3):end);
% Z_iso=Z(1:nSub(1):end,1:nSub(2):end,1:nSub(3):end);
% L_iso=L(1:nSub(1):end,1:nSub(2):end,1:nSub(3):end);

%Derive isosurface
[F,V] = isosurface(X_iso,Y_iso,Z_iso,L_iso,contourLevel);
F=F(:,[3 2 1]); %Flip face order so normal is outward
C=ones(size(F,1),1); %Add color data

%%
%Derive caps
if capOpt==2
    [Fc,Vc] = isocaps(X_iso,Y_iso,Z_iso,L_iso,contourLevel); %Get caps
    if ~isempty(Fc)
        Fc=fliplr(Fc); %Flip face order so normal is outward
        [F,V,C]=joinElementSets({F,Fc},{V,Vc},{C,2*ones(size(Fc,1),1)}); %Join sets
        [F,V]=mergeVertices(F,V); %merge isosurface and cap nodes
    end
end

%% Other improvements

[F,V,C]=triSurfRemoveThreeConnect(F,V,C);
[F,V]=mergeVertices(F,V);
F=patchRemoveCollapsed(F);
[F,V]=patchCleanUnused(F,V);

%%
varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
 
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
