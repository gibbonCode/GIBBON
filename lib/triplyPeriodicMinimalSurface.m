function [varargout]=triplyPeriodicMinimalSurface(inputStruct)

% -----------------------------------------------------------------------
% This function generates a triply periodic minimal surface
%
% Input structure and default values:
%
% inputStruct.L=1; % characteristic length
% inputStruct.Ns=80; % number of sampling points
% inputStruct.isocap=1; %Option to cap the isosurface
% inputStruct.surfaceCase='g'; %Surface type
% inputStruct.numPeriods=[1 1 1]; %Number of periods in each direction
% inputStruct.levelset=0.5; %Isosurface level
%
% Change log: 
% 2020/11/05: Created
% -----------------------------------------------------------------------

%% Parse input

%Create default structure
defaultInputStruct.L=1; % characteristic length
defaultInputStruct.Ns=80; % number of sampling points
defaultInputStruct.isocap=1; %Option to cap the isosurface
defaultInputStruct.surfaceCase='g';
defaultInputStruct.numPeriods=[1 1 1]; 
defaultInputStruct.levelset=0.5;
defaultInputStruct.surfaceSide=1;

%Complete input with default if incomplete
[inputStruct]=structComplete(inputStruct,defaultInputStruct,1); %Complement provided with default if missing or empty

%Get parameters from input structure
L = inputStruct.L; % characteristic length
Ns = inputStruct.Ns; % number of sampling points
isocap= inputStruct.isocap; 
surfaceCase=inputStruct.surfaceCase; 
numPeriods=inputStruct.numPeriods; 
levelset=inputStruct.levelset;

if numel(numPeriods)==1
    numPeriods=numPeriods*ones(1,3);
end

%% Generation of points

%Define coordinate limits
xMin=0;
xMax=xMin+2*pi*numPeriods(1);
yMin=0;
yMax=yMin+2*pi*numPeriods(2);
zMin=0;
zMax=zMin+2*pi*numPeriods(3);

%Create coordinates
xRange=linspace(xMin,xMax,Ns);
yRange=linspace(yMin,yMax,Ns);
zRange=linspace(zMin,zMax,Ns);
[X,Y,Z]=meshgrid(xRange,yRange,zRange);

%Calculate 3D image data
S=triplyPeriodicMinimal(X(:),Y(:),Z(:),surfaceCase);        
S=reshape(S,size(X));

%Scaling coordinates
X=((X./abs(xMax-xMin)).*L);
Y=((Y./abs(yMax-yMin)).*L);
Z=((Z./abs(zMax-zMin)).*L);

%% Generation of level sets and exporting

%Iso-surface
switch inputStruct.surfaceSide
    case 0
        [f,v,c]=getSurface(X,Y,Z,S,levelset);
        if isocap==1
            [fc1,vc1,cc1]=getCaps(X,Y,Z,S,levelset);
            [fc2,vc2,cc2]=getCaps(X,Y,Z,-S,-levelset);
            [f,v,c]=joinElementSets({f,fc1,fc2},{v,vc1,vc2},{c,cc1,cc2+6});
        end
    case 1 
        [f,v,c]=getSurface(X,Y,Z,S,levelset);        
        if isocap==1
            [fc,vc,cc]=getCaps(X,Y,Z,S,levelset);
        end
        %Join sets
        [f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc});
    case -1
        [f,v,c]=getSurface(X,Y,Z,-S,-levelset);        
        if isocap==1 
            [fc,vc,cc]=getCaps(X,Y,Z,-S,-levelset);
        end
        %Join sets
        [f,v,c]=joinElementSets({f,fc},{v,vc},{c,cc});
end


%% Merge nodes and clean-up mesh 

%Merge nodes
[f,v]=mergeVertices(f,v); 

%Check for unique faces
[~,indUni,~]=unique(sort(f,2),'rows');
f=f(indUni,:); %Keep unique faces
c=c(indUni);

%Remove collapsed faces
[f,logicKeep]=patchRemoveCollapsed(f); 
c=c(logicKeep);

%Remove unused points
[f,v]=patchCleanUnused(f,v); 

%Invert faces
f=fliplr(f); 

%% Collect output
varargout{1}=f;
varargout{2}=v;
varargout{3}=c;
varargout{4}=S;

end
%%

function [f,v,c]=getSurface(X,Y,Z,S,levelset)
    [f,v] = isosurface(X,Y,Z,S,levelset);
    c=zeros(size(f,1),1);
end

function [fc,vc,cc]=getCaps(X,Y,Z,S,levelset)
    [fc,vc] = isocaps(X,Y,Z,S,levelset);     %Compute isocaps
    
    nc=patchNormal(fc,vc);
    cc=zeros(size(fc,1),1);
    cc(nc(:,1)<-0.5)=1;
    cc(nc(:,1)>0.5)=2;
    cc(nc(:,2)<-0.5)=3;
    cc(nc(:,2)>0.5)=4;
    cc(nc(:,3)<-0.5)=5;
    cc(nc(:,3)>0.5)=6;    
end

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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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


