function [varargout]=triBox(varargin)

% function [F,V,faceBoundaryMarker]=triBox(boxDim,pointSpacing,numDigitsMerge)
%-------------------------------------------------------------------------
%
%
%-------------------------------------------------------------------------

%%

switch nargin
    case 1
        boxDim=varargin{1};
        pointSpacing=mean(boxDim)/10;
        numDigitsMerge=[];
    case 2
        boxDim=varargin{1};
        pointSpacing=varargin{2};
        numDigitsMerge=[];
    case 3
        boxDim=varargin{1};
        pointSpacing=varargin{2};
        numDigitsMerge=varargin{3};
end

if isempty(numDigitsMerge)    
   numDigitsMerge=6-numOrder(pointSpacing);
end

%% Create quadrilateral box

boxEl=round(boxDim./pointSpacing); %Number of elements per direction 
[Fq,Vq,Cq]=quadBox(boxDim,boxEl); %Create quadrilateral box

%% Remesh 3 side face groups with triangles

%Mesh 3 side faces
for q=[1 3 5]        
    [Eb]=patchBoundary(Fq(Cq==q,:),Vq); %Get boundary of side face set
    [indList]=edgeListToCurve(Eb); %Indiced defining closed curve       
    Vc=Vq(indList(1:end-1),:); %Points defining closed curve
    
    switch q
        case 1 %Side face 1
            [Ft,Vt]=regionTriMesh2D({Vc(:,[2 3])},pointSpacing,0,0);
            Vt=[zeros(size(Vt,1),1) Vt(:,[1 2])];
            Vt(:,1)=mean(Vc(:,1));
            Ft=fliplr(Ft);
            F1=Ft; V1=Vt;
        case 3 %Side face 3
            [Ft,Vt]=regionTriMesh2D({Vc(:,[1 3])},pointSpacing,0,0);
            Vt=[Vt(:,1) zeros(size(Vt,1),1) Vt(:,2)];
            Vt(:,2)=mean(Vc(:,2));            
            F3=Ft; V3=Vt;
        case 5 %Side face set 5
            [Ft,Vt]=regionTriMesh2D({Vc(:,[1 2])},pointSpacing,0,0);            
            Vt(:,3)=mean(Vc(:,3));
            Ft=fliplr(Ft);
            F5=Ft; V5=Vt;
    end    
end

%% Copy, invert faces, and offset to generate the other 3

F2=fliplr(F1);
V2=V1; 
V2(:,1)=max(Vq(:,1));

F4=fliplr(F3);
V4=V3; 
V4(:,2)=max(Vq(:,2));

F6=fliplr(F5);
V6=V5; 
V6(:,3)=max(Vq(:,3));

%% Join face sets
[F,V,faceBoundaryMarker]=joinElementSets({F1,F2,F3,F4,F5,F6},{V1,V2,V3,V4,V5,V6});

%% Merge face sets

[F,V]=mergeVertices(F,V,numDigitsMerge);

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=faceBoundaryMarker;

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
