function [F1c,F2c,F1,F2,C1,C2,VT]=truncatedOctahedronMesh(varargin)

% function [F1c,F2c,F1,F2,C1,C2,VT]=truncatedOctahedronMesh(r,nCopies)
% ------------------------------------------------------------------------
%
% 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2019/02/08 Created
% 2019/10/13 
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        r=varargin{1};
        nCopies=2;
    case 2
        r=varargin{1};
        nCopies=varargin{2};
end

if isempty(nCopies)
    nCopies=2;
end

if numel(nCopies)==1
    nCopies=nCopies*ones(1,3);
end

%%
%Get truncated octahedron
[Fs,Vs,~]=truncatedOctahedron(r);
Fs1=Fs{1}; %hexagon faces
Fs2=Fs{2}; %quadrilateral faces

%Rotate
[R,~]=euler2DCM([0 0 0.25*pi]); %Define rotation
Vs=Vs*R; %Rotate coordinates

%Derive face centre points for offsets
% Vn=patchCentre(Fs,Vs);
Xi=Vs(:,1); Yi=Vs(:,2); Zi=Vs(:,3);
Vn=[mean(Xi(Fs2),2) mean(Yi(Fs2),2) mean(Zi(Fs2),2)];

%%
nTotal=prod(nCopies); %Total number of copies
offsetDirs=[1 2 6]; %Offset direction N.B. varying these affects the offsets/signs below

%Create cell indices for vertices
indC=ones(size(Vs,1),1)*(1:1:nTotal);
indC=indC(:);
[I,J,K] = ind2sub(nCopies,indC);
I=I-1; J=J-1; K=K-1;

%Create cell indices for quad faces
indF1=ones(size(Fs1,1),1)*(1:1:nTotal);
indF1=indF1(:);
indF1=(indF1-1);

indF2=ones(size(Fs2,1),1)*(1:1:nTotal);
indF2=indF2(:);
indF2=(indF2-1);

%Defining the quad faces matrix
F1=repmat(Fs1,nTotal,1)+size(Vs,1).*indF1(:,ones(1,size(Fs1,2))); 
F2=repmat(Fs2,nTotal,1)+size(Vs,1).*indF2(:,ones(1,size(Fs2,2))); 

%Defining offsets
sK=iseven(K); %Shift is adjusted according to z coordinate to create a "cube"
D1=I*2*Vn(offsetDirs(1),:)+(sK*Vn(offsetDirs(1),:)); % X offsets
D2=J*2*Vn(offsetDirs(2),:)+(sK*Vn(offsetDirs(2),:)); % Y offsets
D3=K*Vn(offsetDirs(3),:); % Z offsets
 
D=D1+D2+D3;

%Defining vertices matrix
VT=repmat(Vs,nTotal,1)+D;

%Merge points
[F1,VT,~,ind2]=mergeVertices(F1,VT); 
F2=ind2(F2);

C1=indF1+1; %Index or color number
C2=indF2+1; %Index or color number

%Split up face matrix in to cell groups
F1c=mat2cell(F1,size(Fs1,1)*ones(1,nTotal),size(Fs1,2));
F2c=mat2cell(F2,size(Fs2,1)*ones(1,nTotal),size(Fs2,2));

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
