function [Fc_Q,Fc_T,Ft_Q,Ft_T,Ct_Q,Ct_T,Vt]=rhombicDodecahedronMesh(r,nCopies)

% function [Fc_Q,Fc_T,Ft_Q,Ft_T,Ct_Q,Ct_T,Vt]=rhombicDodecahedronMesh(r,nCopies)
% ------------------------------------------------------------------------
% Creates a rhombic dodecahedron mesh where r sets the radias and nCopies
% (a 1x3 vector) sets the number of copies in the x, y, and z direction.
% The output consists of:
%
% Fc_Q, Fc_T: the quadrilateral and triangular face cell arrays (1 cell
% entry per element). 
%
% Ft_Q, Ft,T: the quadrilateral and triangular face arrays
%
% Ct_Q, Ct,T: color/label data for the face arrays
%
% Vt: the vertex array
%
% 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2019/02/08 Created
% ------------------------------------------------------------------------

%%

%Get rhombic dodecahedron
[Fs,Vs]=rhombicDodecahedron(r);
Fst=[Fs(:,[1 2 3]); Fs(:,[3 4 1])]; %Triangular faces

%Rotate
[R,~]=euler2DCM([0 0 0.25*pi]); %Define rotation
Vs=Vs*R; %Rotate coordinates

%Derive face centre points for offsets
Xi=Vs(:,1); Yi=Vs(:,2); Zi=Vs(:,3);
Vn=[mean(Xi(Fs),2) mean(Yi(Fs),2) mean(Zi(Fs),2)];

%%
nTotal=prod(nCopies); %Total number of copies
offsetDirs=[3 4 6]; %Offset direction N.B. varying these affects the offsets/signs below

%Create cell indices for vertices
indC=ones(size(Vs,1),1)*(1:1:nTotal);
indC=indC(:);
[I,J,K] = ind2sub(nCopies,indC);
I=I-1; J=J-1; K=K-1;

%Create cell indices for quad faces
indF_Q=ones(size(Fs,1),1)*(1:1:nTotal);
indF_Q=indF_Q(:);
indF_Q=(indF_Q-1);

%Defining the quad faces matrix
Ft_Q=repmat(Fs,nTotal,1)+size(Vs,1).*indF_Q(:,ones(1,size(Fs,2))); 

%Create cell indices for triangular faces
indF_T=ones(size(Fst,1),1)*(1:1:nTotal);
indF_T=indF_T(:);
indF_T=(indF_T-1);

%Defining the tri faces matrix
Ft_T=repmat(Fst,nTotal,1)+size(Vs,1).*indF_T(:,ones(1,size(Fst,2))); 

%Defining offsets
sK=~iseven(K); %Shift is adjusted according to z coordinate to create a "cube"
D1=I*2*Vn(offsetDirs(1),:)-(K*Vn(offsetDirs(1),:))+(sK*Vn(offsetDirs(1),:)); % X offsets
D2=J*2*Vn(offsetDirs(2),:)+(K*Vn(offsetDirs(2),:))+(sK*Vn(offsetDirs(2),:)); % Y offsets
D3=K*2*Vn(offsetDirs(3),:); % Z offsets

%Defining vertices matrix
Vt=repmat(Vs,nTotal,1)+(D1+D2+D3);

%Merging points
[~,IND_V,IND_IND]=unique(round(Vt.*(1e5)),'rows');
Vt=Vt(IND_V,:);
Ft_Q=IND_IND(Ft_Q);
Ft_T=IND_IND(Ft_T);

Ct_T=indF_T+1; %Index or color number
Ct_Q=indF_Q+1; %Index or color number

%% Creating cell type output

%Split up face matrix in to cell groups
Fc_Q=mat2cell(Ft_Q,size(Fs,1)*ones(1,nTotal),size(Fs,2));
Fc_T=mat2cell(Ft_T,size(Fst,1)*ones(1,nTotal),size(Fst,2));

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
