function [Fc_Q,Fc_T,Ft_Q,Ft_T,Ct_Q,Ct_T,Vt]=rhombicDodecahedronMesh(r,nCopies)

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

%Defining the quad faces matrix
Ft_T=repmat(Fst,nTotal,1)+size(Vs,1).*indF_T(:,ones(1,size(Fst,2))); 

%Defining offsets
sK=~iseven(K); %Shift is adjusted according to z coordinate to create a "cube"
D1=I*2*Vn(offsetDirs(1),:)-(K*Vn(offsetDirs(1),:))+(sK*Vn(offsetDirs(1),:)); % X offsets
D2=J*2*Vn(offsetDirs(2),:)+(K*Vn(offsetDirs(2),:))+(sK*Vn(offsetDirs(2),:)); % Y offsets
D3=K*2*Vn(offsetDirs(3),:); % Z offsets

%Defining vertices matrix
Vt=repmat(Vs,nTotal,1)+(D1+D2+D3);

%% REMOVING DOUBLE POINTS
%N.B. removal of points here based on rounding!
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
