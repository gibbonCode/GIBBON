function [E,Ev,Fe,Ve]=patch2edgeIND(F,V)

%% DERIVE EDGE MATRIX

% E=[]; %Can be improved through memory allocation and indexing in loop
% for i=1:1:size(F,2)-1;
%   E=[E; F(:,i) F(:,i+1)]; 
% end
% E=[E; F(:,end) F(:,1)]; 

%Format of column index in F
EColumnInd=[(1:size(F,2)); (1:size(F,2))];
EColumnInd=[EColumnInd(2:end) EColumnInd(1)];

%Derive edges matrix
E=F(:,EColumnInd)'; %Use index into F to create edges matrix
E=reshape(E,2,numel(E)/2)'; 

%%
E=sort(E,2); %Sort edge order
[E,~,ind2] = unique(E,'rows'); %Removing double edges, i.e. [1  4] = [4  1]

%%

Fe=reshape(1:numel(F),size(F,2),size(F,1))';
Fe=ind2(Fe);

ind_E=(1:size(E,1))'*ones(1,2);

Ve=sort((sparse(E(:),ind_E(:),ind_E(:),size(V,1),size(E,1))),2);
[~,J]=find(Ve);
Ve=full(Ve(:,min(J):end));

A=E(Ve(Ve>0),:);
B=zeros(size(Ve));
B(Ve>0)=A(:,1);
C=zeros(size(Ve));
C(Ve>0)=A(:,2);
D=(1:size(Ve,1))'*ones(1,size(Ve,2));
B(B==D)=0;
C(C==D)=0;
Ev=B+C;


 
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
