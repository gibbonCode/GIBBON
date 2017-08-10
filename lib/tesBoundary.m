function [indBounary]=tesBoundary(F,V)

if numel(V)==1
    numPoints=V;
else
    numPoints=size(V,1);
end

Fbs=sort(F,2);

sizVirt=numPoints*ones(1,size(F,2));

ind_F=subMat2ind(sizVirt,Fbs);

[~,indUni1,~]=unique(Fbs,'rows'); %Get indices for unique faces
F_uni=F(indUni1,:);

ind_F_uni=ind_F(indUni1,:);

ind=1:1:size(F,1);
ind=ind(~ismember(ind,indUni1));
ind_Fb_cut=ind_F(ind,:);
L_uni=~ismember(ind_F_uni,ind_Fb_cut);

indBounary=indUni1(L_uni,:);
 
%% <-- GIBBON footer text --> 
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
