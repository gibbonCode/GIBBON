function [F_uni,V_uni,C_uni,IND_V,IND_F,F_count]=unique_patch(F,V,C,numDigitKeep)

numFacesIni=size(F,1);

%Removing unused vertices
[F,V]=removeNotIndexed(F,V);

%Removing double vertices
try
    [~,IND_V,IND_IND]=unique(round(V,numDigitKeep,'significant'),'rows');
catch
    [~,IND_V,IND_IND]=unique(sround(V,numDigitKeep),'rows');
end
V_uni=V(IND_V,:);
F=IND_IND(F); %Fix indices in F

%Removing double FACES
[F_uni,IND_F,IND_F_2]=uniqueIntegerRow(F);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Fixing face colors, shared faces now obtain mean colour
if ~isempty(C)
    sharedColourMatrixSparse=sparse(IND_F_2,1:numFacesIni,C,numFacesUni,numFacesIni,numFacesIni);    
    C_uni=full(sum(sharedColourMatrixSparse,2))./F_count;
else
    C_uni=[];
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
