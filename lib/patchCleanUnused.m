function [varargout]=patchCleanUnused(F,V)

logicValid =F>0;% %Treat 0,NaN,inf

numPoints=size(V,1);
indUni=unique(F(logicValid)); %Unique indices of used vertices
Vc=V(indUni,:); %Select relevant points

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(numPoints,1);
indFix2(indUni)=indFix1;
Fc=F; 
Fc(logicValid)=indFix2(F(logicValid));

%Output
varargout{1}=Fc;
varargout{2}=Vc;
varargout{3}=indFix2;
 
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
