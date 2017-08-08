function [varargout]=meshBoundary(varargin)

switch nargin
    case 1
        E=varargin{1};
        [F,CF]=element2patch(E,[]); %Get mesh faces
    case 2
        E=varargin{1};
        elementType=varargin{2};
        [F,CF]=element2patch(E,[],elementType); %Get mesh faces
    case 3        
        E=varargin{1};        
        elementType=varargin{2};        
        C=varargin{3};
        [F,CF]=element2patch(E,C,elementType); %Get mesh faces
    otherwise
        error('Wrong number of input arguments');
end

if isempty(CF)
    CF=zeros(size(F,1),1);
end

%Check shared faces faces
numFacesIni=size(F,1);
[F_uni,indF,IND_F_2]=uniqueIntegerRow(F);

CF=CF(indF,:);
numFacesUni=size(F_uni,1);

%Get face counts
logicColourMatrixEntry=sparse(IND_F_2,1:numFacesIni,1,numFacesUni,numFacesIni,numFacesIni);
F_count=full(sum(logicColourMatrixEntry,2));

%Compose boundary set from faces that are used once
logicBoundary=F_count==1; 
Fb=F_uni(logicBoundary,:);
CFb=CF(logicBoundary,:);

varargout{1}=Fb;
varargout{2}=F_count;
varargout{3}=CFb;
 
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
