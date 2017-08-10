function [V,F]=voronoiDiagramEven2D(varargin)

%% DERIVE STANDARD 2D VORONOI DIAGRAM

%Get Delaunay tesselation
DT=varargin{1};

%Computer standard Voronoi output containing vertices matrix and faces cell
[Vv,Fv_cell] =voronoiDiagram(DT);
numCell=numel(Fv_cell);

%Replace inf entries by NaN
Vv(isinf(Vv))=NaN;

%% Creating even Voronoi cell sampling

%Get number of points 
switch nargin    
    case 2
        np=varargin{2};
        interpMethod='linear';
    case 3
        np=varargin{2};
        interpMethod=varargin{3};
    otherwise
        np=[];
        interpMethod='linear';
end

%Use maximum occuring number of points if np=[]
if isempty(np)
    np=max(cellfun(@numel,Fv_cell));
end

%Resampling Voronoi cells
V=zeros(numCell*np,2); %Vertex matrix
F=reshape([1:size(V,1)]',np,size(V,1)/np)'; %Face matrix
for q=1:1:numCell
    
    %Get current Voronoi cell vertices   
    vv=Vv(Fv_cell{q},:);
    
    %Resample cell all vertices are valid
    if any(isnan(vv))
        vv=NaN;
    else        
        [vv]=evenlySampleCurve(vv,np,interpMethod,1);
    end
    
    %Store new points in array
    V(F(q,:),:)=vv;        
    
end

end
 
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
