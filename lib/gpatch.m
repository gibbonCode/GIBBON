function [hp]=gpatch(varargin)

switch nargin
    case 1
        error('Not enough input arguments, provide at least faces and vertices');
    case 2
        F=varargin{1};
        V=varargin{2};
        C='g';
        CE='k';
        A=1;
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE='k';
        A=1;
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE=varargin{4};
        A=1;
    case 5
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE=varargin{4};
        A=varargin{5};
    otherwise
        error('Wrong number of input arguments');
end

if isa(F,'cell') %Assume all entries are cells defining multiple patch data sets
    for q=1:1:numel(F)
        hp(q)=plotPatch(F{q},V{q},C{q},CE{q},A{q});
    end
else
    hp=plotPatch(F,V,C,CE,A);
end

end

%%
function hp=plotPatch(F,V,C,CE,A)

argInPatch.Faces=F;
argInPatch.Vertices=V;
argInPatch.EdgeColor=CE;

if ischar(C) %Plain single color
    argInPatch.FaceColor=C;
elseif size(C,2)==1
    argInPatch.FaceColor='flat';
    %     if size(C,1)>size(F,1) %Assume vertex shading
    %         argInPatch.FaceVertexCData=C;
    %     else %Assume face shading
    argInPatch.CData=double(C);
    %     end
elseif size(C,2)==3  %Assume RGB type
    argInPatch.FaceColor='flat';
    argInPatch.FaceVertexCData=double(C);
else
    error('Invalid color data input');
end

if numel(A)==1 %Plain single alpha
    argInPatch.FaceAlpha=A;
elseif size(A,2)==1 %Alpha mapping
    argInPatch.FaceAlpha='flat';
    argInPatch.FaceVertexAlphaData=A;
else
    error('Invalid alpha data input');
end

hp=patch(argInPatch);

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
