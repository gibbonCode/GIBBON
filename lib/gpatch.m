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
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
