function [varargout]=gpatch(varargin)

% function [h]=gpatch(F,V,C,CE,A,L)
% ------------------------------------------------------------------------
% This function is a short-hand version of the |patch| command. The inputs
% for |gpatch| are the faces (F), the vertices (V), the color description
% (C), the edge color description CE, the transparancy (A), and the edge
% width (L). 
% The color data descriptions C (or equivalently CE for edges) can be: 
% 1) A string such as 'g' for green
% 2) A triplet of RGD values e.g. [1 0 0] is blue
% 3) A nx1 or a mx1 array of colormapped colors (where n=size(F,1) or
% m=size(V,1)) 
% 4) (simiarl to 3) A nx3 or a mx3 RGB color value array for the faces or
% vertices respectively. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2017 
% 2018/02/07 Added support for colormapped edges
% 2019/07/03 Added handling of empty alpha data
% 2021/05/12 Added "material dull" style lighting settings by default to
% avoid white color reflectance distorting colormapped visualizations. 
%------------------------------------------------------------------------

switch nargin
    case 1
        error('Not enough input arguments, provide at least faces and vertices');
    case 2
        F=varargin{1};
        V=varargin{2};
        C='gw';
        CE='k';
        faceAlpha=1;
        lineWidth=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE='k';
        faceAlpha=1;
        lineWidth=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE=varargin{4};
        faceAlpha=1;
        lineWidth=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE=varargin{4};
        faceAlpha=varargin{5};
        lineWidth=[];
    case 6
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        CE=varargin{4};
        faceAlpha=varargin{5};
        lineWidth=varargin{6};
    otherwise
        error('Wrong number of input arguments');
end

if isempty(C)
    C='g';
end

if isempty(CE)
    C='k';
end

if isempty(faceAlpha)
    faceAlpha=1;
end

if isa(F,'cell') %Assume all entries are cells defining multiple patch data sets
    hp=gobjects(numel(F),1);
    for q=1:1:numel(F)
        f=F{q};
        
        if isa(V,'cell')
            v=V{q};
        else
            v=V;
        end

        if isa(C,'cell')
            c=C{q};
        else
            c=C;
        end
        
        if isa(CE,'cell')
            ce=CE{q};
        else
            ce=CE;
        end
        
        if isa(faceAlpha,'cell')
            a=faceAlpha{q};
        else
            a=faceAlpha;
        end
        hp(q)=plotPatch(f,v,c,ce,a,lineWidth);
    end
else
    hp=plotPatch(F,V,C,CE,faceAlpha,lineWidth);
end

if nargout==1
    varargout{1}=hp;
end

end

%%
function hp=plotPatch(F,V,C,CE,faceAlpha,lineWidth)

argInPatch.Faces=F;
argInPatch.Vertices=V;
argInPatch.EdgeColor=CE;

%Avoid reflectance (similar to material dull)
argInPatch.DiffuseStrength=0.6;
argInPatch.AmbientStrength=0.4;
argInPatch.SpecularExponent=10;
argInPatch.SpecularStrength=0;
argInPatch.SpecularColorReflectance=0;

if ischar(C) %Plain single color
    argInPatch.FaceColor=C;
    if strcmp(C,'kw')
        argInPatch.FaceColor=grayColor(0.5);
    end
    if strcmp(C,'rw')
        argInPatch.FaceColor=[1 0.5 0.5];
    end
    if strcmp(C,'gw')
        argInPatch.FaceColor=[0.5 1 0.5];
    end
    if strcmp(C,'bw')
        argInPatch.FaceColor=[0.5 0.5 1];
    end
    if strcmp(C,'yw')
        argInPatch.FaceColor=[1 1 0.5];
    end
    if strcmp(C,'cw')
        argInPatch.FaceColor=[0.5 1 1];
    end
    if strcmp(C,'mw')
        argInPatch.FaceColor=[1 0.5 1];
    end
    if strcmp(C,'o') %Orange
        argInPatch.FaceColor=orange;
    end
elseif size(C,2)==1
    argInPatch.FaceColor='flat';
    argInPatch.CData=double(C);
elseif size(C,2)==3 && size(C,1)==1 %Assume single RGB level
    argInPatch.FaceColor=double(C);
elseif size(C,2)==3 && size(C,1)>1 %Assume RGB array
    argInPatch.FaceColor='flat';
    argInPatch.FaceVertexCData=double(C);
else
    error('Invalid face-vertex color data input');
end

if ischar(CE) %Plain single color
    argInPatch.EdgeColor=CE;    
    if strcmp(CE,'kw')
        argInPatch.EdgeColor=grayColor(0.5);
    end
    if strcmp(CE,'rw')
        argInPatch.EdgeColor=[1 0.5 0.5];
    end
    if strcmp(CE,'gw')
        argInPatch.EdgeColor=[0.5 1 0.5];
    end
    if strcmp(CE,'bw')
        argInPatch.EdgeColor=[0.5 0.5 1];
    end
    if strcmp(CE,'yw')
        argInPatch.EdgeColor=[1 1 0.5];
    end
    if strcmp(CE,'cw')
        argInPatch.EdgeColor=[0.5 1 1];
    end
    if strcmp(CE,'mw')
        argInPatch.EdgeColor=[1 0.5 1];
    end
    if strcmp(CE,'o') %Orange
        argInPatch.EdgeColor=orange;
    end
elseif size(CE,2)==1
    if size(CE,1)>1
        if size(CE,1)==size(F,1)
            [CE]=faceToVertexMeasure(F,V,CE);
            argInPatch.EdgeColor='flat';
            argInPatch.CData=double(CE);
        end
        if size(CE,1)==size(V,1)            
            argInPatch.EdgeColor='flat';
            argInPatch.CData=double(CE);
        end
    else
        argInPatch.EdgeColor='flat';
        argInPatch.CData=double(CE)*ones(size(V,1),1);
    end    
elseif size(CE,2)==3 && size(CE,1)==1 %Assume single RGB level
    argInPatch.EdgeColor=double(CE);
elseif size(CE,2)==3 && size(CE,1)>1 %Assume RGB array
    argInPatch.EdgeColor='flat';
    argInPatch.FaceVertexCData=double(CE);
else
    error('Invalid edge color data input');
end

if numel(faceAlpha)==1 %Plain single alpha
    argInPatch.FaceAlpha=double(faceAlpha);
elseif size(faceAlpha,2)==1 %Alpha mapping
    argInPatch.FaceAlpha='flat';
    argInPatch.FaceVertexAlphaData=double(faceAlpha);
else
    error('Invalid alpha data input');
end

if ~isempty(lineWidth)
    argInPatch.LineWidth=lineWidth;
end

hp=patch(argInPatch);

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
