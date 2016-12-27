function [h]=gpatch(varargin)

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
end

%%

argInPatch.Faces=F;
argInPatch.Vertices=V;
argInPatch.EdgeColor=CE;

if ischar(C) %Plain single color
    argInPatch.FaceColor=C; 
elseif size(C,2)==1        
    argInPatch.FaceColor='flat'; 
    argInPatch.CData=C;
elseif size(C,2)==3
    argInPatch.FaceColor='flat';
    argInPatch.FaceVertexCData=C;
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

h=patch(argInPatch);

