function [C]=polyContourThick(varargin)

switch nargin
    case 1
        V=varargin{1};
        k=[];
        v=[];
        pointSpacing=[];
        resampleMethod='linear';
    case 2
        V=varargin{1};
        k=varargin{2};
        v=[];
        pointSpacing=[];
        resampleMethod='linear';
    case 3
        V=varargin{1};
        k=varargin{2};
        v=varargin{3};
        pointSpacing=[];
        resampleMethod='linear';
    case 4
        V=varargin{1};
        k=varargin{2};
        v=varargin{3};
        pointSpacing=varargin{4};
        resampleMethod='linear';
    case 5        
        V=varargin{1};
        k=varargin{2};
        v=varargin{3};
        pointSpacing=varargin{4};
        resampleMethod=varargin{5};
end

%%

if size(V,2)==3
   warning('polyContourThick is for 2D polygons. 3D detected, ignoring 3rd dimension');  
   V=V(:,[1 2]);
end

D=sqrt(sum(diff(V,1,1).^2,2));
if isempty(k)
    k=mean(D)*2;
end

if isempty(v)
    v=mean(D)/3;
end

if isempty(pointSpacing)
    pointSpacing=mean(D);
end

%%
%Create coordinate matrices
minV=min(V,[],1)-2*k;
maxV=max(V,[],1)+2*k;
range_V=maxV-minV;
n=round(range_V./v);
[X,Y]=meshgrid(linspace(minV(1),maxV(1),n(1)),linspace(minV(2),maxV(2),n(2)));
Vg=[X(:) Y(:)];

%Derive distance based level-set image
D=minDist(Vg,V); %Nearest point distance between image grid points and polygon
M=reshape(D,size(X)); %Level-set image

%Compute contour
[C]=gcontour(X,Y,M,k,pointSpacing,resampleMethod);

