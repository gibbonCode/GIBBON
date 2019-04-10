function [Vn]=polyLineIntersect(varargin)

% function [Vn]=polyLineIntersect(V,n_cut,cutLevel,isClosed)

switch nargin
    case 2
        V=varargin{1};
        n_cut=varargin{2};
        cutLevel=0;
        isClosed=0;
    case 3
        V=varargin{1};
        n_cut=varargin{2};
        cutLevel=varargin{3};
        isClosed=0;
    case 4
        V=varargin{1};
        n_cut=varargin{2};
        cutLevel=varargin{3};
        isClosed=varargin{4};
end

%%

if isClosed==1
    E=[(1:size(V,1))' [2:size(V,1) 1]'];
else
    E=[(1:size(V,1)-1)' (2:size(V,1))'];
end

VV=V(E(:,2),:)-V(E(:,1),:);
N=vecnormalize(VV);
P=V(E(:,1),:);

N_cut=n_cut(ones(size(P,1),1),:);

d=cutLevel-dot(P,N_cut,2);

s=dot(N,N_cut,2);
f=d./s;

a=(f.*N);
b=dot(N,a,2);
c=dot(N,VV,2);

f(b>c)=NaN;
f(f<0)=NaN;

Vn=P+(f.*N);

