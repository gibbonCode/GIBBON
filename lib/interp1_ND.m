function [Yi]=interp1_ND(X,Y,Xi,interpDim,interpMethod)

% function [Yi]=interp1_ND(X,Y,Xi,interpDim,interpMethod)
% ------------------------------------------------------------------------
% This function is similar to interp1. However it can perform 1D
% interpolation for multidimensional arrays. E.g. Time interpolation for 3D
% image data varying in time. The direction of interpolation is specified
% by interpDim, and the method by interpMethod
% ('nearest','linear','cubic','pchip').
%
%
% Reference for pchip method:
%   F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
%   Interpolation", SIAM J. Numerical Analysis 17, 1980, 238-246.
%   David Kahaner, Cleve Moler and Stephen Nash, Numerical Methods
%   and Software, Prentice Hall, 1988.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/03/04
%------------------------------------------------------------------------

%% Check input

%Check number of dimensions
nDimsX = ndims(X);
nDimsXi = ndims(Xi);

if nDimsX~=nDimsXi
    if size(Xi,nDimsX)~=1
        error('Input arrays have different number of dimensions');
    end
end

%Check column/row form
if iscolumn(X)
    X=X';
end

if iscolumn(Xi)
    Xi=Xi';
end

if isvector(X)
    interpDim=2;
end

%% Cast in 2D array format

siz=size(X);
siz_i=size(Xi);

dimInd=1:nDimsX;
logicNotInterpDim=dimInd~=interpDim;
permOrder=[sort(dimInd(logicNotInterpDim)) interpDim];

%Permuting array order so that interpolation dimension is last
X = permute(X,permOrder);
Y = permute(Y,permOrder);
Xi = permute(Xi,permOrder);

%Convert to 2D array with interpolation performed allong column direction
siz1=prod(siz(logicNotInterpDim));
X=reshape(X(:),[siz1, numel(X)/siz1]);
Y=reshape(Y(:),[siz1, numel(Y)/siz1]);
Xi=reshape(Xi(:),[siz1, numel(Xi)/siz1]);

%%

siz_2D=size(X);
siz_i_2D=size(Xi);

switch interpMethod
    case {'cubic','pchip'}
        diff_Y=diff(Y,1,2);
        diff_X=diff(X,1,2);        
        K=diff_Y./diff_X;
end

switch interpMethod
    case 'cubic'
        M=zeros(siz_2D);
        M(:,1)=K(:,1);
        M(:,end)=K(:,end);
        M(:,2:end-1)=0.5*(K(:,1:end-1)+K(:,2:end));
    case 'pchip'
        M = pchipSlopes(X,Y,K);        
end

%%

I=(1:siz_2D(1))';

%Initialize
Yi=nan(siz_i_2D);
p0=zeros(size(X,1),1); p1=p0;
switch interpMethod
    case 'cubic'
        m0=p0; m1=p0;
    case 'pchip'
        m0=p0; m1=p0;
end

for q=1:1:siz_i_2D(end)
    
    %Find close point pair to define interval
    D=abs(X-Xi(:,q*ones(1,size(X,2))));
    
    [~,J1]=min(D,[],2);
    indMin1=sub2ind(siz_2D,I,J1);
    D(indMin1)=inf;
    [~,J2]=min(D,[],2);
    indSort=[J1 J2];
    indSort=sort(indSort(:,1:2),2); %Sorted indices of closest two
    
    ind1=sub2ind(siz_2D,I,indSort(:,1));
    ind2=sub2ind(siz_2D,I,indSort(:,2));
    
    w=X(ind2)-X(ind1); %interval width
    t=(Xi(:,q)-X(ind1))./w; %t coordinate in domain
    p0(1:end)=Y(ind1); %First point
    p1(1:end)=Y(ind2); %End point
    switch interpMethod
        case {'cubic','pchip'}
            m0(1:end)=M(ind1); %First derivative
            m1(1:end)=M(ind2); %Final derivative            
            
            h00=2.*t.^3-3.*t.^2+1;
            h10=t.^3-2*t.^2+t;
            h01=-2*t.^3+3*t.^2;
            h11=t.^3-t.^2;
            pt=h00.*p0 + w.*h10.*m0 + h01.*p1 + w.*h11.*m1;                        
        case 'linear'
            pt=t.*p1+(1-t).*p0;
        case 'nearest'
            logicPrev=t<=0.5;
            pt(logicPrev)=p0(logicPrev);
            pt(~logicPrev)=p1(~logicPrev);
    end
    
    Yi(:,q)=pt; %Store in array
    
end

%% Cast back in original form

% Xi=reshape(Xi(:),siz_i(permOrder));
% Xi=ipermute(Xi,permOrder);

Yi=reshape(Yi(:),siz_i(permOrder));
Yi=ipermute(Yi,permOrder);

end
%%
function M = pchipSlopes(X,Y,dYdX)
%PCHIPSLOPES  Derivative values for shape-preserving Piecewise Cubic Hermite
% Interpolation.
% d = pchipslopes(x,y,del) computes the first derivatives, d(k) = P'(x(k)).

n = size(X,2);

%  Slopes at interior points.
%  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
%  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

M = zeros(size(Y));

[I,J]= find(sign(dYdX(:,1:n-2)).*sign(dYdX(:,2:n-1)) > 0);
ind=sub2ind(size(X),I,J);
ind1=sub2ind(size(X),I,J+1);

h = diff(X,1,2);
hs = h(ind)+h(ind1);
w1 = (h(ind)+hs)./(3*hs);
w2 = (hs+h(ind1))./(3*hs);
dmax = max(abs(dYdX(ind)), abs(dYdX(ind1)));
dmin = min(abs(dYdX(ind)), abs(dYdX(ind1)));
M(ind1) = dmin./conj(w1.*(dYdX(ind)./dmax) + w2.*(dYdX(ind1)./dmax));

%  Slopes at end points.
%  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.

M(:,1) = ((2*h(:,1)+h(:,2)).*dYdX(:,1) - h(:,1).*dYdX(:,2))./(h(:,1)+h(:,2));
L = (sign(M(:,1)) ~= sign(dYdX(:,1)));
M(L,1) = 0;

L = ((sign(dYdX(:,1)) ~= sign(dYdX(:,2))) & (abs(M(:,1)) > abs(3*dYdX(:,1))));
M(L,1) = 3*dYdX(L,1);

M(:,n) = ((2*h(:,n-1)+h(:,n-2)).*dYdX(:,n-1) - h(:,n-1).*dYdX(:,n-2))./(h(:,n-1)+h(:,n-2));

L= sign(M(:,n)) ~= sign(dYdX(:,n-1));
M(L,n) = 0;
    
L= (sign(dYdX(:,n-1)) ~= sign(dYdX(:,n-2))) & (abs(M(:,n)) > abs(3*dYdX(:,n-1)));
M(L,n) = 3*dYdX(L,n-1);

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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
