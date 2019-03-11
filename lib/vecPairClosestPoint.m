function [varargout]=vecPairClosestPoint(varargin)

%%

switch nargin
    case 2
        a=varargin{1};
        b=varargin{2};
        A=zeros(size(a));
        B=zeros(size(a));
        forceOpt=1;
    case 3
        a=varargin{1};
        b=varargin{2};
        A=zeros(size(a));
        B=zeros(size(a));
        forceOpt=varargin{3};
    case 4
        a=varargin{1};
        b=varargin{2};
        A=varargin{3};
        B=varargin{4};
        forceOpt=1;
    case 5
        a=varargin{1};
        b=varargin{2};
        A=varargin{3};
        B=varargin{4};
        forceOpt=varargin{5};
end

%% Compute closest points

c=B-A;

tt=(-dot(a,b,2).*dot(b,c,2))+(dot(a,c,2).*dot(b,b,2));
bb=( dot(a,a,2).*dot(b,b,2))-(dot(a,b,2).*dot(a,b,2));
dd=(tt./bb);
AX=A+a.*dd(:,ones(size(A,2),1));

tt=( dot(a,b,2).*dot(a,c,2))-(dot(b,c,2).*dot(a,a,2));
dd=(tt./bb);
BX=B+b.*dd(:,ones(size(A,2),1));

%% Compute logics reflecting if points are on line segment if needed

if forceOpt==1 || nargout>2
    AXa=dot(a,AX-A,2);
    aa=dot(a,a,2);
    logicBefore_A=AXa<0;
    logicAfer_Aa=AXa>aa;
    logicOnSegment_A=~logicBefore_A & ~logicAfer_Aa;
    
    BXb=dot(b,BX-B,2);
    bb=dot(b,b,2);
    logicBefore_B=BXb<0;
    logicAfer_Bb=BXb>bb;
    logicOnSegment_B=~logicBefore_B & ~logicAfer_Bb;
end

%% Force output points onto line segment if needed

if forceOpt==1
    Aa=A+a;
    AX(logicBefore_A,:)=A(ones(nnz(logicBefore_A),1),:);
    AX(logicAfer_Aa,:)=Aa(ones(nnz(logicAfer_Aa),1),:);
    
    Bb=B+b;
    BX(logicBefore_B,:)=B(ones(nnz(logicBefore_B),1),:);
    BX(logicAfer_Bb,:)=Bb(ones(nnz(logicAfer_Bb),1),:);
end

%% Collect output 

varargout{1}=AX;
varargout{2}=BX;
if nargout>2
    varargout{3}=logicOnSegment_A;
    varargout{4}=logicOnSegment_B;
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
