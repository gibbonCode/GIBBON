function [varargout]=grid2patch(varargin)

% function [F,V,C,CV]=grid2patch(X,Y,Z,C,perdiocOpt);
% ------------------------------------------------------------------------
%
%
% 2019/12/16: Created as alternative to surf2patch (missing in Octave)
% ------------------------------------------------------------------------

%% parse input

switch nargin 
    case 2
        X=varargin{1};
        Y=varargin{2};
        Z=zeros(size(X));
        C=[];
        perdiocOpt=false(1,2);
    case 3
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        C=[];
        perdiocOpt=false(1,2);
    case 4
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        C=varargin{4};        
        perdiocOpt=false(1,2);
    case 5
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        C=varargin{4};
        perdiocOpt=varargin{5};
end

if isempty(C)
    C=Z(:);
end
siz=size(X);

%%

V=[X(:) Y(:) Z(:)]; %Vertex set

%Create row of faces
f=[1 1+siz(1) 2+siz(1) 2]; %First element 
q=(0:1:siz(1)-2)';
qq=q(:,ones(4,1));
ff=f(ones(siz(1)-1,1),:)+qq; %Row of faces by copying first

if perdiocOpt(1)==1
    ff(end+1,:)=[siz(1) siz(1)+siz(1) 1+siz(1) 1];
end

%Create grid of faces by copying column
q=0:siz(1):(siz(1)*(siz(2)-2));
if perdiocOpt(1)==1
    Q=q(ones(siz(1),1),:);
else
    Q=q(ones(siz(1)-1,1),:);
end
Q=Q(:);
Q=Q(:,ones(1,4));
F=repmat(ff,[siz(2)-1 1])+Q; %Grid of faces

if perdiocOpt(1)==1 && perdiocOpt(2)==1
    t=siz(1)*(siz(2)-1);
    f=[2 2+t 1+t 1]; %First element
    q=(0:1:siz(1)-2)';
    qq=q(:,ones(4,1));
    ff=f(ones(siz(1)-1,1),:)+qq; %Row of faces by copying first
    ff(end+1,:)=[1 1+t siz(1)+t siz(1)];
    F=[F;ff];
elseif perdiocOpt(2)==1
    t=(siz(1)*(siz(2)-1));
    f=[2 2+t 1+t 1]; %First element
    q=(0:1:siz(1)-2)';
    qq=q(:,ones(4,1));
    ff=f(ones(siz(1)-1,1),:)+qq; %Row of faces by copying first
    F=[F;ff];
end

%% Collect output

varargout{1}=F; 
varargout{2}=V; 
if nargout==3
    varargout{3}=vertexToFaceMeasure(F,C(:));
end
if nargout==4
    A=C;
    varargout{4}=A(:);
end

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
