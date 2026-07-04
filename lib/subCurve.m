function [Vs]=subCurve(varargin)

% function [Vs]=subCurve(V,n,closeLoopOpt)
%------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/12/08 Updated to included closed loop option
%------------------------------------------------------------------------

%% Parse Input

switch nargin
    case 1
        V=varargin{1};
        n=1;
        closeLoopOpt=0;
    case 2
        V=varargin{1};
        n=varargin{2};
        closeLoopOpt=0;
    case 3
        V=varargin{1};
        n=varargin{2};
        closeLoopOpt=varargin{3};
    otherwise
        error('Wrong number of input arguments')
end
%%

if closeLoopOpt
    V(end+1,:)=V(1,:); %Add first point to end to close loop
end
    
if n==0 %No subdevision
    Vs=V;     
elseif n>0 %Subdevision of segments
    Vs=zeros(size(V,1)+(size(V,1)-1)*n,size(V,2));
    for q=1:1:size(V,2)
        X=V(:,q);
        XX=linspacen(X(1:end-1),X(2:end),n+2);
        XX=XX(:,1:end-1)';
        Vs(1:end-1,q)=XX(:);
        Vs(end,q)=V(end,q);
    end
elseif n<1 %Throw error
    error('n should be >=0')
end

if closeLoopOpt
    Vs=Vs(1:end-1,:); %Take away last point
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
