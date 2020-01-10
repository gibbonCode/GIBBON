function [varargout]=quiverCurve(varargin)

%% Parse input

switch nargin
    case 1
        Vg=varargin{1};
        dirOpt=1;
        colorOpt='k';
        alphaLevel=1;
    case 2
        Vg=varargin{1};
        dirOpt=varargin{2};
        colorOpt='k';
        alphaLevel=1;
    case 3
        Vg=varargin{1};
        dirOpt=varargin{2};
        colorOpt=varargin{3};        
        alphaLevel=1;
    case 4
        Vg=varargin{1};
        dirOpt=varargin{2};
        colorOpt=varargin{3};
        alphaLevel=varargin{4};
end

%%

if size(Vg,2)==2
    Vg(:,3)=0;
end

switch dirOpt
    case 1 %Forward
        V=Vg(1:end-1,:);
        U=diff(Vg,1,1);        
    case 2 %Backward
        V=Vg(2:end,:);
        U=flipud(diff(flipud(Vg),1,1));        
    case 3 %Central
        V=Vg;
        Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
        Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
        U=Uf;
        U(:,:,2)=Ub;
        U=gnanmean(U,3);        
end
D=sqrt(sum(U.^2,2));
a=[min(D(:)) max(D(:))];


[F,V,C]=quiver3Dpatch(V(:,1),V(:,2),V(:,3),U(:,1),U(:,2),U(:,3),[],a);

if isempty(colorOpt) %If empty use colormapping
    h=gpatch(F,V,C,'none',alphaLevel);
else %else use specified which could be 'k'
    h=gpatch(F,V,colorOpt,'none',alphaLevel);
end

varargout{1}=h;

 
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
