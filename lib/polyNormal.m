function [varargout]=polyNormal(varargin)

% [N,VN]=polyNormal(V_poly)
% ------------------------------------------------------------------------
% Normals directions are derived based on cross product of curve segments
% with the outware z-direction. A horizontal and clockwise curve segment
% therefore results in an upward normal. The normals are provided for each
% segment or for each point on the curve. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2017/03/16 Created
%------------------------------------------------------------------------

%% Parse input

defaultOptionStruct.closeLoopOpt=false;
defaultOptionStruct.type='vertex'; % 'vertex'/'node' or 'segment'/'edge' 
defaultOptionStruct.zDir=[0 0 1];

switch nargin
    case 1
        V=varargin{1};
        optionStruct=defaultOptionStruct;
    case 2
        V=varargin{1};
        optionStruct=varargin{2};
end

[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

zDir=vecnormalize(optionStruct.zDir);
if size(zDir,1)<size(V,1)
    zDir=zDir(ones(size(V,1),1),:); % Extend to match V
end
closeLoopOpt=optionStruct.closeLoopOpt;
typeOpt=optionStruct.type;

%%
nd=size(V,2);
if nd==2 %Cope with 2D
    V(:,3)=0;
end

%%

switch typeOpt
    case {'segment','edge'}
        if closeLoopOpt
            V(end+1,:)=V(1,:);
            zDir(end+1,:)=zDir(1,:);
        end
        U  = V(2:end,:)-V(1:end-1,:);
        ZN = (zDir(1:end-1,:)+zDir(2:end,:))/2;
        N  = vecnormalize(cross(ZN,U));
        VN = (V(1:end-1,:)+V(2:end,:))/2;
    case {'node','vertex'}
        if closeLoopOpt                   
            U = [V(2:end,:)-V(1:end-1,:); V(1,:)-V(end,:)];
            U = [U(1,:)+U(end,:); (U(1:end-1,:)+U(2:end,:))]./2; %Central difference         
        else
            U = V(2:end,:)-V(1:end-1,:);
            U = [U(1,:); (U(1:end-1,:)+U(2:end,:))/2; U(end,:)]; %Central difference         
        end
        N = cross(zDir,U);
        N = vecnormalize(N);
        VN=V;
end

if nd==2 %Remove 3rd dimension from vectors and coordinates for 2D input
    N=N(:,[1 2]);
    VN=VN(:,[1 2]);
end

%% Collect output
varargout{1}=N; 
if nargout==2
    varargout{2}=VN;
end

end

% function [V,U]=polyDerivative(Vg,dirOpt)
% 
% switch dirOpt
%     case 1 %Forward
%         V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
%         U=diff(Vg,1,1);        
%     case 2 %Backward
%         V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
%         U=flipud(diff(flipud(Vg),1,1));        
%     case 3 %Central
%         V=Vg;
%         Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
%         Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
%         U=Uf;
%         U(:,:,2)=Ub;
%         U=gnanmean(U,3);        
% end
% U=vecnormalize(U);
% 
% end
 

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
