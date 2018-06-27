function [varargout]=tri3_tri6(varargin)

% function [TRI6,V6,VX6C]=tri3_tri6(TRI3,V3,VXC)
% 
% This function converts 3 node (e.g. linear) triangular elements into 6
% node (e.g. quadratic) triangular elements compatible with FEBio. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 09/12/2013
% 2018/06/26
%----------------------------------------------------------------------
%%

switch nargin
    case 2
        TRI3=varargin{1};
        V3=varargin{2};
        VXC={};
    case 3
        TRI3=varargin{1};
        V3=varargin{2};
        VXC=varargin{3};
end

%%
%Collect nodes
V_1_3=V3;
V_4 =0.5.*(V3(TRI3(:,1),:)+V3(TRI3(:,2),:));
V_5 =0.5.*(V3(TRI3(:,2),:)+V3(TRI3(:,3),:));
V_6 =0.5.*(V3(TRI3(:,3),:)+V3(TRI3(:,1),:));
V6=[V_1_3;V_4;V_5;V_6;];


%Define elements
numTris=size(TRI3,1);
num_V_1_3=size(V_1_3,1);
indTris=(1:numTris)';
TRI6=[TRI3...                            % 1-3
       indTris+num_V_1_3+(numTris*(1-1))... % 4
       indTris+num_V_1_3+(numTris*(2-1))... % 5
       indTris+num_V_1_3+(numTris*(3-1))... % 6       
       ];  
   
%Removing double points
[TRI6,V6,ind1,ind2]=mergeVertices(TRI6,V6); 

%Derive VX6C
if ~isempty(VXC)    
    VX6C=VXC;
    for q=1:1:numel(VXC)
        VX=VXC{q};
        VX_1_3=VX;
        VX_4 =0.5.*(VX(TRI3(:,1),:)+VX(TRI3(:,2),:));
        VX_5 =0.5.*(VX(TRI3(:,2),:)+VX(TRI3(:,3),:));
        VX_6 =0.5.*(VX(TRI3(:,3),:)+VX(TRI3(:,1),:));
        VX6=[VX_1_3;VX_4;VX_5;VX_6;];
        VX6C{q}=VX6(ind1,:);
    end
else
    VX6C={};
end

varargout{1}=TRI6;
varargout{2}=V6;
varargout{3}=VX6C;

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
