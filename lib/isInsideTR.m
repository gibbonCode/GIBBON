function [varargout]=isInsideTR(varargin)

if nargin<2
    error('Not enough input arguments');
end

TR=varargin{1};
numTR=size(TR.ConnectivityList,1);
QP=varargin{2};

if nargin==3
    ti=varargin{3};
else
    ti=1:numTR;
end
    
%Get barycentric coordinates of points
baryCoords=cartesianToBarycentric(TR,ti(:),QP);

logicFoundEnclosing=all(baryCoords>0,2);

varargout{1}=logicFoundEnclosing;
if nargout==2
    varargout{2}=baryCoords;
end



 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
