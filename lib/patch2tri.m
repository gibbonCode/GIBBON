function [varargout]=patch2tri(varargin)

% function [Ft,Vt,Ct]=patch2tri(F,V,C)
% ------------------------------------------------------------------------
% This function converts patch data to triangular patch data. It is assumed
% that the centroid (mean of patch vertices) lies on the patch and
% connection of the patch edges to the centroid yields appropriate
% triangles. 
% 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/12/22 Created based on quad2tri
%------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};        
        C=[];
    case 3        
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};    
    otherwise
        error('False number of input arguments');
end

%%

V_new=zeros(size(F,1),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    V_new(:,q)=mean(X(F),2);    
end

columnCopy=repmat((1:size(F,2))',1,2)'; columnCopy=[columnCopy(2:end) 1];
Ft=F(:,columnCopy);
Ft=reshape(Ft',2,numel(Ft)/2)';
indAdd=size(V,1)+(1:size(V_new,1));
Vt=[V;V_new];
indAddRep=repmat(indAdd(:)',size(F,2),1);
Ft(:,3)=indAddRep(:);

%Color vector now assumes all schemes allow for simple copying of the color
%information
if ~isempty(C)
    Ct=C(Ft(:,3)-size(V,1));
else
    Ct=[];
end

%% Collect output
varargout{1}=Ft;
varargout{2}=Vt;
varargout{3}=Ct;
 
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
