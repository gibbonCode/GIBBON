function [varargout]=patchCleanUnused(F,V)

% function [Fc,Vc,indFix2,logicValid,indUni]=patchCleanUnused(F,V)
%-------------------------------------------------------------------------
%
%
% 2010 Created
% 2019/03/04 Added basic description header to function
%-------------------------------------------------------------------------

%%

logicValid =F>0;% %Treat 0,NaN,inf

numPoints=size(V,1);
indUni=unique(F(logicValid)); %Unique indices of used vertices
Vc=V(indUni,:); %Select relevant points

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(numPoints,1);
indFix2(indUni)=indFix1;
Fc=F; 
Fc(logicValid)=indFix2(F(logicValid));

%Output
varargout{1}=Fc;
varargout{2}=Vc;
varargout{3}=indFix2;
varargout{4}=logicValid;
varargout{5}=indUni;
 
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
