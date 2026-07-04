function [varargout]=patchCleanUnused(F,V)

% function [Fc,Vc,indFix2,logicValid,indUni]=patchCleanUnused(F,V)
%-------------------------------------------------------------------------
%
%
% 2010 Created
% 2019/03/04 Added basic description header to function
% 2019/04/22 Added handling of cell input 
%-------------------------------------------------------------------------

%%

if isa(F,'cell')
    indUni=[];    
    for q=1:1:numel(F)
        f=F{q};
        indUni=[indUni;f(:)];
    end    
else
    indUni=unique(F(:));
end
indUni=unique(indUni(~isnan(indUni)));

%%

numPointsOriginal=size(V,1); %Number of original points
numPointsNew=numel(indUni); %Number of points in new set
Vc=V(indUni,:); %Select relevant points

%Fix indices in faces matrix
indFix1=1:1:numPointsNew;
indFix2=zeros(numPointsOriginal,1);
indFix2(indUni)=indFix1;

Fc=F;
if isa(F,'cell')
    for q=1:1:numel(F)
        fc=F{q};
        logicValid=~isnan(fc);    
        fc(logicValid)=indFix2(fc(logicValid));
        Fc{q}=fc;
    end
else
    logicValid=~isnan(F);    
    Fc(logicValid)=indFix2(F(logicValid));
end

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
