function [minVal,minInd]=sparseMin(S,L,minDim)

% clear all; close all; clc;
% S=rand(6,7);
% indRand=randperm(numel(S));
% indRand=indRand(1:round(numel(S)/2));
% S(indRand)=0;
% minDim=2;
% 
% S=sparse(S);
% L=S~=0;
%

%%

warning('This function is depricated. Use spmin instead');
[minVal,minInd]=spmin(S,[],minDim,'includenan',L);
minVal=full(minVal(:));
minInd=full(minInd(:));

%% OLD 
% switch nargin
%     case 1 %TYPE min(S(:))
%         L=S~=0;
%         indSparse=find(L);
%         [minVal,minInd_L]=min(full(S(L)));
%         minInd=indSparse(minInd_L);
%     case 2 %L is specified
%         indSparse=find(L);
%         [minVal,minInd_L]=min(full(S(L)));
%         minInd=indSparse(minInd_L);
%     case 3 %L and minDim specified
%         if isempty(L)
%             L=S~=0;
%         end
%         if isempty(minDim)
%             minDim=1;
%         end  
%         
%         %If all zeros are encountered in the direction of minDim then they
%         %are set to NaN
%         L_allZeros=sum(L,minDim)==0;
%         if minDim==1
%             S(1,L_allZeros)=NaN;
%             L(1,L_allZeros)=1;
%         elseif minDim==2
%             S(L_allZeros,1)=NaN;
%             L(L_allZeros,1)=1;
%         end
%             
%         S(L)=S(L)+1; %ADD 1 TO TREAT ZEROS
%         if nargout==1 %Skip indMin
%             if minDim==2
%                 S=S'; %Transpose because sorting provides column permutation orders irrespective of minDim
%             end
%             Ss=sort(S,1);
%             Ls=Ss~=0;
%             L_min=cumsum(Ls,1)==1;
%             minVal=Ss(L_min);           
%         else
%             if minDim==2
%                 S=S'; %Transpose because sorting provides column permutation orders irrespective of minDim
%             end
%             [Ss,indSort]=sort(S,1); %Warning indSort is full!
%             Ls=(Ss~=0);
%             L_min=cumsum(Ls,1)==1;
%             minVal=full(Ss(L_min));
%             minInd=indSort(L_min);            
%         end
%         minVal=minVal-1; %SUBTRACT 1 AGAIN
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
