function [linearIndexVoigt,linearIndexFourthOrder]=tensor2voigtMap(c)

siz2=3*ones(1,2);
siz4=3*ones(1,4);
if ndims(c)==4 %4d array e.g. 4th order tensor
    if all(size(c)==siz4) %4th order tensor
        %Creating indices for tensors
        sizeFourthOrderTensor=siz4; %size of fourth order tensor
        sizeVoigtTensor=6.*ones(1,2); %size
        indexCombinations=gcombvec(1:3,1:3,1:3,1:3)';
        I=indexCombinations(:,1); J=indexCombinations(:,2); K=indexCombinations(:,3); L=indexCombinations(:,4);
        
        kroneckerDelta_IJ=double(I==J);
        kroneckerDelta_KL=double(K==L);
        
        iVoigt=(I.*kroneckerDelta_IJ)+((1-kroneckerDelta_IJ).*(9-I-J)); %Voigt mapping
        jVoigt=(K.*kroneckerDelta_KL)+((1-kroneckerDelta_KL).*(9-K-L)); %Voigt mapping
        
        linearIndexFourthOrder=sub2ind(sizeFourthOrderTensor,I,J,K,L); %Linear version of subscript indices
        linearIndexVoigt=sub2ind(sizeVoigtTensor,iVoigt,jVoigt); %Linear version of subscript indices

    else
        %TO DO add warning here
    end
elseif ismatrix(c) %2D array
    if all(size(c)==siz2) %2nd order tensor
        linearIndexVoigt=[1:6 1:6];        
        linearIndexFourthOrder=sub2ind(siz2,[1 2 3 2 1 1 1 2 3 3 3 2],[1 2 3 3 3 2 1 2 3 2 1 1]);
    else
         %TO DO add warning here
    end
else
    %TO DO add warning here
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
