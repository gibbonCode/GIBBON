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

 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
