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
    
switch nargin
    case 1 %TYPE min(S(:))
        L=S~=0;
        indSparse=find(L);
        [minVal,minInd_L]=min(full(S(L)));
        minInd=indSparse(minInd_L);
    case 2 %L is specified
        indSparse=find(L);
        [minVal,minInd_L]=min(full(S(L)));
        minInd=indSparse(minInd_L);
    case 3 %L and minDim specified
        if isempty(L)
            L=S~=0;
        end
        if isempty(minDim)
            minDim=1;
        end  
        
        %If all zeros are encountered in the direction of minDim then they
        %are set to NaN
        L_allZeros=sum(L,minDim)==0;
        if minDim==1
            S(1,L_allZeros)=NaN;
            L(1,L_allZeros)=1;
        elseif minDim==2
            S(L_allZeros,1)=NaN;
            L(L_allZeros,1)=1;
        end
            
        S(L)=S(L)+1; %ADD 1 TO TREAT ZEROS
        if nargout==1 %Skip indMin
            if minDim==2
                S=S'; %Transpose because sorting provides column permutation orders irrespective of minDim
            end
            Ss=sort(S,1);
            Ls=Ss~=0;
            L_min=cumsum(Ls,1)==1;
            minVal=Ss(L_min);           
        else
            if minDim==2
                S=S'; %Transpose because sorting provides column permutation orders irrespective of minDim
            end
            [Ss,indSort]=sort(S,1); %Warning indSort is full!
            Ls=(Ss~=0);
            L_min=cumsum(Ls,1)==1;
            minVal=full(Ss(L_min));
            minInd=indSort(L_min);            
        end
        minVal=minVal-1; %SUBTRACT 1 AGAIN
end

%


 
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
