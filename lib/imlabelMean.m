function [Q_mean]=imlabelMean(M,ML)

% function [Q_mean]=imlabelMean(M,ML)
% ------------------------------------------------------------------------
%
% This function takes the mean for each of the labeled groups (NaN's
% ignored) in ML according to the intensities in M.
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/01
%------------------------------------------------------------------------

defaultMethod=1; %Hard coded for now
switch defaultMethod    
    case 1 %SPARSE ARRAY BASED (more efficient for large arrays)
        labelSet=unique(ML(:));
        labelSet=labelSet(~isnan(labelSet)); %Remove the nan labelled group
        
        numLabels=max(labelSet);
        
        logic_ML_not_nan=~isnan(ML);
        nnzQ=nnz(logic_ML_not_nan);
        Iq=ML(logic_ML_not_nan);
        Jq=1:nnzQ;
        Sq=M(logic_ML_not_nan);
        sizQ=[numLabels,nnzQ];
        nnzQ=nnz(logic_ML_not_nan);
        Q=sparse(Iq,Jq,Sq,sizQ(1),sizQ(2),nnzQ);
        
        Q_sum=full(sum(Q,2));
        Q_voxelCount=full(sum(spones(Q),2));
        Q_mean=Q_sum./Q_voxelCount;
    case 2 %REGIONPROPS BASED
        ML(isnan(ML))=0;
        A=regionprops(ML,M,'MeanIntensity');
        Q_mean=[A.MeanIntensity];
        Q_mean=Q_mean(:);
end


 
%% 
% ********** _license boilerplate_ **********
% 
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
