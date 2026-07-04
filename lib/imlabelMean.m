function [varargout]=imlabelMean(varargin)

% function [Q_mean]=imlabelMean(M,ML,methodOpt)
% ------------------------------------------------------------------------
%
% This function takes the mean for each of the labeled groups (NaN's
% ignored) in ML according to the intensities in M.
%
%
% Kevin Mattheus Moerman
% 2014/04/01
%------------------------------------------------------------------------

%%

switch nargin
    case 2
        M=varargin{1};
        ML=varargin{2};
        methodOpt=1; 
    case 3
        M=varargin{1};
        ML=varargin{2};
        methodOpt=varargin{3};
end

%%
logic_ML_not_nan=~isnan(ML); %Logic for nan entries
[labelSet,~,ind2]=unique(ML(logic_ML_not_nan)); %Get the unique label set and ind2 (reconstructor)
numLabels=numel(labelSet); %The number of unique labels
labelIndices=1:1:numLabels; %Label indices (used to replace arbitary labels e.g. -1,0,pi to 1,2,3)

switch methodOpt    
    case 1 %SPARSE ARRAY BASED (more efficient for large arrays)        
        nnzQ=nnz(logic_ML_not_nan);
        Iq=labelIndices(ind2);
        Jq=1:nnzQ;
        Sq=M(logic_ML_not_nan);
        sizQ=[numLabels,nnzQ];
        nnzQ=nnz(logic_ML_not_nan);
        Q=sparse(Iq,Jq,Sq,sizQ(1),sizQ(2),nnzQ);
        
        Q_sum=full(sum(Q,2));
        Q_voxelCount=full(sum(spones(Q),2));
        Q_mean=Q_sum./Q_voxelCount;
    case 2 %REGIONPROPS BASED
        M_id=nan(size(ML));
        M_id(logic_ML_not_nan)=labelIndices(ind2);
        M_id(isnan(M_id))=0;
        A=regionprops(M_id,M,'MeanIntensity');
        Q_mean=[A.MeanIntensity];
        Q_mean=Q_mean(:);
end

%%

varargout{1}=Q_mean;
if nargout==2
    varargout{2}=labelSet;
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
