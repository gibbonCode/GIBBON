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
