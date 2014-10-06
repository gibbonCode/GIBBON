function [IND]=box_indices(siz)

% function [IND]=box_indices(siz)
% ------------------------------------------------------------------------
% This function returns the indices of the boundary entries (e.g. pixels or
% voxels) for an array with size siz. The boundary entries form a box at
% the edge of the matrix. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

numDim=numel(siz); 
logicSides=false(siz); 
permuteOrder=1:numDim;
for q=1:1:numDim
    if q>1
        logicSides=permute(logicSides,permuteOrder);
    end
    
    switch numDim
        case 2
            logicSides(1,:)=1;
            logicSides(end,:)=1;
        case 3
            logicSides(1,:,:)=1;
            logicSides(end,:,:)=1;
        case 4
            logicSides(1,:,:,:)=1;
            logicSides(end,:,:,:)=1;
    end
       
    if q>1
        logicSides=ipermute(logicSides,permuteOrder);
    end
    
    permuteOrder=[permuteOrder(2:end) permuteOrder(1)]; %Shift order
   
end

IND=find(logicSides);