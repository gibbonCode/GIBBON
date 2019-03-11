function [M_low]=voxelate(M,reduction_factor)

% function [M_low]=voxelate(M,reduction_factor)
% ------------------------------------------------------------------------
% This function lowers the resulution of a 3D image based on
% 'reduction_factor'.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 25/05/2008
% ------------------------------------------------------------------------

%%
num_dims=ndims(M);

switch num_dims
    
    case 2
        %Setting up cell shape
        cell_shape_1=reduction_factor*ones(size(M,1)/reduction_factor,1);
        cell_shape_2=reduction_factor*ones(size(M,2)/reduction_factor,1);
        
        %Converting matrix to cell
        M_cell=mat2cell(M, cell_shape_1, cell_shape_2);
        
        %Calculating average for each cell
        M_cell=cellfun(@mean, cellfun(@mean, M_cell, 'UniformOutput',0), 'UniformOutput',0);
        
    case 3
        %Setting up cell shape
        cell_shape_1=reduction_factor*ones(size(M,1)/reduction_factor,1);
        cell_shape_2=reduction_factor*ones(size(M,2)/reduction_factor,1);
        cell_shape_3=reduction_factor*ones(size(M,3)/reduction_factor,1);
        
        %Converting matrix to cell
        M_cell=mat2cell(M, cell_shape_1, cell_shape_2, cell_shape_3);
        
        %Calculating average for each cell
        M_cell=cellfun(@mean, cellfun(@mean, cellfun(@mean, M_cell, 'UniformOutput',0), 'UniformOutput',0), 'UniformOutput',0);
        
end


%Converting cell to matrix
M_low=cell2mat(M_cell);
%% OLD VERSION
%
% % Setting field of view FOV size
% FOV=[size(M,1)*voxeldim_high, size(M,2)*voxeldim_high, size(M,3)*voxeldim_high];
% no_elements=FOV(1)/voxeldim_low;
%
% reduction_factor=size(M,1)/no_elements;
% cell_shape=reduction_factor*ones(size(M,1)/reduction_factor,1);
%
% for k=1:size(M,3)
%     M_slice=M(:,:,k);
%     cell_slice=mat2cell(M_slice, cell_shape, cell_shape);
%     column_means = cellfun(@mean, cell_slice, 'UniformOutput',0);
%     total_means = cellfun(@mean, column_means, 'UniformOutput',0);
%     M_means(:,:,k)=cell2mat(total_means);
% end
%
% i=1;
% for k=1:size(M,3)/reduction_factor
%     M_low(:,:,k)=sum(M_means(:,:,(i:(i+reduction_factor-1))),3) /reduction_factor;
%     i=i+reduction_factor;
% end
%
% %% END
 
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
