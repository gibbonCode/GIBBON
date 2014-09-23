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