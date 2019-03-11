function [IMAGE_OBJECTS]=imobjects(IND,IND_start,IND_adj)

% function [IMAGE_OBJECTS]=imobjects(IND,IND_start,IND_adj)
% ------------------------------------------------------------------------
% This function finds objects in the image M using adjacency. 
%   [IMAGE_OBJECTS]=imageobjects(M,IND)
% The seach will only apply to the voxels specified by the linear indices
% of the vector IND. IND could be generated to exclude voxels in the
% image of a certain intensity. 
%   [IMAGE_OBJECTS]=imageobjects(M,IND,IND_start)
% The seach will only apply to the voxels specified by the linear indices
% of the vector IND. However groups are only grown at the sites specified
% by IND_start. IND_start could be a product of an image processing
% procedure, a starting guess to the objects of interest.
% Together IND and IND_start allow for control of intensity and
% aproximate location of the objects of interest. 
% Groups are stored as vectors containing linear indices in the cell
% IMAGE_OBJECTS.
%
% 07/04/2011
% ------------------------------------------------------------------------

%%

IND_start_grouped=IND;
IND_start_grouped(~ismember(IND,IND_start))=0;

IND_adj_grouped=IND_adj;
IMAGE_OBJECTS={};
group_no=1;
while any(IND_start_grouped)    
    new_start_index=find(IND_start_grouped>0,1); % Next non-zero entry in IND_start_grouped
    new_points=IND_adj(new_start_index,:); new_points=new_points(new_points>0);
    object_group=[IND_start_grouped(new_start_index); new_points(:)];    
    num_group=1;
    group_full=0;
    while group_full==0
        object_group=unique(object_group); %Removing doubles NB and sorts!
        L_adj=any(ismembc(IND_adj_grouped,object_group),2); %Second input needs to be sorted, see help ismember
        IND_adj_found=IND_adj_grouped(L_adj,:); %Linear indices for all voxels touching group elements
        IND_adj_found=IND_adj_found(IND_adj_found>0);
        object_group=[object_group;IND_adj_found(:)]; %Adding found to group 
        
        object_group=unique(object_group); %Removing doubles NB and sorts!
        
        IND_adj_grouped(L_adj,:)=0; %Setting indices found to zero in adjacency matrix        
        if num_group==numel(object_group); %BREAK for while loop
            
            IMAGE_OBJECTS{group_no}=object_group; %Adding found elements to cell            
            L_start=ismembc(IND_start_grouped,object_group);
            IND_start_grouped(L_start)=0; %Setting indices found to zero in start vector
            group_no=group_no+1; %Group counter
            group_full=1; 
        end
        num_group=numel(object_group);
    end
end

%%


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
