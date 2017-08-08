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
