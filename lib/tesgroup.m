function varargout=tesgroup(F)

% function [G,G_iter]=tesgroup(F)
% ------------------------------------------------------------------------
%
% This function finds groups in the tesselation defined by F. F may
% represent patch type faces or for instances node indices for
% tetrehedrons, hexahedrons. Row entries in F (e.g. tetrahedron vertex
% indices) which are "connected" (sharing vertex indices with other row
% entries in F) are grouped together. The output G is a logic matrix of
% size(F,1) rows and "number of groups" columns. Each column represents a
% group and ones appear in the column for each face belonging to the group.
%
% % %EXAMPLE:
% clear all; close all; clc;
% 
% 
% %% Simulating some isosurface data
% [X,Y,Z]=meshgrid(linspace(-5,5,35));
% phi=(1+sqrt(5))/2;
% M=2 - (cos(X + phi*Y) + cos(X - phi*Y) + cos(Y + phi*Z) + cos(Y - phi*Z) + cos(Z - phi*X) + cos(Z + phi*X));
% M=M./max(M(:)); 
% [F,V] = isosurface(X,Y,Z,M,0.1);
% 
% %% Normal isosurface plot showing seperate patch objects
% figure; 
% h=patch('faces',F,'vertices',V); 
% set(h,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
% view(3);light; grid on; axis vis3d;
% 
% %% Iso surface plots showing grouped patch objects
% 
% G=tesgroup(F); %Logic array for patch groups
% pcolors=jet(size(G,2));
% figure; 
% for i=1:1:size(G,2);    
%     hg=patch('faces',F(G(:,i),:),'vertices',V); %Plotting individual group
%     set(hg,'FaceColor',pcolors(i,:),'EdgeColor','none','FaceAlpha',0.8);
% end
% view(3);light; grid on; axis vis3d;
% colormap(pcolors); colorbar; caxis([0 size(G,2)]);
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/07/2010
%------------------------------------------------------------------------

IND_F=(1:1:size(F,1))';
IND_F_search=IND_F;
G=false(size(F,1),1);
L=false(size(F,1),1);
L_previous=false(size(F,1),1);

G_ind=nan(size(F,1),1);
v_search=[ ];
done=0;
num_v_search=0;
group_found=1;
group_n=0;
q=1; %Counter
while done==0; 
    L=false(size(F,1),1);
    if group_found==1;
        indNext=find(IND_F_search>0,1); %next un-grouped element
        L(indNext)=1;
        IND_F_search(L)=0; %Setting found to zero
        group_found=0;
    else
        L = any(ismember(F,v_search), 2);
        IND_F_search(L)=0; %Setting found to zero
    end
    v_new=F(L,:); 
    v_search=unique([v_search; v_new(:)]); %Growing number of search vertices    
    
    G_ind(L&~L_previous)=q;
    L_previous=L; 
    
    if numel(v_search)==num_v_search; %If the group has not grown
        group_found=1;        
        group_n=group_n+1;
        G(:,group_n)=L;        
        v_search=[ ];
    end
    
    if all(IND_F_search==0);
        done=1;
        group_found=1;        
        group_n=group_n+1;
        if ~all(any(G,2))%if not all points are grouped keep remainder
            G(:,group_n)=L;
        end        
        v_search=[ ];
    end    
    num_v_search=numel(v_search);
    
    q=q+1; %Increment counter
end

%%

G_iter=G_ind(:,ones(size(G,2),1));
G_iter(~G)=NaN;
G_iter_min=nanmin(G_iter,[],1);
G_iter=G_iter-G_iter_min(ones(size(G_iter,1),1),:)+1;

%% Prepare output
switch nargout
    case 1
        varargout{1}=G;
    case 2
        varargout{1}=G;
        varargout{2}=G_iter;
end

%%
% function G=tesgroup(F)
% 
% % function G=tesgroup(F)
% % ------------------------------------------------------------------------
% %
% % This function finds groups in the tesselation defined by F. F may
% % represent patch type faces or for instances node indices for
% % tetrehedrons, hexahedrons. Row entries in F (e.g. tetrahedron vertex
% % indices) which are "connected" (sharing vertex indices with other row
% % entries in F) are grouped together. The output G is a logic matrix of
% % size(F,1) rows and "number of groups" columns. Each column represents a
% % group and ones appear in the column for each face belonging to the group.
% %
% % % %EXAMPLE:
% % clear all; close all; clc;
% % 
% % 
% % %% Simulating some isosurface data
% % [X,Y,Z]=meshgrid(linspace(-5,5,35));
% % phi=(1+sqrt(5))/2;
% % M=2 - (cos(X + phi*Y) + cos(X - phi*Y) + cos(Y + phi*Z) + cos(Y - phi*Z) + cos(Z - phi*X) + cos(Z + phi*X));
% % M=M./max(M(:)); 
% % [F,V] = isosurface(X,Y,Z,M,0.1);
% % 
% % %% Normal isosurface plot showing seperate patch objects
% % figure; 
% % h=patch('faces',F,'vertices',V); 
% % set(h,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
% % view(3);light; grid on; axis vis3d;
% % 
% % %% Iso surface plots showing grouped patch objects
% % 
% % G=tesgroup(F); %Logic array for patch groups
% % pcolors=jet(size(G,2));
% % figure; 
% % for i=1:1:size(G,2);    
% %     hg=patch('faces',F(G(:,i),:),'vertices',V); %Plotting individual group
% %     set(hg,'FaceColor',pcolors(i,:),'EdgeColor','none','FaceAlpha',0.8);
% % end
% % view(3);light; grid on; axis vis3d;
% % colormap(pcolors); colorbar; caxis([0 size(G,2)]);
% %
% % Kevin Mattheus Moerman
% % kevinmoerman@hotmail.com
% % 15/07/2010
% %------------------------------------------------------------------------
% 
% IND_F=(1:1:size(F,1))';
% IND_F_search=IND_F;
% 
% G=false(size(F,1),1);
% v_search=[ ];
% L=ones(size(IND_F));
% done=0;
% num_v_search=[ ];
% group_found=1;
% group_n=0;
% while done==0;
%     if group_found==1;
%         L=find(IND_F_search>0,1); %next un-grouped triangle
%         v_new=F(L,:); v_new=v_new(:);
%         v_search=[v_search; v_new]; v_search=unique(v_search(:)); %Growing number of search vertices
%         
%         group_found=0;
%     else
%         L = any(ismember(F,v_search), 2);
%         IND_F_search=IND_F_search.*(L==0); %Setting found to zero
%         v_new=F(L,:); v_new=v_new(:);
%         v_search=[v_search; v_new]; v_search=unique(v_search(:)); %Growing number of search vertices
%     end
% 
%     if numel(v_search)==num_v_search; %If the group has not grown
%         group_found=1;
%         group_n=group_n+1;
%         G(:,group_n)=L;
%         v_search=[ ];
%     end
% 
%     num_v_search=numel(v_search);
%     if all(IND_F_search==0);
%         done=1;
%         group_found=1;
%         group_n=group_n+1;
%         if any(G)==0
%             G(:,group_n)=L;
%         end
%         v_search=[ ];
%     end
% 
% end
 
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
