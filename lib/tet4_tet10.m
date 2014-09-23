function [TET10,V10,VX10C]=tet4_tet10(TET4,V4,VXC)

% function [TET10,V10]=tet4_tet10(TET4,V4)
% 
% This function converts 4 node (e.g. linear) tetrahedral elements into 10
% node (e.g. quadratic) tetrahedral elements compatible with FEBio. 
%
% 09/12/2013, Kevin Mattheus Moerman
% kevinmoerman@hotmail.com.com
%----------------------------------------------------------------------
%%
%Collect nodes
V_1_4=V4;
V_5 =0.5.*(V4(TET4(:,1),:)+V4(TET4(:,2),:));
V_6 =0.5.*(V4(TET4(:,2),:)+V4(TET4(:,3),:));
V_7 =0.5.*(V4(TET4(:,3),:)+V4(TET4(:,1),:));
V_8 =0.5.*(V4(TET4(:,1),:)+V4(TET4(:,4),:));
V_9 =0.5.*(V4(TET4(:,2),:)+V4(TET4(:,4),:));
V_10=0.5.*(V4(TET4(:,3),:)+V4(TET4(:,4),:));
V10=[V_1_4;V_5;V_6;V_7;V_8;V_9;V_10];

%Define elements
numTets=size(TET4,1);
num_V_1_4=size(V_1_4,1);
indTets=(1:numTets)';
TET10=[TET4...                            % 1-4
       indTets+num_V_1_4+(numTets*(1-1))... % 5
       indTets+num_V_1_4+(numTets*(2-1))... % 6
       indTets+num_V_1_4+(numTets*(3-1))... % 7
       indTets+num_V_1_4+(numTets*(4-1))... % 8
       indTets+num_V_1_4+(numTets*(5-1))... % 9
       indTets+num_V_1_4+(numTets*(6-1))];  % 10
   
%Removing double points
% [~,ind_uni_1,ind_uni_2]=uniqueEps(Vs,'rows',4); %N.B.
[~,ind_uni_1,ind_uni_2]=unique(pround(V10,5),'rows');
V10=V10(ind_uni_1,:);
TET10=ind_uni_2(TET10); %Changing indices in faces matrix

if size(TET4,1)==1 %Transpose in this case due to MATLAB behaviour
    TET10=TET10';
end

%Derive VX10
if ~isempty(VXC)    
    VX10C=VXC;
    for q=1:1:numel(VXC)
        VX=double(VXC{q});
        VX_1_4=VX;
        VX_5 =0.5.*(VX(TET4(:,1),:)+VX(TET4(:,2),:));
        VX_6 =0.5.*(VX(TET4(:,2),:)+VX(TET4(:,3),:));
        VX_7 =0.5.*(VX(TET4(:,3),:)+VX(TET4(:,1),:));
        VX_8 =0.5.*(VX(TET4(:,1),:)+VX(TET4(:,4),:));
        VX_9 =0.5.*(VX(TET4(:,2),:)+VX(TET4(:,4),:));
        VX_10=0.5.*(VX(TET4(:,3),:)+VX(TET4(:,4),:));
        VX10=[VX_1_4;VX_5;VX_6;VX_7;VX_8;VX_9;VX_10];        
        VX10C{q}=VX10(ind_uni_1,:);
    end
else
    VX10C={};
end




   
