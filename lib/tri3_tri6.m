function [TRI6,V6,VX6C]=tri3_tri6(TRI3,V3,VXC)

% function [TRI6,V6,VX6C]=tri3_tri6(TRI3,V3,VXC)
% 
% This function converts 3 node (e.g. linear) triangular elements into 6
% node (e.g. quadratic) triangular elements compatible with FEBio. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 09/12/2013
%----------------------------------------------------------------------
%%

%Collect nodes
V_1_3=V3;
V_4 =0.5.*(V3(TRI3(:,1),:)+V3(TRI3(:,2),:));
V_5 =0.5.*(V3(TRI3(:,2),:)+V3(TRI3(:,3),:));
V_6 =0.5.*(V3(TRI3(:,3),:)+V3(TRI3(:,1),:));
V6=[V_1_3;V_4;V_5;V_6;];

%Derive VX10
if ~isempty(VXC)    
    VX6C=VXC;
    for q=1:1:numel(VXC)
        VX=VXC{q};
        VX_1_3=VX;
        VX_4 =0.5.*(VX(TRI3(:,1),:)+VX(TRI3(:,2),:));
        VX_5 =0.5.*(VX(TRI3(:,2),:)+VX(TRI3(:,3),:));
        VX_6 =0.5.*(VX(TRI3(:,3),:)+VX(TRI3(:,1),:));
        VX6=[VX_1_3;VX_4;VX_5;VX_6;];
        VX6C{q}=VX6;
    end
else
    VX6C={};
end

%Define elements
numTris=size(TRI3,1);
num_V_1_3=size(V_1_3,1);
indTris=(1:numTris)';
TRI6=[TRI3...                            % 1-3
       indTris+num_V_1_3+(numTris*(1-1))... % 4
       indTris+num_V_1_3+(numTris*(2-1))... % 5
       indTris+num_V_1_3+(numTris*(3-1))... % 6       
       ];  % 10
   
%Removing double points
numDigitKeep=6; 
try
    [~,ind_uni_1,ind_uni_2]=unique(round(V6,numDigitKeep,'significant'),'rows');
catch
    [~,ind_uni_1,ind_uni_2]=unique(sround(V6,numDigitKeep),'rows');
end
V6=V6(ind_uni_1,:);
TRI6=ind_uni_2(TRI6); %Changing indices in faces matrix

   
