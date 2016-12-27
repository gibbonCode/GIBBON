function [M,linearInd]=subImage(M,nSub)

% function M=subImage(M,nSub)
% ------------------------------------------------------------------------
%
%This function refines the input image M by splitting its voxels in half during nSub times
%which contains all folders and sub-folders within the folder specified by
%the input pathName.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 18/04/2013
%------------------------------------------------------------------------

%%
%TO DO: ADD WARNINGS AND INPUT PARSING HERE
% nSub=round(nSub);

%%

if numel(nSub)==1   
    nSub=nSub.*ones(1,3);
end

stepSize_I=1/nSub(1);
stepSize_J=1/nSub(2);
stepSize_K=1/nSub(3);

%Define image coordinates of new voxel centres within image coordinate
%system of original coordinate mesh 
I_range=linspace((0.5+(stepSize_I/2)),(size(M,1)-(stepSize_I/2))+0.5,size(M,1)*nSub(1));
J_range=linspace((0.5+(stepSize_J/2)),(size(M,2)-(stepSize_J/2))+0.5,size(M,2)*nSub(2));
K_range=linspace((0.5+(stepSize_K/2)),(size(M,3)-(stepSize_K/2))+0.5,size(M,3)*nSub(3));  

[In,Jn,Kn]=ndgrid(I_range,J_range,K_range); %Coordinate matrices

%Rounding the coordinates so they "snap" to the mother voxel indices they
%are found in.
In=round(In); Jn=round(Jn); Kn=round(Kn);

%Convert these subscript indiced to linear indices in the image
linearInd =sub2ind(size(M),In,Jn,Kn);

%Get intensity values
Mc=M(:); %Image as column
M=Mc(linearInd); %The refined image

