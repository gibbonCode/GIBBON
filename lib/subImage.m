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
