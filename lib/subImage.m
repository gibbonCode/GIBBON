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
