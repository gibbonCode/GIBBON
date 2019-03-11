function [P2i]=pointSetDistMap(P1,P1i,P2,Fw)

% function [P2i]=pointSetDistMap(P1,P1i,P2,Tw)
% ------------------------------------------------------------------------
% This function aims to use a simple distance based weighting to map the
% points P1i (ideally inside or on the cloud defined by) P1 within the
% point cloud defined by P2. A point-to-point correspondence between the
% sets P1 and P2 is assumed such that point displacements are defined
% simple as P2-P1. 
% A weights matrix is derived based on: one over the distance from each
% point in P1 to each point in P1i. Hence this function is inefficient for
% large data sets. The weights are used to set the displacement of the
% points in P1i. The closer a point in P1i is to a certain point in P1 the
% higher its associated weight is for that point and thus the more its
% displacement is based on the displacement of this close point. For zero
% distances points in P1i actually coincide with a point in P1, in this
% case the weights cause them to move to the corresponding point in P2
% exactly. The factor Fw is used to alter the weighting such that closer
% points are favoured. The weights matrix W is adjusted in this fashion
% W=W.^Fw. Whereby the initial W is simply distance based and so Fw=2
% results in a squared distances based mapping. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 12/08/2013
%------------------------------------------------------------------------

%% Derive distance based weights matrix 

%Calulate displacement
pointDisp=P2-P1;

DW=dist(P1,P1i'); %Distances

%Express distance in termps of the average point spacing
DP=pathLength(P1);
pointSpacing=max(DP(:))./size(P1,1);
DW=DW./pointSpacing; %Distances in multiples of pointSpacing

W=1./DW; %The weights matrix initialized as the inverse of the distance metric

%Treating possible INF entries
L=isinf(W);
W(:,any(L,1))=0; %Set rows with INF values to zero
W(L)=1; %Set INF locations back to 1

%Raise power of the weights to bias closer points
W=W.^Fw;

%Normalise weights to sum
sumW=sum(W,1);
W=W./sumW(ones(1,size(W,1)),:);

%% Map points

P2i=P1i; %Initialize output
for q=1:1:size(P1,2) %Loop for all dimensions
    dispN=pointDisp(:,q*ones(size(W,2),1));
    P2i(:,q)=P1i(:,q)+sum(dispN.*W,1)';
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
