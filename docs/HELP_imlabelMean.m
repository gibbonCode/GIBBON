%% imlabelMean
% Below is a demonstration of the features of the |imlabelMean| function

%%
clear; close all; clc;

%% Syntax
% |[Q_mean]=imlabelMean(M,ML);|

%% Description 
% This function takes the mean for each of the labeled groups (NaN's
% ignored) in ML according to the intensities in M.

%% Examples 
% 

%%
% Create example labelled image

% Defining a multi boundary set
r=2; %Sphere radius
rc=3; %Central radius
nr=15;
nc=25;
ptype='quad';
[F1,V1]=patchTorus(r,nr,rc,nc,ptype);
[F2,V2]=quadSphere(2,r,2);
V2(:,2)=V2(:,2)-5;
[F3,V3]=quadSphere(2,r/2,2);
V3(:,2)=V3(:,2)-5;
[F4,V4]=quadSphere(3,r/2,2);
V4(:,1)=V4(:,1)+2;
V4(:,2)=V4(:,2)+2;

[F,V,C]=joinElementSets({F1,F2,F3,F4},{V1,V2,V3,V4});

%Convert to label image
[~,~,ML]=patch2Im(F,V,[],0.1*ones(size(F,1),1));

%Create image to compute stats for these labelled regions
M=nan(size(ML));
c=0;
for q=min(ML(~isnan(ML))):1:max(ML(~isnan(ML)))
    if nnz(ML==q)>0
        c=c+1;
    end
M(ML==q)=c+q/10; %Force data to be <group number>.<label>
end

%%
% Compute the mean for each region

[Q_mean,labelSet]=imlabelMean(M,ML)

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
