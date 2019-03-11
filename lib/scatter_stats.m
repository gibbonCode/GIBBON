function [STATS]=scatter_stats(D)

%Fitting Gaussian mixture (1) model
GM=gmdistribution.fit(D,1); 
MU=GM.mu; %Means
MU_total=rms(MU); %Total mean
VAR=GM.Sigma; %Variance
STD=sqrt(VAR); %Standard deviations
[VAR_vec,VAR_eig] = eig(VAR); %Get principal components
STD_total=sqrt(trace(VAR_eig)./3); %Standard deviations

%Creating output structure
STATS.GM=GM;
STATS.MU=MU;
STATS.MU_total=MU_total;
STATS.VAR=VAR;
STATS.STD=STD;
STATS.STD_total=STD_total;

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
