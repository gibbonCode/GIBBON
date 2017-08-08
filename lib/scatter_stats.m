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
