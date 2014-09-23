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