function [varargout]=edgeVec(E,V)

% function [N,Vp]=edgeVec(E,V)
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------

%% Compute edge vectors
Vp=V(E(:,1),:); %Edge vector origin
N=V(E(:,2),:)-Vp; %Edge vector

%% Collect output
varargout{1}=N;
varargout{2}=Vp;

