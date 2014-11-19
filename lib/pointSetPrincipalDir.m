function [varargout]=pointSetPrincipalDir(X)


%% 

%Centre on own mean
MU=mean(X,1); %Point set mean
X=X-MU(ones(size(X,1),1),:); %Centre points around mean

%Compute singular value decomposition to get principal directions
[U,S,V]=svd(X,0); 

%% Collect output
switch nargout
    case 1
        varargout{1}=V; 
    case 2
        varargout{1}=V;
        varargout{2}=S; 
    case 3
        varargout{1}=V;
        varargout{2}=S;
        varargout{3}=U; 
end

