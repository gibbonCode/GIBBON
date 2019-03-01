function [varargout]=patchRemoveCollapsed(F)

% function [F,logicKeep]=patchRemoveCollapsed(F)
% ------------------------------------------------------------------------
% 
% ------------------------------------------------------------------------

%% Detect repeated indices for all faces
% Create logic defining collapsed faces

F_sort=sort(F,2); %Sort faces in 2nd direction so 2 1 2 -> 1 2 2
d=diff(F_sort,[],2); %Difference in 2nd direction so 1 2 2 -> 1 0
logicKeep=~any(d==0,2); %Logic for faces without zeros in difference measure of sorted faces

%% Create a new face set
F=F(logicKeep,:); %Selecting faces without repeated indices

%% Collect output

varargout{1}=F;
varargout{2}=logicKeep;

