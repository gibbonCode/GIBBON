function [F,V]=coarsepatch(F,V,f)

%%
% Run reducepatch
[F,V] = reducepatch(F,V,f,'fast');

%%
% Avoid/fix reducepatch bugs
[F,V]=triSurfRemoveThreeConnect(F,V); 
[F,V]=mergeVertices(F,V); %Merge vertices
F=patchRemoveCollapsed(F); %remove collapsed (edges)
F=uniqueIntegerRow(F); %Removing double faces
[F,V]=patchCleanUnused(F,V); %Remove unused points