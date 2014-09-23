function [F,C]=hex2patch(E,C)

%%
% ------------------------------------------------------------------------
% This function generates FACE data for patch grapics disply of hexahedral
% elements
%
%N.B. This function assumes a node order whereby 1 2 3 4 described the top
%face and 5 6 7 8 the bottom face
% ------------------------------------------------------------------------


%% 

% F =[E(:,[1 2 3 4]);... %top
%     E(:,[5 6 7 8]);... %bottom
%     E(:,[1 2 6 5]);... %side 1
%     E(:,[3 4 8 7]);... %side 2
%     E(:,[2 3 7 6]);... %front
%     E(:,[1 4 8 5]);]; %back

F =[E(:,[4 3 2 1]);... %top
    E(:,[5 6 7 8]);... %bottom
    E(:,[1 2 6 5]);... %side 1
    E(:,[3 4 8 7]);... %side 2
    E(:,[2 3 7 6]);... %front
    E(:,[5 8 4 1]);]; %back

C=repmat(C,6,1);

%%
end