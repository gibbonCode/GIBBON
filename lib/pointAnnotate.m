function [varargout]=pointAnnotate(V,nodeIndices,varargin)

% function [ht]=pointAnnotate(V,nodeIndices,varargin)
% ------------------------------------------------------------------------
% This function annotates the points definced by V using the indices
% nodeIndices. Optional additional inputs are those associated with
% MATLAB's text function. 
% 
% Change log: 
% 2019/08/06 Created
% ------------------------------------------------------------------------

%% Get node indices 
if isempty(nodeIndices)
    nodeIndices=1:1:size(V,1);
end

%% Create text data 
t = sprintfc('%i',nodeIndices); %Text cell for node id's

%% Plot text at coordinates
varargout{1}=text(V(:,1),V(:,2),V(:,3),t,varargin{:}); %Plot text at nodes

