function [C]=gray2RGBColorMap(varargin)

% function [C]=gray2RGBColorMap(G,cMap,cLim)
% ------------------------------------------------------------------------
% This function maps the monochrome data vector G, containing for instance 
% gray scale values, to a colormapped RGB array C using the colormap cMap
% and the color limits cLim. If C is not provided the jet colormap is used.
% If cLim is not provided or empty the colorlimits are adjusted according
% to the extreme values in G. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        G=varargin{1};
        cMap=jet(250);
        cLim=[];
    case 2
        G=varargin{1};
        cMap=varargin{2};
        cLim=[];
    case 3
        G=varargin{1};
        cMap=varargin{2};
        cLim=varargin{3};
end

if isempty(cLim)
   cLim=[min(G(:)) max(G(:))]; 
end

G=G(:); %G in column form

%% Create colormapped RGB array

%Snap to color limits
G(G<cLim(1))=cLim(1);
G(G>cLim(2))=cLim(2);

%Normalized to range [0-1]
G=G-min(G(:)); 
max_G=max(G(:));
if max_G>0
    G=G./max(G(:));
end
        
numColorLevels=size(cMap,1); %Number of color levels
linearIndex_G=1+round(G.*(numColorLevels-1));
C=cMap(linearIndex_G,:);

