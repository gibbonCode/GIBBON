function [F,V]=quadPlate(varargin)

% function [F,V]=quadPlate(plateDim,plateEl)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2020/07/27
%------------------------------------------------------------------------

%% Parse input 

switch nargin
    case 1
        boxDim=varargin{1};
        boxEl=[10 10 10];        
    case 2
        boxDim=varargin{1};
        boxEl=varargin{2};                
end

%%

dX=boxDim(1); 
dY=boxDim(2); 
nX=boxEl(1); 
nY=boxEl(2);

[X,Y] = meshgrid(linspace(-dX/2,dX/2,nX+1),linspace(-dY/2,dY/2,nY+1)); %Grid

[F,V] = grid2patch(X,Y,zeros(size(X))); %Convert to patch data (quadrilateral faces)

