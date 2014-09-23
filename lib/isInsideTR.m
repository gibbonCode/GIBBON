function [varargout]=isInsideTR(varargin)

if nargin<2
    error('Not enough input arguments');
end

TR=varargin{1};
numTR=size(TR.ConnectivityList,1);
QP=varargin{2};

if nargin==3
    ti=varargin{3};
else
    ti=1:numTR;
end
    
%Get barycentric coordinates of points
baryCoords=cartesianToBarycentric(TR,ti(:),QP);

logicFoundEnclosing=all(baryCoords>0,2);

varargout{1}=logicFoundEnclosing;
if nargout==2
    varargout{2}=baryCoords;
end



