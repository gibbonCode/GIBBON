function [L]=isapprox(varargin)
% function [L]=isapprox(A,B,tolLevel)
% ------------------------------------------------------------------------
% Check if A is approximately equal to B (to within tolLevel). 
%
% ------------------------------------------------------------------------
%%

switch nargin
    case 1
        A=varargin{1};
        B=zeros(size(A));
        tolLevel=eps(1);
    case 2
        A=varargin{1};
        B=varargin{2};
        tolLevel=eps(1);
    case 3
        A=varargin{1};
        B=varargin{2};
        tolLevel=varargin{3};
end

if numel(B)==1
    B=B.*ones(size(A));
end

L=abs(A-B)<tolLevel;

