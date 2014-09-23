function [b,m,n,c,IND_MAP]=unique_map(varargin)

[b,m,n]=unique(varargin{1:end});
IND_MAP=sparse(1:numel(n),n,1:numel(n),numel(n),max(n(:)),numel(n));
c=full(sum(IND_MAP>0,1))';