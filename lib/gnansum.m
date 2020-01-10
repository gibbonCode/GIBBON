function [y]=gnansum(varargin)
% function [y]=gnansum(A,[],dimDir)
%------------------------------------------------------------------------
% this function is a replacement for nansum which is part of the
% statistics and machine learning toolbox. 
% The use of sum with the omitnan flag is used or if this does not work
% (for old MATLAB versions) a custom implementation is used. 
%
% Change log:
% 2020/01/09 Created
%------------------------------------------------------------------------

%%
if nargin==1
    varargin{2}=[];
end

try 
    %Use 'omitnan' flag    
    if nargin==1
        varargin{2}=[];
    end
    y = sum(varargin{:},'omitnan');    
catch
    A=varargin{1}; 
    A(isnan(A))=0;
    varargin{1}=A;
    y = sum(varargin{:});    
end

end

