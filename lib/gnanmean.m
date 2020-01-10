function [y]=gnanmean(varargin)
% function [y]=gnanmean(A,dimDir)
%------------------------------------------------------------------------
% this function is a replacement for nanmean which is part of the
% statistics and machine learning toolbox. 
% The use of mean with the omitnan flag is used or if this does not work
% (for old MATLAB versions) a custom implementation is used. 
%
% Change log:
% 2020/01/09 Created
%------------------------------------------------------------------------

%%
try 
    %Use 'omitnan' flag
    y = mean(varargin{:},'omitnan');
catch
    %Use custom version
    A=varargin{1};
    if nargin==2
        dimDir=varargin{2};
    else
        dimDir=1;
    end
    L=isnan(A);
    Q=sum(~L,dimDir);
    B=A;
    B(L)=0;
    S=sum(B,dimDir);
    y=S./Q;
    y(all(L,dimDir))=NaN;
end

end

