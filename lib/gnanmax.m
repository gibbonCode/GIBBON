function [varargout]=gnanmax(varargin)
% function [y,ind]=gnanmax(A,[],dimDir)
%------------------------------------------------------------------------
% this function is a replacement for nanmax which is part of the
% statistics and machine learning toolbox. 
% The use of max with the omitnan flag is used or if this does not work
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
    if nargout==2
        [y,ind] = max(varargin{:},'omitnan');
    else
        y = max(varargin{:},'omitnan');
    end
catch
    %Use custom version     
    A=varargin{1}; %First array
    B=varargin{2}; %Second array or number or empty
    
    if nargin==3
        dimDir=varargin{3};
    else
        dimDir=1;
    end
    
    LA=isnan(A);    
    minVal=min(A(~LA));
    
    if ~isempty(B)
        LB=isnan(B);
        minVal=minVal+min(B(~LB));
    end
    if isempty(minVal)
        minVal=0;
    end
    
    AV=A;
    AV(LA)=minVal-1;    
    if ~isempty(B)
        if nargout==2
            error('If two matrices are compared two output arguments are not supported.');
        end
        BV=B;
        BV(LB)=1+minVal;                
        y=max(AV,BV); 
        y(LA&LB)=NaN;
        y(~LA&LB)=A(~LA&LB);
        y(LA&~LB)=B(LA&~LB);
    else
        [y,ind]=max(AV,[],dimDir);            
        y(all(LA,dimDir))=NaN; %replace all NaN directions
    end
    
end

varargout{1}=y;
if nargout==2
    varargout{2}=ind;
end

end

