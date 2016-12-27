function [varargout]=gcontour(varargin)

% function [C,cSiz,cLevel]=gcontour(X,Y,Z,k,pointSpacing,resampleMethod)
% ------------------------------------------------------------------------
%
% This function computes contour data for the data Z 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/10/18
%
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 4
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        k=varargin{4};
        pointSpacing=[]; 
    case 5
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        k=varargin{4};
        pointSpacing=varargin{5};
        resampleMethod='pchip';
    case 6
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        k=varargin{4};
        pointSpacing=varargin{5};
        resampleMethod=varargin{6};
    otherwise
        error('Wrong number of input arguments');
end
        
%%
%Use contourc to derive contours
x=X(1,:)';
y=Y(:,1);
v=abs([x(2)-x(1) y(2)-y(1)]);

%Convert if it is a logical array
if isa(Z,'logical')
    Z=double(Z);
end

if numel(k)==1
    Q = contourc(x,y,Z,[k k]);
else
    Q = contourc(x,y,Z,k);
end

%Get contour indices
c=0;
ind2=0;
n=0; 
IND1=[];
IND2=[];
while 1       
    numPoints=Q(2,ind2+1);
    n=n+numPoints;
    ind1=1+(ind2+1);
    ind2=1+n+c;    
    c=c+1;     
    IND1(c)=ind1; %#ok<AGROW>
    IND2(c)=ind2; %#ok<AGROW>    
    if ind2==size(Q,2)
        break        
    end
end

numContours=numel(IND1); 
C=cell(numContours,1);
cSiz=zeros(numContours,1);
cLevel=zeros(numContours,1);
for q=1:1:numContours
    
    Vg=[Q(1,IND1(q):IND2(q))' Q(2,IND1(q):IND2(q))']; %Get coordinates
    cLevel(q)=Q(1,IND1(q)-1); %Get level
    
    if size(Vg,1)>1
        %Check if contour is closed loop or not
        d=sqrt(sum((Vg(1,:)-Vg(end,:)).^2)); %Distance between first and last point
        if d<(min(v)/10) %Smaller than 1/10th of a pixel dimension
            Vg=Vg(1:end-1,:);
            closeLoopOpt=1;
        else
            closeLoopOpt=0;
        end
        
        %Upsample if desired
        if ~isempty(pointSpacing)
            D=pathLength(Vg);               
            [Vg] = evenlySampleCurve(Vg,n,resampleMethod,closeLoopOpt);
        end        
    end
    
    %Store contour in cell array
    C{q}=Vg;
    cSiz(q)=size(Vg,1);
end

varargout{1}=C;
varargout{2}=cSiz;
varargout{3}=cLevel;

end

