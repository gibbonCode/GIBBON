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

if ~isempty(Q)
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
        
        %Check for non-unique points        
        numDigitsMerge=6-numOrder(mean(sum(diff(Vg,[],1).^2,2)));
        [~,indUni,~]=unique(pround(Vg,numDigitsMerge),'rows');
        logicKeep=false(size(Vg,1),1);
        logicKeep(indUni)=1;
        Vg=Vg(logicKeep,:);
        
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
else
    C={};
    cSiz=[];
    cLevel=[];
end
varargout{1}=C;
varargout{2}=cSiz;
varargout{3}=cLevel;

end

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
