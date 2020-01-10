function varargout=maxDist(varargin)

% function [D1,maxIND]=maxDist(V1,V2,maxVarSize,selfAvoid)

%% Parse input
if nargin<2
    error('Insufficient input arguments');
end

V1=varargin{1};
V2=varargin{2};
switch nargin
    case 2        
        maxVarSize=[]; %Empty will force calucation below
        selfAvoid=0; 
    case 3
        maxVarSize=varargin{3};
        selfAvoid=0; 
    case 4
        maxVarSize=varargin{3};
        selfAvoid=varargin{4};        
end

%Get max variable size available        
if isempty(maxVarSize)
    [numFreeBytes]=freeMemory;
    maxVarSize=numFreeBytes/2;
end

%Derive class dependent variable size
[~,b1]=maxnumel(V1(1));
[~,b2]=maxnumel(V2(1));
b=max([b1 b2]);
numelVar=numel(V1)*numel(V2);
varSize=numelVar*b;

numSteps=ceil(varSize/maxVarSize);
indSteps=round(linspace(0,size(V1,1),numSteps));
indSteps=sort(unique(indSteps));
numSteps=numel(indSteps);

if numSteps>1 %In steps
    D1=zeros(size(V1,1),1);
    maxIND=zeros(size(V1,1),1);
    for q=1:1:numSteps-1
        v1=V1(indSteps(q)+1:indSteps(q+1),:);
        try 
            d=dist(v1,V2'); %dist from Neural network toolbox
        catch
            d=distND(v1,V2); %GIBBON's dist function
        end
        if selfAvoid
            %Set "diagonal" to something too large so self is avoided in
            %minimum 
            I=1:size(v1,1);
            J=indSteps(q)+1:indSteps(q+1);
            ind=sub2ind(size(d),I,J); %Indices of selfies            
            d(ind)=1+max(d(:)); %Overide selfies
        end
        
        [max_d,max_ind]=max(d,[],2);
        D1(indSteps(q)+1:indSteps(q+1))=max_d;
        maxIND(indSteps(q)+1:indSteps(q+1))=max_ind;        
    end
else %In one go
    try
        D=dist(V1,V2'); %dist from Neural network toolbox
    catch
        D=distND(V1,V2); %GIBBON's dist function
    end    
    if selfAvoid
        %Set "diagonal" to something too large so self is avoided in
        %minimum 
        L=eye(size(D))>0;
        D(L)=1+max(D(:));
    end
    [D1,maxIND]=max(D,[],2);         
    D1=D1(:);
    maxIND=maxIND(:);
end

switch nargout
    case 1
        varargout{1}=D1;
    case 2
        varargout{1}=D1;
        varargout{2}=maxIND;
    otherwise
        error('wrong number of output arguments');
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
