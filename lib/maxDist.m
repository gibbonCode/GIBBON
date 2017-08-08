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
            %minimum (could use NaN and nanmin but the latter is a toolbox
            %function)
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
        %minimum (could use NaN and nanmin but the latter is a toolbox
        %function) 
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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
