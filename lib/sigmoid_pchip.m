function [s]=sigmoid_pchip(optStruct)

%%

fieldNames={'c1','c2','r1','r2','n','nLin','rMode'};
defaultVal={0,0,0.25,0.25,100,2,1};
for q=1:1:numel(fieldNames)    
    if ~isfield(optStruct,fieldNames{q})
        optStruct.(fieldNames{q})=defaultVal{q};
%     elseif isempty(optStruct.(fieldNames{q}))
%         optStruct.(fieldNames{q})=defaultVal{q};
    end
end

c1=optStruct.c1;
c2=optStruct.c2;
r1=optStruct.r1;
r2=optStruct.r2;
nLin=optStruct.nLin;
n=optStruct.n;
rMode=optStruct.rMode;

%%

if rMode==1
    x1=r1;
    x2=r2;
else
    x1=sqrt((r1^2)/(c1^2+1));
    x2=sqrt((r2^2)/(c2^2+1));
end
y1=c1*x1;
y2=(c2*x2);

x2=1-x2;
y2=1-y2;

x1=linspace(0,x1,nLin);
y1=linspace(0,y1,nLin);
x2=linspace(x2,1,nLin);
y2=linspace(y2,1,nLin);

if numel(n)==1
    xi=linspace(0,1,n);
elseif numel(n)>1
    xi=n;
else %n must be empty
    xi=[]; %Make xi empty will force pp form output
end

if ~isempty(xi)
    yi = pchip([x1(:);x2(:)],[y1(:);y2(:)],xi);
else
    pp = pchip([x1(:);x2(:)],[y1(:);y2(:)]);    
end

if numel(n)==1
    s=[xi(:) yi(:)];
elseif numel(n)>1
    s=yi;
else
    s=pp;
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
