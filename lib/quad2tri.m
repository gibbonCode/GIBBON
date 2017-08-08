function [varargout]=quad2tri(varargin)

% function [Ft,Vt,Ct]=quad2tri(Fq,Vq,triType,Cq)
% ------------------------------------------------------------------------
% This function converts quadrilateral patch data to triangular patch data.
% 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/11/03
%------------------------------------------------------------------------


%% Parse input
switch nargin
    case 2
        Fq=varargin{1};
        Vq=varargin{2};
        triType=[];
        Cq=[];
    case 3        
        Fq=varargin{1};
        Vq=varargin{2};
        triType=varargin{3};
        Cq=[];
    case 4        
        Fq=varargin{1};
        Vq=varargin{2};
        triType=varargin{3};
        Cq=varargin{4};        
    otherwise
        error('False number of input arguments');
end

if isempty(triType)
    triType='a';
end

%%
switch triType
    case 'b' %backward slash
        Ft=[Fq(:,[1 2 3]);Fq(:,[3 4 1]);];
        Vt=Vq;
        Ct=repmat(Cq,[2,1]);
    case 'f' %Forward slash
        Ft=[Fq(:,[1 2 4]);Fq(:,[2 3 4]);];        
        Vt=Vq;        
        Ct=repmat(Cq,[2,1]);
    case 'e'        
        d1=sum((Vq(Fq(:,1),:)-Vq(Fq(:,3),:)).^2,2);
        d2=sum((Vq(Fq(:,2),:)-Vq(Fq(:,4),:)).^2,2);
        L1=(d1<d2);
        Ft=[Fq(L1,[1 2 3]);Fq(L1,[3 4 1]); Fq(~L1,[1 2 4]);Fq(~L1,[2 3 4])];
        Vt=Vq;
        if ~isempty(Cq)
            Ct=[Cq(L1,:); Cq(L1,:); Cq(~L1,:); Cq(~L1,:)];
        else
            Ct=Cq;
        end
    case 'a'        
        F1=[Fq(:,[1 2 3]);Fq(:,[3 4 1])];
        F2=[Fq(:,[1 2 4]);Fq(:,[2 3 4])];
        [a1]=patchEdgeAngles(F1,Vq);
        a1=reshape(a1,[size(Fq,1),6]);
        d1=sum(abs(a1-(pi/3)),2);
        
        [a2]=patchEdgeAngles(F2,Vq);
        a2=reshape(a2,[size(Fq,1),6]);
        d2=sum(abs(a2-(pi/3)),2);
        L1=(d1<d2);
        Ft=[Fq(L1,[1 2 3]);Fq(L1,[3 4 1]); Fq(~L1,[1 2 4]);Fq(~L1,[2 3 4])];
        Vt=Vq;
        if ~isempty(Cq)
            Ct=[Cq(L1,:); Cq(L1,:); Cq(~L1,:); Cq(~L1,:)];
        else
            Ct=Cq;
        end
    case 'x' %Cross type
        Vm=zeros(size(Fq,1),size(Vq,2));
        for q=1:1:size(Vq,2)
            X=Vq(:,q);
            FX=X(Fq);
            if size(Fq,1)==1 %Treat special case of single face
                FX=FX';
            end
            Vm(:,q)=mean(FX,2);
        end
        %Join point sets
        Vt=[Vq;Vm];
        
        indVm=(size(Vq,1)+1):size(Vt,1);
        %Create faces
        Ft=[Fq(:,1) Fq(:,2) indVm(:);...
            Fq(:,2) Fq(:,3) indVm(:);...
            Fq(:,3) Fq(:,4) indVm(:);...
            Fq(:,4) Fq(:,1) indVm(:)];
        Ct=repmat(Cq,[4,1]);
end
    
varargout{1}=Ft;
varargout{2}=Vt;
varargout{3}=Ct;

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
