function [varargout]=patch2tri(varargin)

% function [Ft,Vt,Ct]=patch2tri(F,V,C)
% ------------------------------------------------------------------------
% This function converts patch data to triangular patch data. It is assumed
% that the centroid (mean of patch vertices) lies on the patch and
% connection of the patch edges to the centroid yields appropriate
% triangles. 
% 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/12/22 Created based on quad2tri
%------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};        
        C=[];
    case 3        
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};    
    otherwise
        error('False number of input arguments');
end

%%

V_new=zeros(size(F,1),size(V,2));
for q=1:1:size(V,2);
    X=V(:,q);
    V_new(:,q)=mean(X(F),2);    
end

columnCopy=repmat((1:size(F,2))',1,2)'; columnCopy=[columnCopy(2:end) 1];
Ft=F(:,columnCopy);
Ft=reshape(Ft',2,numel(Ft)/2)';
indAdd=size(V,1)+(1:size(V_new,1));
Vt=[V;V_new];
indAddRep=repmat(indAdd(:)',size(F,2),1);
Ft(:,3)=indAddRep(:);

%Color vector now assumes all schemes allow for simply copying of the color
%information
if ~isempty(C)
    Ct=C(Ft(:,3)-size(V,1));
else
    Ct=[];
end

varargout{1}=Ft;
varargout{2}=Vt;
varargout{3}=Ct;

 
%% <-- GIBBON footer text --> 
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
