function [varargout]=triPolyDualRefine(varargin)

% function [Ft,Vt,Ct,indIni]=triPolyDualRefine(F,V)
% ------------------------------------------------------------------------
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/05/01 Updated for GIBBON
% 2015/05/01 Added varargin/varargout
%
% TO DO allow multiple iterations with "book keeping" of color/index data
%------------------------------------------------------------------------

%% Parse input

% switch nargin
%     case 2
        F=varargin{1}; %Faces
        V=varargin{2}; %Vertices
%         CF=ones(size(F,1),1); %Color data
%     case 3        
%         F=varargin{1}; %Faces
%         V=varargin{2}; %Vertices
%         CF=varargin{3}; %Color data
% end

%%

[Vd,~,Fds]=patch_dual(V,F);

Cd=(1:1:size(Fds,1))';
Cds=Cd(:,ones(size(Fds,2),1));

Fds_t=Fds'; 
Cds_t=Cds';

Vt=[Vd;V];

Q1=Fds_t(:);
L1=Q1>0; 
Q1=Q1(L1);
Ct=Cds_t(:);
Ct=Ct(L1);

Q2=Fds(:,[2:size(Fds,2) 1])';
Q2=Q2(:);
L2=Q2>0; 
Q2=Q2(L2);

q=1:size(Fds,1);
Q3=q(ones(size(Fds,2),1),:)+size(Vd,1);
Q3=Q3(:);
A=Q3;
Q3=Q3(L1);
Ft=[Q1 Q2 Q3];

indIni=(size(Vd,1)+1):(size(Vd,1)+size(V,1));

%% Collect output

varargout{1}=Ft;
varargout{2}=Vt;
varargout{3}=Ct;
varargout{4}=indIni;



 
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
