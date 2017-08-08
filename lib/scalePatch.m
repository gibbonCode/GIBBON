function [Fs,Vs,Cs]=scalePatch(F,V,C,scaleFactor)

% --------------------------------------------------------------------
% function [Fs,Vs,Cs]=scalePatch(F,V,C,scaleFactor)
%
% CHANGE LOG: 
% 19/12/2013 Fixed error related to single face entry, see if statement
% related to size(F,1)
%
% --------------------------------------------------------------------
%%

%If vertices are shared they need to be disconnected by adding vertices
% if numel(F)~=numel(V)
Vn=zeros(numel(F),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q); %The coordinate set
    XF=X(F); %The coordinates for each face vertex
    XFt=XF'; %transpose to prepare for column array
    Vn(:,q)=XFt(:); %Add column as part of new vertex matrix
end

%Fix indices in face matrix
indVn=1:numel(F);
Fs=reshape(indVn(:),size(F,2),size(F,1))';
Cs=C;

%Derive face means to shift vertices around mean
Vs=zeros(size(Vn));
for q=1:1:size(V,2)
    X=Vn(:,q);
    XF=X(Fs); 
    if size(F,1)==1
        XF=XF';
    end   
    meanXF=mean(XF,2)*ones(1,size(F,2));     
    meanXFt=meanXF'; 
    meanV(:,q)=meanXFt(:);   
end
Vn_shift=Vn-meanV; %Shift vertices
Vn_shift_scale=Vn_shift*scaleFactor; %Scale shifted vertices with respect to face center 
Vs=Vn_shift_scale+meanV; %Shift scaled vertices back
 
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
