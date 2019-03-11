function [varargout]=triEdgeSplit(varargin)

% function [Fs,Vs,Cs,CVs,CFs]=triEdgeSplit(F,V,E,CF,CV)
%-------------------------------------------------------------------------
%
%
% 
%-------------------------------------------------------------------------

%% Parse input 
% 

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        CF=[];
        CV=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        CF=varargin{4};
        CV=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        CF=varargin{4};
        CV=varargin{5};
end

if isempty(CF)
    CF=ones(size(F,1),1);
end

if isempty(CV)
    CV=ones(size(V,1),1);
end

%%

logicFacesNow=sum(ismember(F,E),2)==2; %Logic for faces to split

F_keep=F(~logicFacesNow,:); %Faces to keep
CF_keep=CF(~logicFacesNow,:); %Colors to keep

F_remove=F(logicFacesNow,:); %Faces to split
CF_remove=CF(logicFacesNow,:); %Colors of faces to split

Ns=patchNormal(F_remove,V); %Normals of faces to split

%Create and append mid edge point to vertex list
Vs=[V; mean(V(E,:),1)];
CVs=[CV; mean(CV(E,:),1)];
indNewPoint=size(Vs,1); %Index of new point

%Create new faces
F_new=[];
CF_new=[];
for q=1:1:size(F_remove,1)
    F_now=F_remove(q,:); %The current face for splitting
    N_now=Ns(q,:); %Current correct face normal;
    C_now=CF_remove(q,:); %The current face for splitting
    indNotEdge=F_now(~ismember(F_now,E)); %Index of "other point" (not member of edge)
    f_new=[indNotEdge E(1)  indNewPoint; indNotEdge indNewPoint E(2)]; %New faces
    N_new=patchNormal(f_new,Vs); %Current normals of new faces    
    logicFlip=dot(N_new,N_now(ones(size(N_new,1),1),:),2)<0; %Logic for faces to flip
    
    f_new(logicFlip,:)=fliplr(f_new(logicFlip,:)); %Flip faces to correct normals    
    F_new=[F_new; f_new]; %Append new faces
    
    c_new=[C_now(ones(size(f_new,1),1),:)];
    CF_new=[CF_new; c_new]; %Append new colors
end
Fs=[F_keep;F_new]; %Compose full face set
CFs=[CF_keep;CF_new]; %Compose full face set
Cs=[ones(size(F_keep,1),1); 2*ones(size(F_new,1),1);]; %Colors for labelling

%% Collect output

varargout{1}=Fs;
varargout{2}=Vs;
switch nargout
    case 3
        varargout{1}=Fs;
        varargout{2}=Vs;        
        varargout{3}=Cs;        
    case 4
        varargout{1}=Fs;
        varargout{2}=Vs;        
        varargout{3}=Cs;
        varargout{4}=CFs;
    case 5
        varargout{1}=Fs;
        varargout{2}=Vs;        
        varargout{3}=Cs;
        varargout{4}=CFs;
        varargout{5}=CVs;
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
