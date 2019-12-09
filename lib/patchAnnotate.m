function [varargout]=patchAnnotate(F,V,d,varargin)

% function [hv,hf]=patchAnnotate(F,V,d,varargin)
% ------------------------------------------------------------------------
%
% Change log: 
% 2019/08/06 Created
% 
% ------------------------------------------------------------------------

%% Parse input

% Force 3D coordinates
if size(V,2)==2
    V(:,3)=0; 
end

%% Get node indices 

nodeIndices=1:1:size(V,1);
faceIndices=1:1:size(F,1);

%% Get coordinates for text segments

if isempty(d) 
    d=mean(patchEdgeLengths(F,V))/10; %
end

[NF,VF,NV]=patchNormal(F,V);
VV=V+NV.*d;
VF=VF+NF.*d;

%% Annotate nodes and faces

hv=pointAnnotate(VV,nodeIndices,varargin{:},'FontWeight','Normal','FontAngle','italic'); %Annotate Nodes
hf=pointAnnotate(VF,faceIndices,varargin{:},'FontWeight','Bold'); %Annotate faces

%% Collect output

varargout{1}=hv;
varargout{2}=hf;

end
