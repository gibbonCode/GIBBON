function [hp]=patchNormPlot(varargin)

% function [hp]=patchNormPlot(F,V,a,pathType)
%
% Simple function to plot surface normal vectors with magnitude a for the
% patch data specified by the faces F and vertices V. The option third
% input patchType sets the type of surface normal to use, i.e. the face
% normal if patchType='f' or the vertex normal if patchType='v'; 
%
% Quiver patch options can be altered using the output handle "hp".
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 29/06/2013
%
% Change log: 
% 2014/06/02 added vertex normal code
% ------------------------------------------------------------------------
%%

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        a=varargin{3};
        patchType='f'; 
    case 4
        F=varargin{1};
        V=varargin{2};
        a=varargin{3};
        patchType=varargin{4};         
    otherwise
        error('Wrong numer of input arguments!');
end

%Get face normals
switch patchType
    case 'f'       
        [N,Vn]=patchNormal(F,V);
    case 'v'
        [~,~,N]=patchNormal(F,V);
        Vn=V;
    otherwise
        error('Wrong patchType specified!');
end

%Derive quiver patch data
[Fni,Vni,~]=quiver3Dpatch(Vn(:,1),Vn(:,2),Vn(:,3),N(:,1),N(:,2),N(:,3),[],[a a]);

%Patch plotting the face normals
hp=patch('Faces',Fni,'Vertices',Vni);

%Set defaults
set(hp,'EdgeColor','none','FaceColor','k');
