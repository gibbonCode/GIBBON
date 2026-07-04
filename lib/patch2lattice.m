function [Fn,Vn,Cn]=patch2lattice(varargin)

% function [Fn,Vn,Cn]=patch2lattice(E,V,cPar)

% 2017/04/26 fixed bug for older MATLAB versions in relation to subtracting
% columns from matrices e.g. subtracting an nx1 column from all columns in
% an nxm matrix.
% 2018/05/07 Updated handling of default input (no change in functionality)
% 2018/05/07 Added default boundary face use if not provided at all. 

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        cPar=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        cPar=varargin{3};
end

%The default control parameters
cParDefault.shrinkFactor=0.25;
cParDefault.latticeSide=1;
cParDefault.meshType='tri';

%Complement input structure with default
[cPar]=structComplete(cPar,cParDefault,0);

numNode=size(F,2);
switch numNode
    case 3
        inputMeshType='tri';
    case 4
        inputMeshType='quad';
    otherwise 
        inputMeshType='other';
end

%%
switch cPar.latticeSide
    case 1 
        s=1-cPar.shrinkFactor; 
    case 2
        s=cPar.shrinkFactor; 
end
[Fc,Vc]=scalePatch(F,V,s); %Get shrunk elements
E=patchEdges(F,0); %Input face edges
Ec=patchEdges(Fc,0); %Shrunk face edges

switch cPar.latticeSide
    case 1                                
        switch cPar.meshType
            case 'tri'                
                Ec=Ec+size(V,1);
                Fn=[E           Ec(:,2);...
                    Ec(:,[2 1]) E(:,1);...
                    ]; %New face set
                Vn=[V;Vc]; %Combined node set
                Cn=[];
            case 'quad'
                Fn=[E fliplr(Ec)+size(V,1)]; %New face set
                Vn=[V;Vc]; %Combined node set
                Cn=[zeros(size(E,1),1)];
        end        
    case 2
        [Ecc,Vcc]=scalePatch(E,V,s); %Get contracted edge coordinates
        
        switch cPar.meshType
            case 'tri'
                switch inputMeshType
                    case 'tri'
                        Ecc=Ecc+size(Vc,1);
                        Fn=[Fc;...
                            Ec(:,2)  Ec(:,1) Ecc(:,1);...
                            Ecc(:,1) Ecc(:,2) Ec(:,2)]; %New face set
                        Vn=[Vc;Vcc]; %Combined node set
                        Cn=[zeros(size(Fc,1),1); ones(size(Ec,1),1); ones(size(Ecc,1),1);];
                    case 'quad'
                        [Fct,~]=quad2tri(Fc,Vc,'f');
                        Ecc=Ecc+size(Vc,1);
                        Fn=[Fct;...
                            Ec(:,2)  Ec(:,1) Ecc(:,1);...
                            Ecc(:,1) Ecc(:,2) Ec(:,2)]; %New face set
                        Vn=[Vc;Vcc]; %Combined node set
                        Cn=[zeros(size(Fct,1),1); ones(size(Ec,1),1); ones(size(Ecc,1),1);];
                    otherwise
                        Vm=patchCentre(Fc,Vc);
                        Fcm=[Ec (size(Vc,1)+size(Vcc,1)+1)+repmat((1:size(Fc,1))',size(Fc,2),1) ];
                        
                        Ecc=Ecc+size(Vc,1);
                        Fn=[Fcm;...
                            Ec(:,2)  Ec(:,1) Ecc(:,1);...
                            Ecc(:,1) Ecc(:,2) Ec(:,2)]; %New face set
                        Vn=[Vc;Vcc;Vm]; %Combined node set
                        Cn=[];
                end
            case 'quad'
                switch inputMeshType
                    case 'tri'
                        Fq=[fliplr(Ec) Ecc+size(Vc,1)]; %New face set
                        Vq=[Vc;Vcc]; %Combined node set
                        [Fqq,Vqq]=subQuad(Fq,Vq,1);
                        [Fqt,Vqt]=tri2quad(Fc,Vc,1);
                        Fn=[Fqq;Fqt+size(Vqq,1)];
                        Vn=[Vqq;Vqt];
                        Cn=[zeros(size(Fqq,1),1); ones(size(Fqt,1),1); ];
                    case 'quad'                        
                        Fn=[Fc; fliplr(Ec) Ecc+size(Vc,1)]; %New face set
                        Vn=[Vc; Vcc]; %Combined node set                        
                        Cn=[];
                    otherwise
%                         Vm=patchCentre(Fc,Vc);
                        Fn=[fliplr(Ec) Ecc+size(Vc,1)]; %New face set
                        Vn=[Vc; Vcc]; %Combined node set
                        Cn=[];
                end                      
        end
end

%Removing double vertices
[Fn,Vn]=mergeVertices(Fn,Vn);

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
