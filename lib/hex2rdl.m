function [varargout]=hex2rdl(varargin)

% function [Ep,Et,Vs]=hex2rdl(E,V,inputStruct)
% ------------------------------------------------------------------------
% This function can generate a rhombicdodecahedron lattice or a diamond
% lattice, based on the input hexahedral element mesh specified by E, (the
% element array), and V (the vertex array). 
%
% The optional input structure containts the following defaults: 

% inputStructDefault.shrinkFactor=0.25; 
% inputStructDefault.latticePhaseType=1; % 1 = "bubble" centred, 2 = vertex centred, 3 = nested
% inputStructDefault.latticeType=1; % rhombic-dodecahedron (1) or Diamond (2) 
% inputStructDefault.removeNonConnected=1; %Removing non-connected tetrahedra

% ------------------------------------------------------------------------
%%

%The default control parameters
inputStructDefault.shrinkFactor=0.25;
inputStructDefault.latticePhaseType=1; % 1 = "bubble" centred, 2 = vertex centred, 3 = nested
inputStructDefault.latticeType=1; % rhombic-dodecahedron (1) or Diamond (2) 
inputStructDefault.removeNonConnected=1; %Removing non-connected tetrahedra

%Get inputs
E=varargin{1};
V=varargin{2};
switch nargin
    case 2
        inputStruct=inputStructDefault;
    case 3
        inputStruct=varargin{3};
end

%Complement input structure with default
inputStruct=structComplete(inputStruct,inputStructDefault,1);


%Access control parameters
latticePhaseType=inputStruct.latticePhaseType;

switch latticePhaseType
    case 3 %Nested
        inputStruct.latticePhaseType=1;
        [Ep1,Et1,Vs1]=hex2rdl(E,V,inputStruct);

        inputStruct.latticePhaseType=2;
        [Ep2,Et2,Vs2]=hex2rdl(E,V,inputStruct);

        Ep=[Ep1; Ep2+size(Vs1,1)];
        Et=[Et1; Et2+size(Vs1,1)];
        Vs=[Vs1;Vs2];

    otherwise

        shrinkFactor=inputStruct.shrinkFactor;
        latticeType=inputStruct.latticeType;
        removeNonConnected=inputStruct.removeNonConnected;

        %% subdevide hex

        C=repmat((1:8),size(E,1),1);
        C=C(:);

        [E,V]=subHex(E,V,1);

        %%

        switch latticePhaseType
            case 1 % "Bubble" centered
                logicFlip=ismember(C,[1 3 6 8]);
                E(logicFlip,:)=E(logicFlip,[4 3 7 8 1 2 6 5 ]);
            case 2 % Vertex centered
                logicFlip=ismember(C,[2 4 5 7]);
                E(logicFlip,:)=E(logicFlip,[4 3 7 8 1 2 6 5 ]);
        end

        tetInd=[5 1 8 6; 7 3 6 8; 2 1 6 3; 4 1 3 8; 6 8 3 1];
        a=tetInd';
        a=a(:)';
        A=E(:,a);
        Et=reshape(A',4,5.*size(E,1))';

        logicFlip_Et=repmat(logicFlip',5,1);
        logicFlip_Et=logicFlip_Et(:);
        [Ftp,logicFlip_Ftp]=element2patch(Et,logicFlip_Et,'tet4');
        Cs=repmat((1:1:size(Et,1))',4,1);

        Ftp_sort=sort(Ftp,2);
        [~,~,faceIdUni]=unique(Ftp_sort,'rows');
        faceId=(1:1:size(Ftp,1))';

        logic1 = rem(Cs,5)==0;
        faceIdUni1=faceIdUni(logic1);

        logic2 = ~logic1 & ismember(faceIdUni,faceIdUni1);
        faceIdUni2=faceIdUni(logic2);

        [~,indSort1]=sort(faceIdUni1);
        [~,indSort2]=sort(faceIdUni2);

        faceId1=faceId(logic1);
        faceId1=faceId1(indSort1);
        faceId2=faceId(logic2);
        faceId2=faceId2(indSort2);

        logicFlip_Ep=logicFlip_Ftp;
        logicFlip_Ep=logicFlip_Ep(logic1);
        logicFlip_Ep=logicFlip_Ep(indSort1);

        %%

        V_scale=V(Et(:,1),:);
        Et_ind=(1:1:size(Et,1))';
        logic1_Es=rem(Et_ind,5)==0;
        for q=1:1:size(V,2)
            X=V(:,q);
            if nnz(logic1_Es)==1
                V_scale(logic1_Es,q)=mean(X(Et(logic1_Es,:)),1);
            else
                V_scale(logic1_Es,q)=mean(X(Et(logic1_Es,:)),2);
            end
        end
        [Et,Vs]=scalePatch(Et,V,shrinkFactor,V_scale);

        %%

        switch latticeType
            case 1 %Rhombic dodecahedron
                logicPickTet=true(size(logic1_Es));
                logicPickPenta=true(size(logicFlip_Ep));
            case 2 %Diamond lattice
                logicPickTet= ~logic1_Es | (logicFlip_Et==0 & logic1_Es);
                logicPickPenta=logicFlip_Ep==0;
        end
        
        Fs=element2patch(Et,logic1_Es,'tet4'); %  [Fs,logic1_Es_Fs]=element2patch(Et,logic1_Es,'tet4');
        Ep=[Fs(faceId1(logicPickPenta),:) fliplr(Fs(faceId2(logicPickPenta),:))];
        Et=Et(logicPickTet,:);

        % Merge nodes
        [~,Vs,~,indFix]=mergeVertices(Fs,Vs);        
        Et=indFix(Et);
        Ep=indFix(Ep);              

        if latticeType==2 && removeNonConnected==1            
            Et=Et(any(ismember(Et,Ep),2),:);
        end

        %Remove unused nodes
        indUsed=unique([Et(:);Ep(:)]);
        [~,Vs,indFix]=patchCleanUnused(indUsed,Vs);
        Et=indFix(Et);
        Ep=indFix(Ep);
        
end

%% Collect output
varargout{1}=Ep;
varargout{2}=Et;
varargout{3}=Vs;
% varargout{4}=logicFlip_Ep(logicPickPenta);
% varargout{5}=logic1_Es(logicPickTet,:);

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
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
