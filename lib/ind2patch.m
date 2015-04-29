function [F,V,C]=ind2patch(IND,M,ptype)

% function [F,V,C]=ind2patch(IND,M,ptype)
% ------------------------------------------------------------------------
%
% This function generates patch data (faces 'F', vertices 'V' and color
% data 'C') for 3D images. The patches are only generated for the voxels
% specified by the linear indices in 'IND'. The variable 'ptype' indicates
% the type of patch:
%
% 'v'               Voxel patch data with unshared vertices and faces
%                  such that each voxel has 8 unshared vertices and 6
%                  unshared faces (may be faster than 'vu' and 'vb'
%                  options which require UNIQUE costly computations).
% 'vu'             Voxel patch data such that where possible voxels share
%                  vertices and faces (making patch data computation slower
%                  but plotting more memory efficient).
%                  FaceColor data is averaged for shared faces.
% 'vb'             Voxel patch data faces are exported for unique unshared
%                  faces e.g. only boundary for enclosed volume (plotted
%                  data is visually equivalent to 'v' and 'vu' options when
%                  FaceAlpha is 1)
% 'si', 'sj', 'sk'    Mid-voxel slice patch data for i, j and k direction
%                  respectively
% 'siu', 'sju', 'sku' Same as 'si', 'sj', 'sk' but with double points
%                  removed.
% 'h'              Creates a hexahedral element description instead (e.g
%                  nx8) element data.
% 'hu'             Same as 'h' but with shared unique nodes.
%
%%% EXAMPLE
% clear all; close all; clc; 
% 
% %% Simulating 3D image
% [X,Y,Z]=meshgrid(linspace(-4.77,4.77,25));
% phi=(1+sqrt(5))/2;
% M=2 - (cos(X + phi*Y) + cos(X - phi*Y) + cos(Y + phi*Z) + cos(Y - phi*Z) + cos(Z - phi*X) + cos(Z + phi*X));
% M=M./max(M(:)); %Normalise, not required
% 
% figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
% hold on; xlabel('X-J','FontSize',20);ylabel('Y-I','FontSize',20);zlabel('Z-K','FontSize',20);
% 
% %% Creating and plotting patch data
% 
% %Setting up indices for X direction slices
% S=round(size(M,2)./2); %Selection of middle slice
% L_plot=false(size(M)); L_plot(:,S,:)=1;
% IND=find(L_plot);
% [F,V,C]=ind2patch(IND,M,'sj'); %Creating patch data for x mid-voxel slices
% hs=patch('Faces',F,'Vertices',V,'EdgeColor','none', 'CData',C,'FaceColor','flat','FaceAlpha',0.75);
% 
% %Setting up indices for Y direction slices
% S=round(size(M,1)./2); %Selection of middle slice
% L_plot=false(size(M)); L_plot(S,:,:)=1;
% IND=find(L_plot);
% [F,V,C]=ind2patch(IND,M,'si'); %Creating patch data for y mid-voxel slices
% hs=patch('Faces',F,'Vertices',V,'EdgeColor','none', 'CData',C,'FaceColor','flat','FaceAlpha',0.75);
% 
% %Setting up indices for Z direction slices
% S=round(size(M,3)./2); %Selection of middle slice
% L_plot=false(size(M)); L_plot(:,:,S)=1;
% IND=find(L_plot);
% [F,V,C]=ind2patch(IND,M,'sk'); %Creating patch data for z mid-voxel slices
% hs=patch('Faces',F,'Vertices',V,'EdgeColor','none', 'CData',C,'FaceColor','flat','FaceAlpha',0.75);
% 
% %Setting up indices for voxels 
% IND=find(M>-0.2 & M<=0); % A selection of low intensity voxels
% [F,V,C]=ind2patch(IND,M,'v'); %Creating patch data for selection of low voxels
% hs=patch('Faces',F,'Vertices',V,'EdgeColor','k', 'CData',C,'FaceColor','flat','FaceAlpha',1);
% 
% %Setting up indices for voxels to plot
% IND=find(M>0.9); % A selection of high intensity voxels. 
% [F,V,C]=ind2patch(IND,M,'v'); %Creating patch data for selection of high voxels
% hs=patch('Faces',F,'Vertices',V,'EdgeColor','k', 'CData',C,'FaceColor','flat','FaceAlpha',1);
% 
% axis equal; view(3); axis tight; colormap jet; colorbar; caxis([0 1]); grid on;
% set(gca,'FontSize',20);
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

%% PARSING INPUT

%Get indices to patch
if islogical(IND) %treated as a logic index and converted to linear indices
    IND=find(IND); 
end
IND=IND(:);

%Deal with complex data
if ~isreal(M)
    M=real(M); 
    disp('Warning: Complex data, imaginary components were ignored')
end

%%

numIND=numel(IND);

[I,J,K] = ind2sub(size(M),IND); %Convert to I,J,K coordinates

switch ptype
    case {'si','siu'} %I midvoxel slice
        i_shift=ones(size(IND))*[ 0    0    0    0   ];
        j_shift=ones(size(IND))*[-0.5 -0.5  0.5  0.5 ];
        k_shift=ones(size(IND))*[-0.5  0.5  0.5 -0.5 ];
        faceOrder=[1 2 3 4];
        numNodes=4;
    case {'sj','sju'} %J midvoxel slice
        i_shift=ones(size(IND))*[-0.5  0.5  0.5 -0.5 ];
        j_shift=ones(size(IND))*[ 0    0    0    0   ];
        k_shift=ones(size(IND))*[-0.5 -0.5  0.5  0.5 ];
        faceOrder=[1 2 3 4];
        numNodes=4;
    case {'sk','sku'} %K midvoxel slice
        i_shift=ones(size(IND))*[-0.5  0.5  0.5 -0.5 ];
        j_shift=ones(size(IND))*[-0.5 -0.5  0.5  0.5 ];
        k_shift=ones(size(IND))*[ 0    0    0    0   ];
        faceOrder=[1 2 3 4];
        numNodes=4;
    case {'v','vu','vb'} %Voxels, same with unique faces and vertices, only boundary (unshared) faces and vertices
        i_shift=ones(size(IND))*[-0.5  0.5  0.5 -0.5   -0.5  0.5  0.5 -0.5];
        j_shift=ones(size(IND))*[-0.5 -0.5  0.5  0.5   -0.5 -0.5  0.5  0.5];
        k_shift=ones(size(IND))*[ 0.5  0.5  0.5  0.5   -0.5 -0.5 -0.5 -0.5];
        faceOrder=[1 2 3 4;...    %Top
            5 6 7 8;...    %Bottom
            1 2 6 5; ...   %Left face
            3 4 8 7; ...   %Right face
            1 4 8 5; ...   %Back
            2 3 7 6];      %Front
        numNodes=8;
    case {'h','hu'} %Hexahedral element
        i_shift=ones(size(IND))*[-0.5  0.5  0.5 -0.5   -0.5  0.5  0.5 -0.5];
        j_shift=ones(size(IND))*[-0.5 -0.5  0.5  0.5   -0.5 -0.5  0.5  0.5];
        k_shift=ones(size(IND))*[ 0.5  0.5  0.5  0.5   -0.5 -0.5 -0.5 -0.5];
        faceOrder=[1 2 3 4 5 6 7 8];
        numNodes=8;
    otherwise
        warning('wrong input for argument ptype, valid inputs are s for surfaces patches and v for voxel patches');
end
numFacesPerVoxel=size(faceOrder,1);
numNodesPerFace=size(faceOrder,2);

VI=I*ones(1,numNodes);
VJ=J*ones(1,numNodes);
VK=K*ones(1,numNodes);

VI=VI+i_shift; VI=VI'; VI=VI(:);
VJ=VJ+j_shift; VJ=VJ'; VJ=VJ(:);
VK=VK+k_shift; VK=VK'; VK=VK(:);

V=zeros(length(VI),3);

V=[VJ VI VK]; %N.B. I and J direction are switched

%Creates faces
Fi=ones(numIND,numNodesPerFace);
F=repmat(Fi,numFacesPerVoxel,1);
b=(numNodes.*((1:1:numIND)'-1))*ones(1,numNodesPerFace);
for q=1:1:numFacesPerVoxel
    Fi=ones(size(IND))*faceOrder(q,:)+b;   
    F((1+numIND*(q-1)):numIND*q,:)=Fi;
end

%Preparing color information
C=M(IND);
C=C(:);

if ~isa(C,'double'); %If the image is not a double, then convert
    C=double(C);
end
if numFacesPerVoxel>1;
    C=repmat(C,[numFacesPerVoxel,1]);
end

switch ptype
    case {'vu','vb'}
        
        %Removing double VERTICES
        [V,~,IND_IND]=unique(V,'rows'); %works well because coordinates are integers +/- 0.5
        F=IND_IND(F);
        
        %Removing double FACES
        Fs=sort(F,2); %Sort so faces with same nodes have the same rows
        [~,IND_F,IND_F_2]=unique(Fs,'rows'); %integer unique operation
        F=F(IND_F,:);
        
        %Get face counts (used for averaging colour)
        numF=size(Fs,1); numFuni=size(F,1);
        logicColourMatrixEntry=sparse(IND_F_2,1:numF,1,numFuni,numF,numF);
        F_count=full(sum(logicColourMatrixEntry,2));
        
        %Fixing face colors, shared faces now obtain mean colour
        if ~isempty(C)
            sharedColourMatrixSparse=sparse(IND_F_2,1:numF,C,numFuni,numF,numF);
            C=full(sum(sharedColourMatrixSparse,2))./F_count;
        end
        
        %Removing non-unique faces for option 'vb'
        if strcmp(ptype,'vb')
            %Only keeping un-shared faces and color data
            F=F(F_count==1,:);
            C=C(F_count==1,:);
            
            %Remove excess points
            indVUni=unique(F(:));                        
            V=V(indVUni,:);            
                        
            %Now fix face matrix entries
            indVFix=nan(size(V,1),1);
            indVFix(indVUni)=1:numel(indVUni);            
            F=indVFix(F);             
        end
    case  {'siu','sju','sku'}
        %Removing double vertices
        [V,~,IND_IND]=unique(V,'rows'); %works well because coordinates are integers +/- 0.5
        F=IND_IND(F);
    case 'hu'
        %Removing double vertices
        [V,~,IND_IND]=unique(V,'rows'); %works well because coordinates are integers +/- 0.5
        F=IND_IND(F);         
end

%Tranposing F if required (possible if numel(find(IND))=1)
if size(F,2)==1 
    F=F';
end

end






