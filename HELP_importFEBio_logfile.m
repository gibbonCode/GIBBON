%% import_FEB
% Below is a demonstration of the features of the |import_FEB| function

%%
clear all; close all; clc;

%%
% Plot settings
figColor='w'; figColorDef='white';
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=0.5;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=50;


%% Importing a FEB file

%Set main folder
defaultFolder = fileparts(mfilename('fullpath'));
pathName=fullfile(defaultFolder,'data','FEB');

febFileNamePart='tetGenModel.feb';
febFileName=fullfile(pathName,febFileNamePart);
[febXML,nodeStruct,elementCell]=import_FEB(febFileName);

%%
%%
% Content:

nodeStruct
elementCell{1}

V=nodeStruct.N;

%% Plotting model

% Plotting the example model surfaces
hf1=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('Visualizing fullmodel','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

uniqueMaterialIndices=[];
for q=1:1:numel(elementCell)
    uniqueMaterialIndices=unique([uniqueMaterialIndices(:); elementCell{q}.E_mat(:)]);
    E=elementCell{q}.E;
    switch elementCell{q}.E_type
        case {'tri3', 'quad4'}
            F=E;
            V=nodeStruct.N;
            C=elementCell{q}.E_mat;
        case {'hex8', 'tet4'}
            E=elementCell{q}.E;
            [F,C]=element2patch(E,elementCell{q}.E_mat); %Creates faces and colors (e.g. stress) for patch based plotting
    end
    
    hp=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','flat','Cdata',C,'FaceAlpha',0.8);
    
    subplot(1,2,1);
    
    title('Full model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    hp=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','flat','Cdata',C,'FaceAlpha',1);
    view(3); axis tight;  axis equal;  grid on;
    colormap(autumn);
    camlight headlight;
    set(gca,'FontSize',fontSize);
    
    subplot(1,2,2);
    
    %Selecting half of the model to see interior
    Y=V(:,3); ZE=mean(Y(E),2);
    L=ZE<mean(Y);
    [Fs,Cs]=element2patch(E(L,:),C(L));
    
    title('Cut view of model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    view(3); axis tight;  axis equal;  grid on;
    colormap(autumn);
    camlight headlight;
    set(gca,'FontSize',fontSize);
    drawnow;
    
end

%% IMPORTING LOG FILE
% Importing nodal displacements from a log file

FEB_outputName=fullfile(pathName,'tetGenModel_node_out.txt');
[~, N_disp_mat,~]=importFEBio_logfile(FEB_outputName); %Nodal displacements

DN=N_disp_mat(:,2:end,end); %Final nodal displacements

%% EXAMPLE CREATING NODE SET IN DEFORMED STATE
V_def=V+DN;

%%
% Plotting the model

%Selecting half of the model to see interior
Z=V(:,3); ZE=mean(Z(E),2);
L=ZE<mean(Z);
[Fs,~]=element2patch(E(L,:),[]);

Cs=sqrt(sum(DN.^2,2)); %Color towards displacement magnitude

hf1=figuremax(figColor,figColorDef);
title('Cut view of deformed model from imported results','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',Fs,'Vertices',V_def,'FaceColor','flat','FaceVertexCData',Cs);

view(3); axis tight;  axis equal;  grid on;
colormap jet; colorbar; shading interp;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;


%%
%
% <<gibbVerySmall.gif>>
%
% GIBBON
%
% Kevin M. Moerman (kevinmoerman@hotmail.com)