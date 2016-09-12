%% import_FEB
% Below is a demonstration of the features of the |import_FEB| function

%%
close all; clc; clear;

%%
% Plot settings
figColor='w'; figColorDef='white';
fontSize=10;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;

%% Importing a FEB file

%Set main folder
defaultFolder = fileparts(mfilename('fullpath'));
pathName=fullfile(defaultFolder,'data','FEB');

testCase=2; 
switch testCase
    case 1
        febFileNamePart='tetGenModel.feb'; %febio_spec 1.2
    case 2
        febFileNamePart='tempModel_2p0.feb'; %febio_spec 2.0
end
febFileName=fullfile(pathName,febFileNamePart);
[febXML,nodeStruct,elementCell]=import_FEB(febFileName);

%%
% Content:

nodeStruct
elementCell{:}

V=nodeStruct.N;

%% Plotting model

% Plotting the example model surfaces
hf1=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('Visualizing fullmodel','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

for q=1:1:numel(elementCell)

    E=elementCell{q}.E;
    C=elementCell{q}.E_mat;
    if numel(C)==1
        C=C.*ones(size(E,1),1);
    end
        
    [F,C]=element2patch(E,C); 
    
    subplot(1,2,1);
    
    title('Full model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    hp=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','flat','Cdata',C,'FaceAlpha',1);
    view(3); axis tight;  axis equal;  grid on;
    colormap(autumn);   
    set(gca,'FontSize',fontSize);
    
    subplot(1,2,2);
    
    %Selecting half of the model to see interior
    Y=V(:,2); YE=mean(Y(E),2);
    L=YE>mean(Y);
    [Fs,Cs]=element2patch(E(L,:),C(L));
    
    title('Cut view of model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    hps=patch('Faces',Fs,'Vertices',V,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    view(3); axis tight;  axis equal;  grid on;
    colormap(autumn);
    set(gca,'FontSize',fontSize);
    
    drawnow;

end

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
