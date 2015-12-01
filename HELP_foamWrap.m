%% geoSphere
% Below is a demonstration of the features of the |geoSphere| function

%clear; close all; clc;

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha=1;
edgeColor=0.1*ones(1,3);
edgeWidth=1;

%Modifying standard figure properties
figStruct.Name='GIBBON'; %Figure name
figStruct.Color='w'; %Figure background color

%Custom figure properties
figStruct.ColorDef='white'; %Setting colordefinitions to black
figStruct.ScreenOffset=0; %Setting sp

%Defining geodesic dome
r=1; %sphere radius

[F,V,~]=geoSphere(2,r);
% [F,V]=parasaurolophus;
% [F,V]=cow;
% [F,V]=graphicsModels(4);
% [F,V]=stanford_bunny;

% nc=6;
cmap=gjet(4);
% [C,logicConverged]=triSurfPermuteColor(F,V,nc);

%%
cFigure(figStruct); 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(F,V,r/3);
patch('Faces',F,'Vertices',V,'FaceColor',[84 22 180]./256,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','k');

camlight headlight; 
colormap(cmap);
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on; axis off; 


%%

% [F,V,C,indIni]=triPolyDualRefine(F,V);

C=[1:size(F,1)]';

n=1; 
for q=1:1:n
    [F,V]=subtri(F,V,1);
    C=repmat(C,[4 1]);
end

Z=V(:,3);
ZF=mean(Z(F),2);
logicTop=ZF>0; 
% C(logicTop)=0; 


%%
cFigure(figStruct); 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor','none');

camlight headlight; lighting phong; 
colormap(cmap);
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on; axis off; 

%%

cPar.n=3; 
cPar.dirFlip=1; 
% cPar.foamThickness=0.02;
cParSmooth.Method='HC';
cParSmooth.n=25;
cPar.Smooth=cParSmooth; 

%%
L_remove=true(size(F,1),1);
[FT,VT,CT]=foamWrap(F,V,C,cPar,logicTop);

%%
cFigure(figStruct);  hold on; 
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

% hp=patch('Faces',FT(CT==1,:),'Vertices',VT,'FaceColor','flat','CData',CT(CT==1,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hp=patch('Faces',FT(CT==1,:),'Vertices',VT,'FaceColor','flat','CData',CT(CT==1),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
hp=patch('Faces',FT(CT==2,:),'Vertices',VT,'FaceColor','flat','CData',CT(CT==2),'FaceAlpha',1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% hp=patch('Faces',FT(CT==1,:),'Vertices',VT,'FaceColor','w','FaceAlpha',faceAlpha,'lineWidth',4,'edgeColor',edgeColor);
% hp=patch('Faces',FT(CT==2,:),'Vertices',VT,'FaceColor',0.75.*[180 64 16]./255,'FaceAlpha',faceAlpha,'lineWidth',4,'edgeColor',edgeColor);
% plotV(VT(indRigid,:),'r.','MarkerSize',35);

camlight headlight; 
% lighting phong;
colormap((gray(2))); 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
axis off;


%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>