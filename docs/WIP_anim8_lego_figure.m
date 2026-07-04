clear; close all; clc;

%%

% Path names
gibbonFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(gibbonFolder,'data','STL');
savePath=fullfile(gibbonFolder,'data','OBJ');
fileName=fullfile(savePath,'lego_figure.obj');

fileNames={'lego_figure_head.stl',...        %1
    'lego_figure_eyes_mouth.stl',...  %2
    'lego_figure_torso.stl',...       %3
    'lego_figure_arm_right.stl'...    %4
    'lego_figure_arm_left.stl'...     %5
    'lego_figure_hand_right.stl'...   %6
    'lego_figure_hand_left.stl'...    %7
    'lego_figure_belt.stl',...        %8
    'lego_figure_leg_right.stl',...   %9
    'lego_figure_leg_left.stl',...    %10
    };

cMap =[  0         0         0;... %1 Eyes/mouth
	     0         0    0.3500;... %2 Belt
         0    0.6000    0.1451;... %3 Arms
         0.2000    0.4118    0.9098;... %4 Legs
         0.8353    0.0588    0.1451;... %5 Torso
         1.0000    1.0000         0];   %6 Head/hands

% Origins and rotation parameters and vectors
rotAngleJoints=0.25*pi;
translationTotal=50;
numStepsAnimateStride=10;

bodyTranslation=[25 0 0];
bodyRotation=euler2DCM([0 0 0.1*pi]);

rotVecArmRight=[0.0000    -0.9828   0.1848];
origin_arm_right=[0.1427   -6.2127    5.7342];

rotVecArmLeft=rotVecArmRight;
rotVecArmLeft(2)=-rotVecArmRight(2);
origin_arm_left=origin_arm_right;
origin_arm_left(2)=-origin_arm_left(2);

rotVecLegRight=[0 -1 0];
origin_leg_right=[0 -1 -8.2];

rotVecLegLeft=[0 1 0];
origin_leg_left=[0 1 -8.2];

%% Import STL files

numSets=numel(fileNames);
F=cell(1,numSets);
V=cell(1,numSets);
for q=1:1:numSets
    TR = stlread(fullfile(pathName,fileNames{q}));
    
    Fn=TR.ConnectivityList; %Faces
    Vn=TR.Points; %Vertices
    [Fn,Vn]=mergeVertices(Fn,Vn); % Merging nodes
    F{q}=Fn;
    V{q}=Vn;
    C{q}=q*ones(size(Fn,1),1); %Color label
end

%Join model sets
[Fm,Vm,C_file]=joinElementSets(F,V,C);

% File numbers          1  2  3  4  5  6  7  8  9  10
fileColormappingVector=[6  1  5  3  3  6  6  2  4  4]'; %Mapping vector
Cm=fileColormappingVector(C_file); %Face colors for colormap driven colors

%%

cFigure; hold on;
gpatch(Fm,Vm,Cm,'none',1);
axisGeom; view(45,30); camlight('headlight');
colormap(cMap); icolorbar;
drawnow;

%% 
%Visualize origins/rotation vectors

cFigure; hold on;
gpatch(Fm,Vm,'w','none',0.1);

plotV(origin_arm_right,'r.','markerSize',25)
quiverVec(origin_arm_right,rotVecArmRight,5,'r');

plotV(origin_arm_left,'b.','markerSize',25)
quiverVec(origin_arm_left,rotVecArmLeft,5,'b');

plotV(origin_arm_left,'b.','markerSize',25)
quiverVec(origin_arm_left,rotVecArmLeft,5,'b');

plotV(origin_leg_right,'r.','markerSize',25)
quiverVec(origin_leg_right,rotVecLegRight,5,'r');

plotV(origin_leg_left,'b.','markerSize',25)
quiverVec(origin_leg_left,rotVecLegLeft,5,'b');

axisGeom; view(45,30); camlight('headlight');
drawnow;

%% Create animation

hf=cFigure; hold on;
gpatch(Fm,Vm,'w','none',0.25);
hp=gpatch(Fm,Vm,Cm,'none',1);
axisGeom; view(45,30); camlight('headlight');
colormap(cMap); %icolorbar;
drawnow;

a1=linspace(0,rotAngleJoints,numStepsAnimateStride);
a=[a1 fliplr(a1(1:end-1)) -a1(2:end) fliplr(-a1(1:end-1)) ];

nSteps=numel(a); %Number of animation steps
t=linspace(0,translationTotal,nSteps);
animStruct.Time=linspace(0,1,nSteps);

V1=V;

for qs=1:1:nSteps
    V1=V; %Copy original coordinates
    
    rotAngle=a(qs); %The current angle

    %Set up motion structure
    motionStruct.body.T=bodyTranslation;
    motionStruct.body.Orientation=bodyRotation;
    motionStruct.body.Origin=[0 0 0];
    
    motionStruct.arm_right.RotAxis=rotVecArmRight;
    motionStruct.arm_right.RotAngle=-rotAngle;
    motionStruct.arm_right.Origin=origin_arm_right;
    
    motionStruct.arm_left.RotAxis=rotVecArmLeft;
    motionStruct.arm_left.RotAngle=-rotAngle;
    motionStruct.arm_left.Origin=origin_arm_left;
    
    motionStruct.leg_right.RotAxis=rotVecLegRight;
    motionStruct.leg_right.RotAngle=rotAngle;
    motionStruct.leg_right.Origin=origin_leg_right;
    
    motionStruct.leg_left.RotAxis=rotVecLegLeft;
    motionStruct.leg_left.RotAngle=rotAngle;
    motionStruct.leg_left.Origin=origin_leg_left;
    
    %Rotate right arm+hand
    for q=[4 6]
        Vn=V{q};
        Vn=Vn-motionStruct.arm_right.Origin(ones(size(Vn,1),1),:);
        R=vecAngle2Rot(motionStruct.arm_right.RotAngle,motionStruct.arm_right.RotAxis);
        Vn=Vn*R;
        Vn=Vn+motionStruct.arm_right.Origin(ones(size(Vn,1),1),:);
        V1{q}=Vn;
    end
    
    %Rotate left arm+hand
    for q=[5 7]
        Vn=V{q};
        Vn=Vn-motionStruct.arm_left.Origin(ones(size(Vn,1),1),:);
        R=vecAngle2Rot(motionStruct.arm_left.RotAngle,motionStruct.arm_left.RotAxis);
        Vn=Vn*R;
        Vn=Vn+motionStruct.arm_left.Origin(ones(size(Vn,1),1),:);
        V1{q}=Vn;
    end
    
    %Rotate right leg
    for q=9
        Vn=V{q};
        Vn=Vn-motionStruct.leg_right.Origin(ones(size(Vn,1),1),:);
        R=vecAngle2Rot(motionStruct.leg_right.RotAngle,motionStruct.leg_right.RotAxis);
        Vn=Vn*R;
        Vn=Vn+motionStruct.leg_right.Origin(ones(size(Vn,1),1),:);
        V1{q}=Vn;
    end
    
    %Rotate left leg
    for q=10
        Vn=V{q};
        Vn=Vn-motionStruct.leg_left.Origin(ones(size(Vn,1),1),:);
        R=vecAngle2Rot(motionStruct.leg_left.RotAngle,motionStruct.leg_left.RotAxis);
        Vn=Vn*R;
        Vn=Vn+motionStruct.leg_left.Origin(ones(size(Vn,1),1),:);
        V1{q}=Vn;
    end
    
    %Translate all
    motionStruct.body.T=[t(qs) 0 0];
    
    for q=1:1:numSets
        Vn=V1{q};
        Vn=Vn+motionStruct.body.T(ones(size(Vn,1),1),:);
        V1{q}=Vn;
    end
    
    %Join coordinate sets for updating model
    [~,Vmp,~]=joinElementSets(F,V1,C);
    
    % Populate animation structure           
    animStruct.Handles{qs}=hp; %Handles of objects to animate
    animStruct.Props{qs}={'Vertices'}; %Properties of objects to animate
    animStruct.Set{qs}={Vmp}; %Property values for to set in order to animate
end

anim8(hf,animStruct);
axis manual; %axis off;



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
