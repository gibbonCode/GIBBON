%% Bardiya Akhbari (04.03.2020)
%% Nataliya Perevoshchikova (04.04.2020)

clc;clear;close all;
fontSize=30;
%% Save Path
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
SavePath=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data_carpal_model','motion');
load_path= fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data_carpal_model');

saveOn=1;
Start=1;
Step=10;

%There are motion data from 5 specimens:
subjNum = 02751; info_txt = [num2str(subjNum,'%05.f') '_F_Y_'];
% subjNum = 02772; info_txt = [num2str(subjNum,'%05.f') '_M_Y_'];
% subjNum = 50029; info_txt = [num2str(subjNum,'%05.f') '_F_O_'];
%subjNum = 50090; info_txt = [num2str(subjNum,'%05.f') '_M_O_'];
%subjNum = 01424; info_txt = [num2str(subjNum,'%05.f') '_M_Y_'];

iBone = 1;
for bone_list = {'cap','sca','lun'}
    bone_str = bone_list{1};
    % myRT = (forViz_RT) * cur_RT{iB,1} * transposefX4(forViz_RT);
    myPath     = fullfile(load_path,'ACS_Models',[info_txt bone_str '_ics.wrl']);
    bone_models{iBone,1} = GB_read_vrml_fast_BA(myPath);
    
    %motion_RT_16x1{iBone,1} = csvread(fullfile(load_path,'Motion_Pose',[info_txt 'RT_' bone_str '_ba_model_16x1.tra']));
    motion_RT_16x1{iBone,1} = csvread(fullfile(load_path,'Motion_Pose_Sinusoidal',[info_txt 'RT_' bone_str '_ba_model_16x1.tra']));
    neut_16x1{iBone,1} = reshape(csvread(fullfile(load_path,'Neutral_Pose',[info_txt bone_str '_neut_pose.tra'])),4,4)';
    %bone_6x1_TraRot{iBone,1} = csvread(fullfile(load_path,'Motion_Pose',[info_txt 'EUL_' bone_str '_ba_model_eul_6x1.tra']));
    bone_6x1_TraRot{iBone,1} = csvread(fullfile(load_path,'Motion_Pose_Sinusoidal',[info_txt 'EUL_' bone_str '_ba_model_eul_6x1.tra']));
    
    iBone = iBone + 1;
end

RAD_Model = GB_read_vrml_fast_BA(fullfile(load_path,'ACS_Models',[info_txt 'rad_acs.wrl']));
nFrames = length(bone_6x1_TraRot{1, 1});

%% Organize Data
for iFrame = Start:Step:nFrames
    for iBone = [1,2,3] % 2: scaphoid, 3: lunate
        motion_RT_4x4{iFrame,iBone} = reshape(motion_RT_16x1{iBone,1}(iFrame,:),[4 4])';
    end
end

% Show Models
figure(1); clf;
set(gcf,'Position',[364 250 560 420],'Color',[1 1 1]);
for iFrame = Start:Step:nFrames
    figure(1); clf;
    GB_plotboneVRML_New(RAD_Model,eye(4,4),'k',0,gca);
    GB_plotboneVRML_New(bone_models{1},motion_RT_4x4{iFrame,1}*neut_16x1{1,1},'b',0,gca);
    GB_plotboneVRML_New(bone_models{2},motion_RT_4x4{iFrame,2}*neut_16x1{2,1},'b',0,gca);
    GB_plotboneVRML_New(bone_models{3},motion_RT_4x4{iFrame,3}*neut_16x1{3,1},'b',0,gca);

    GB_giveMeSomeLights(1);
    view(-90,-90);
    axis equal;
    xlim([-40,30])
    ylim([-30,30]);
    set(gca,'XColor','none','YColor','none','ZColor','none')
    pause(0.2);
end

%% Make it like Autoscoper Code
% bone_list = {'cap','sca','lun'}
for iFrame = Start:Step:nFrames
    SCA_4x4_POSE{iFrame,1} =  motion_RT_4x4{iFrame,2}*neut_16x1{2,1};
    LUN_4x4_POSE{iFrame,1} =  motion_RT_4x4{iFrame,3}*neut_16x1{3,1};
    CAP_4x4_POSE{iFrame,1} =  motion_RT_4x4{iFrame,1}*neut_16x1{1,1};
end

ref_frame = Start;
SCA_4x4_REF =  SCA_4x4_POSE{ref_frame,1};
LUN_4x4_REF =  LUN_4x4_POSE{ref_frame,1};
CAP_4x4_REF =  CAP_4x4_POSE{ref_frame,1};

count=1;
for iFrame = Start:Step:nFrames
    SCA_4x4_MOTION{iFrame,1} = SCA_4x4_POSE{iFrame,1} * GB_transposefX4(SCA_4x4_REF);
    LUN_4x4_MOTION{iFrame,1} = LUN_4x4_POSE{iFrame,1} * GB_transposefX4(LUN_4x4_REF);
    CAP_4x4_MOTION{iFrame,1} = CAP_4x4_POSE{iFrame,1} * GB_transposefX4(CAP_4x4_REF);

    Trans_L=LUN_4x4_MOTION{iFrame,1};
    Trans_S=SCA_4x4_MOTION{iFrame,1};   
    Trans_C=CAP_4x4_MOTION{iFrame,1};
    
    TFL=any(isnan(LUN_4x4_MOTION{iFrame,1}));
    TFS=any(isnan(SCA_4x4_MOTION{iFrame,1}));
    TFC=any(isnan(CAP_4x4_MOTION{iFrame,1}));
    
    
    if any(TFL==1) || any(TFS==1) || any(TFC==1)
        disp('Do something motion matrix has NaN');
        disp(iFrame);
    else
        save(fullfile(SavePath,sprintf('Flex_Ext_%d_Carpal_Model_Lunate.mat',count)), 'Trans_L', '-mat');
        save(fullfile(SavePath,sprintf('Flex_Ext_%d_Carpal_Model_Scaphoid.mat',count)), 'Trans_S', '-mat');
        save(fullfile(SavePath,sprintf('Flex_Ext_%d_Carpal_Model_Capitate.mat',count)), 'Trans_C', '-mat');
        count=count+1;
    end
    
    
    SCA_6x1_TraRot(iFrame,:) = GB_RT_to_eul_4x4(SCA_4x4_MOTION{iFrame,1});
    LUN_6x1_TraRot(iFrame,:) = GB_RT_to_eul_4x4(LUN_4x4_MOTION{iFrame,1});
    CAP_6x1_TraRot(iFrame,:) = GB_RT_to_eul_4x4(CAP_4x4_MOTION{iFrame,1});
end

FC=bone_models{1}.conn(:,1:3)+1;%faces
FL=bone_models{3}.conn(:,1:3)+1;%faces
FS=bone_models{2}.conn(:,1:3)+1;%faces


LUN_MAT = LUN_4x4_REF;
SCA_MAT = SCA_4x4_REF;
CAP_MAT = CAP_4x4_REF;

VL=tform(LUN_MAT,bone_models{3}.pts);
VS=tform(SCA_MAT,bone_models{2}.pts);
VC=tform(CAP_MAT,bone_models{1}.pts);

VR=RAD_Model.pts; % vertices scaphoid
FR=RAD_Model.conn(:,1:3)+1;%faces

%%
% Initiate a visualization in a figure window while storing an object
cFigure;
gpatch(FL,VL,'gw','none',1);
gpatch(FS,VS,'bw','none',1);
gpatch(FR,VR,'w','none',0.8);
gpatch(FC,VC,'y','none',0.8);

camlight headlight;
view(-90,-90);
axis equal;
axisGeom(gca,fontSize);
xlim([-40,30])
ylim([-30,30]);
set(gca,'XColor','none','YColor','none','ZColor','none');
drawnow; 

%% Save output file of first frame
if saveOn==1
    dataStructOut.VL=VL;
    dataStructOut.FL=FL;
    dataStructOut.VS=VS;
    dataStructOut.FS=FS;
    dataStructOut.VR=VR;
    dataStructOut.FR=FR;
    dataStructOut.VC=VC;
    dataStructOut.FC=FC;
    
    saveName=fullfile(SavePath,['modelGeometry_',num2str(subjNum,'%05.f'),'.mat']);
    save(saveName,'-struct','dataStructOut');
end

%% Show Motion time dependent
time_tp = 0:seconds(1/200):seconds((nFrames-1)/200);

CAP_6x1_TraRot = bone_6x1_TraRot{1,1};
SCA_6x1_TraRot = bone_6x1_TraRot{2,1};
LUN_6x1_TraRot = bone_6x1_TraRot{3,1};

cFigure;
subplot(2,3,1); hold on;
plot(time_tp,CAP_6x1_TraRot(:,4),'k','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,5),'k:','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,6),'k--','LineWidth',2);

ylabel('Rotation ({\circ})');

ylim([-70,70]);

set(gca,'FontSize',24); 
grid on; box on;


subplot(2,3,2); hold on;
plot(time_tp,SCA_6x1_TraRot(:,4),'k','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,5),'k:','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,6),'k--','LineWidth',2);

ylim([-70,70]);
set(gca,'FontSize',24); grid on; box on;


subplot(2,3,3); hold on;
hl(1)=plot(time_tp,LUN_6x1_TraRot(:,4),'k','LineWidth',2);
hl(2)=plot(time_tp,LUN_6x1_TraRot(:,5),'k:','LineWidth',2);
hl(3)=plot(time_tp,LUN_6x1_TraRot(:,6),'k--','LineWidth',2);

ylim([-70,70]);

legend(hl,{'Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation',},'Position',[0.52 0.92 0.15 0.0869]);
legend('Orientation','horizontal');
legend('boxoff');

set(gca,'FontSize',24); 
grid on; box on;


subplot(2,3,4); hold on;
plot(time_tp,CAP_6x1_TraRot(:,1),'k','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,2),'k:','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,3),'k--','LineWidth',2);

ylabel('Translation (mm)')

ylim([-15,13]);

set(gca,'FontSize',24); 
grid on; box on;
subplot(2,3,5); hold on;
plot(time_tp,SCA_6x1_TraRot(:,1),'k','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,2),'k:','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,3),'k--','LineWidth',2);


ylim([-15,13]);
set(gca,'FontSize',24); grid on; box on;

subplot(2,3,6); hold on;
plot(time_tp,LUN_6x1_TraRot(:,1),'k','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot(:,2),'k:','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot(:,3),'k--','LineWidth',2);
ylim([-15,13]);

set(gca,'FontSize',24); grid on; box on;

set(gcf,'Color',[1 1 1],'Position',[223 373 1970 682]);

drawnow; 

