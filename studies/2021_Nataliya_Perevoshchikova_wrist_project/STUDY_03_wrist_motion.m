%% Importing and Visualizing Autoscoper Data
%% Bardiya Akhbari (04.03.2020)
clear; clc; close all;
fontSize=30;
%% Generic Path based on folder structure for Griffith College 
my_root_path = fileparts(pwd);
find_reference_trial = 99;  % Keep this 99 to use reference frame, unless put 1 for ignoring reference
frame_ref            = 1;   % Use this frame as reference


%% Save Path
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
SavePath=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','motion');
%SavePath_FILT=fullfile(defaultFolder,'MATLAB','MTP','Data','Brown_University','FILT');
%% Tracked Files
tracked_folder = fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','Tracked_Files');
tracked_dir_all = dir(tracked_folder);

saveOn=1;

nModels = length(tracked_dir_all)-2;
for iModel = 1:nModels
    tracked_path = fullfile(tracked_folder,tracked_dir_all(iModel+2).name);
    
    tracked_16x1{iModel,1} = csvread(tracked_path); % Read file (autoscoper defaults saved values)
    tracked_16x1{iModel,1} = tracked_16x1{iModel,1}(:,1:16); % Make sure data is 16 by 1
  
    nFrames = size(tracked_16x1{iModel,1},1);

    if nFrames > 20 % Smoothing if enough frame numbers
        [tracked_16x1{iModel,1},tracked_16x1_filt_er{iModel,1}] = GB_smooth_aut_kin(tracked_16x1{iModel,1});
    end
end

%% Models in VRML format
models_folder = fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','ACS_Models');
models_dir_all = dir(models_folder);
for iModel = 1:nModels
    models_acs{iModel,1} = GB_read_vrml_fast_BA(fullfile(models_folder,models_dir_all(iModel+2).name));
end

%% Transformation from Autoscoper CS to Anatomical CS
tform_folder = fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','Cropped2ACS_TFM');

tform_dir_all = dir(tform_folder);
for iModel = 1:nModels
    tform_aut2acs{iModel,1}    = GB_readGeoTFM(fullfile(tform_folder,tform_dir_all(iModel+2).name));
end

%% Read Kinematics
% convert to 4x4 matrix and ACS
for iFrame = 1:nFrames
    for iModel = 1:nModels
        XCS_tracked_4x4{iModel,1}{iFrame,1} = reshape(tracked_16x1{iModel,1}(iFrame,:),[4 4])';
        
        ACS_tracked_4x4{iModel,1}{iFrame,1} = XCS_tracked_4x4{iModel,1}{iFrame,1} * ...
            GB_transposefX4(tform_aut2acs{iModel,1});
    end
end

%% Visualize the Model in X-ray Space
bone_col = [227,218,201]./255;
Start=1000;
Step=5;%2

for iFrame = Start:Step:nFrames
%for iFrame = Start
    figure(1); clf; hold on;
    for iModel = [1,2,3,4]
        GB_plotboneVRML_New(models_acs{iModel,1},...
            GB_transposefX4(ACS_tracked_4x4{3,1}{iFrame,1}) * ACS_tracked_4x4{iModel,1}{iFrame,1},...
            bone_col,0,gca,0.4);  
    end    
    GB_giveMeSomeLights(1);
    grid on; axis equal;
    title(['X-Ray Space, Frame #: ' num2str(iFrame)]);
    set(gcf,'Color',[1 1 1],'Position',[1200 200 1037 800])
    xlim([-50,20])
    set(gca,'Color',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],...
        'CameraPositionMode','manual','CameraPosition',[-25 -250 500],'View',[0 60],...
        'CameraUpVector',[-1 0 0],'CameraUpVectorMode','manual')  
    
end

FL = models_acs{2,1}.conn(:,1:3)+1;
FS = models_acs{4,1}.conn(:,1:3)+1;
FC= models_acs{1,1}.conn(:,1:3)+1;
FR= models_acs{3,1}.conn(:,1:3)+1;

LUN_MAT=GB_transposefX4(ACS_tracked_4x4{3,1}{Start,1}) * ACS_tracked_4x4{2,1}{Start,1};
SCA_MAT=GB_transposefX4(ACS_tracked_4x4{3,1}{Start,1}) * ACS_tracked_4x4{4,1}{Start,1};
CAP_MAT=GB_transposefX4(ACS_tracked_4x4{3,1}{Start,1}) * ACS_tracked_4x4{1,1}{Start,1};
RAD_MAT=GB_transposefX4(ACS_tracked_4x4{3,1}{Start,1}) * ACS_tracked_4x4{3,1}{Start,1};

VL=tform(LUN_MAT,models_acs{2,1}.pts);
VS=tform(SCA_MAT,models_acs{4,1}.pts);
VC=tform(CAP_MAT,models_acs{1,1}.pts);
VR=tform(RAD_MAT,models_acs{3,1}.pts);


%% Refrence Pose
%sel_model = 2; % 1: is Capitate, 2: is Lunate, 3: is Radius, 4: is Scaphoid (CHECK models_dir_all variable)

% Relative position/orientation of "iModel" bone for "iFrame" frame in
% Radius Coordinate System (_POSE)
for iFrame = Start:nFrames
    SCA_4x4_POSE{iFrame,1} =  GB_transposefX4(ACS_tracked_4x4{3,1}{iFrame,1}) * ACS_tracked_4x4{4,1}{iFrame,1};
    LUN_4x4_POSE{iFrame,1} =  GB_transposefX4(ACS_tracked_4x4{3,1}{iFrame,1}) * ACS_tracked_4x4{2,1}{iFrame,1};
    CAP_4x4_POSE{iFrame,1} =  GB_transposefX4(ACS_tracked_4x4{3,1}{iFrame,1}) * ACS_tracked_4x4{1,1}{iFrame,1};
    RAD_4x4_POSE{iFrame,1} =  GB_transposefX4(ACS_tracked_4x4{3,1}{iFrame,1}) * ACS_tracked_4x4{3,1}{iFrame,1};
end

% Reference Position (_REF)
ref_frame = Start;
SCA_4x4_REF =  SCA_4x4_POSE{ref_frame,1};
LUN_4x4_REF =  LUN_4x4_POSE{ref_frame,1};
CAP_4x4_REF =  CAP_4x4_POSE{ref_frame,1};
RAD_4x4_REF =  RAD_4x4_POSE{ref_frame,1};


% Calculate the movement relative to this reference position (_MOTION)
count=1;
for iFrame = Start:Step:nFrames
    SCA_4x4_MOTION{iFrame,1} =  SCA_4x4_POSE{iFrame,1} * GB_transposefX4(SCA_4x4_REF);
    LUN_4x4_MOTION{iFrame,1} =  LUN_4x4_POSE{iFrame,1} * GB_transposefX4(LUN_4x4_REF);
    CAP_4x4_MOTION{iFrame,1} =  CAP_4x4_POSE{iFrame,1} * GB_transposefX4(CAP_4x4_REF);
    RAD_4x4_MOTION{iFrame,1} =  RAD_4x4_POSE{iFrame,1} * GB_transposefX4(RAD_4x4_REF);
    
    Trans_L=LUN_4x4_MOTION{iFrame,1};
    Trans_S=SCA_4x4_MOTION{iFrame,1};   
    Trans_C=CAP_4x4_MOTION{iFrame,1};
    Trans_R=RAD_4x4_MOTION{iFrame,1};   
    
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Lunate.mat',count)), 'Trans_L', '-mat');
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Scaphoid.mat',count)), 'Trans_S', '-mat');
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Capitate.mat',count)), 'Trans_C', '-mat');
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Radius.mat',count)), 'Trans_R', '-mat');

    % Convert the 4x4 matrix to 6x1 matrix
    % 1: proximal-distal translation (mm)
    % 2: radial-ulnar translation (mm)
    % 3: volar-dorsal translation (mm)
    % 4: pronation-supination rotation (degree)
    % 5: flexion-extension rotation (degree)
    % 6: radial-ulnar deviation rotation (degree)
    
    SCA_6x1_TraRot(count,:) = GB_RT_to_eul_4x4(SCA_4x4_MOTION{iFrame,1});
    LUN_6x1_TraRot(count,:) = GB_RT_to_eul_4x4(LUN_4x4_MOTION{iFrame,1});
    CAP_6x1_TraRot(count,:) = GB_RT_to_eul_4x4(CAP_4x4_MOTION{iFrame,1});
    RAD_6x1_TraRot(count,:) = GB_RT_to_eul_4x4(RAD_4x4_MOTION{iFrame,1});
    
    count=count+1;
end

%% Showing the rotation of the bone relative to the reference position
time_tp = 0:seconds(1/200)*Step:seconds((nFrames-Start)/200);
        
figure(2); clf; 
subplot(2,3,1); hold on;
plot(time_tp,SCA_6x1_TraRot(:,4),'r','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,5),'g','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,6),'b','LineWidth',2);
ylabel('Scaphoid Rotation ({\circ})')
ylim([-70,70]);
set(gca,'FontSize',18); grid on; box on;

subplot(2,3,2); hold on;
plot(time_tp,LUN_6x1_TraRot(:,4),'r','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot(:,5),'g','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot(:,6),'b','LineWidth',2);
ylabel('Lunate Rotation ({\circ})')
ylim([-70,70]);
%legend('Flexion-Extension','Radial-Ulnar Deviation','location','southeast');
set(gca,'FontSize',18); grid on; box on;

subplot(2,3,3); hold on;
plot(time_tp,CAP_6x1_TraRot(:,4),'r','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,5),'g','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,6),'b','LineWidth',2);
ylabel('Capitate Rotation ({\circ})')
ylim([-70,70]);
legend('Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
%legend('Flexion-Extension','Radial-Ulnar Deviation','location','southeast');
set(gca,'FontSize',18); grid on; box on;


subplot(2,3,4); hold on;
plot(time_tp,SCA_6x1_TraRot(:,1),'r','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,2),'g','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot(:,3),'b','LineWidth',2);
ylabel('Scaphoid Translation (mm)')
ylim([-8,13]);
set(gca,'FontSize',18); grid on; box on;

subplot(2,3,5); hold on;
plot(time_tp,LUN_6x1_TraRot(:,1),'r','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot(:,2),'g','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot(:,3),'b','LineWidth',2);
ylabel('Lunate Translation (mm)')
ylim([-8,13]);
%legend('Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
set(gca,'FontSize',18); grid on; box on;

subplot(2,3,6); hold on;
plot(time_tp,CAP_6x1_TraRot(:,1),'r','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,2),'g','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot(:,3),'b','LineWidth',2);
ylabel('Capitate Translation (mm)')
ylim([-8,13]);
legend('Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
%legend('Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
set(gca,'FontSize',18); grid on; box on;

set(gcf,'Color',[1 1 1],'Position',[223 373 1970 682]);
drawnow; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second round of filteration of translation and rotation
count=1;
for iFrame = Start:Step:nFrames  
    
    % Seperate translations X Y Z
    LUN_translations(count,:) =  LUN_4x4_MOTION{iFrame,1}(:,end);
    SCA_translations(count,:) =  SCA_4x4_MOTION{iFrame,1}(:,end);
    CAP_translations(count,:) =  CAP_4x4_MOTION{iFrame,1}(:,end);
    RAD_translations(count,:) =  RAD_4x4_MOTION{iFrame,1}(:,end);
    
    %  Seperate rotations and convert to quaternion
    LUN_rotations(count,:) =  quaternion(LUN_4x4_MOTION{iFrame,1}(1:3,1:3), 'rotmat', 'frame');
    SCA_rotations(count,:) =  quaternion(SCA_4x4_MOTION{iFrame,1}(1:3,1:3), 'rotmat', 'frame');
    CAP_rotations(count,:) =  quaternion(CAP_4x4_MOTION{iFrame,1}(1:3,1:3), 'rotmat', 'frame');
    RAD_rotations(count,:) =  quaternion(RAD_4x4_MOTION{iFrame,1}(1:3,1:3), 'rotmat', 'frame');
    count=count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = linspace(0.1,1,60);% cahnge the cutoff frequence from 0.1 to 1
%If fcut is set too high, less signal distortion
%occurs, but too much noise is allowed to pass. Conversely, if fc is too low, the
%noise is reduced drastically, but at the expense of increased signal distortion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 10;%sampling frequency
%fs = sampling frequency; fc= cutoff frequency of filter; order = filter order; lowORhigh = type of filter.
count=1;
LUN_4x4_FILT=eye(4,4);
SCA_4x4_FILT=eye(4,4);
CAP_4x4_FILT=eye(4,4);
RAD_4x4_FILT=eye(4,4);
for ifc=fc
    
    %Translation
    LUN_translations_filt= matfilt(fs, ifc, 2, LUN_translations ,'low');
    SCA_translations_filt= matfilt(fs, ifc, 2, SCA_translations ,'low');
    CAP_translations_filt= matfilt(fs, ifc, 2, CAP_translations ,'low');
    RAD_translations_filt= matfilt(fs, ifc, 2, RAD_translations ,'low');
   
    %Residual
    R_LUN(count,:)=sqrt(sum((LUN_translations-LUN_translations_filt).^2)/numel(Start:Step:nFrames));
    R_SCA(count,:)=sqrt(sum((SCA_translations-SCA_translations_filt).^2)/numel(Start:Step:nFrames));
    R_CAP(count,:)=sqrt(sum((CAP_translations-CAP_translations_filt).^2)/numel(Start:Step:nFrames));
    R_RAD(count,:)=sqrt(sum((RAD_translations-RAD_translations_filt).^2)/numel(Start:Step:nFrames)); 
    
    %Rotation
    [FilteredData_LUN] = Filter_Slerp(LUN_rotations,ifc);
    [FilteredData_SCA] = Filter_Slerp(SCA_rotations,ifc);
    [FilteredData_CAP] = Filter_Slerp(CAP_rotations,ifc);
    [FilteredData_RAD] = Filter_Slerp(RAD_rotations,ifc);
   
    %Convert quaternions to rotations matrices
    LUN_rotations_filt=rotmat(FilteredData_LUN, 'frame');
    SCA_rotations_filt=rotmat(FilteredData_SCA, 'frame');
    CAP_rotations_filt=rotmat(FilteredData_CAP, 'frame');
    RAD_rotations_filt=rotmat(FilteredData_RAD, 'frame');
    
   
    for i=1:1:length(LUN_rotations_filt)
        LUN_4x4_FILT(1:3,1:3)=LUN_rotations_filt(:,:,i);
        LUN_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(LUN_4x4_FILT);
        
        SCA_4x4_FILT(1:3,1:3)=SCA_rotations_filt(:,:,i);
        SCA_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(SCA_4x4_FILT);
        
        CAP_4x4_FILT(1:3,1:3)=CAP_rotations_filt(:,:,i);
        CAP_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(CAP_4x4_FILT);
        
        RAD_4x4_FILT(1:3,1:3)=RAD_rotations_filt(:,:,i);
        RAD_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(RAD_4x4_FILT);
    end
    
  
    %Residual
    Rrot_LUN(count,:)=sqrt(sum((LUN_6x1_TraRot(:,4:6)-LUN_6x1_TraRot_FILT(:,4:6)).^2)/numel(Start:Step:nFrames));
    Rrot_SCA(count,:)=sqrt(sum((SCA_6x1_TraRot(:,4:6)-SCA_6x1_TraRot_FILT(:,4:6)).^2)/numel(Start:Step:nFrames));
    Rrot_CAP(count,:)=sqrt(sum((CAP_6x1_TraRot(:,4:6)-CAP_6x1_TraRot_FILT(:,4:6)).^2)/numel(Start:Step:nFrames));
    Rrot_RAD(count,:)=sqrt(sum((RAD_6x1_TraRot(:,4:6)-RAD_6x1_TraRot_FILT(:,4:6)).^2)/numel(Start:Step:nFrames));
    
    count=count+1;
end

%Tangent line at the point to fitted curve: translation.
%Lunate
n=8;
p=polyfit(fc',R_LUN(:,1),n);
p1=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
%Calculate the gradient of the tangent at the A point
ch=5.5*n;
A=[fc(ch) R_LUN(ch,1)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
fc0=linspace(0,1,60);
Ytan=m*(fc0'-A(1))+A(2);
Yhor=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx1,fc_LUN]=polyxpoly(Yhor,fc0',p1,fc');

%Scaphoid
p=polyfit(fc',R_SCA(:,1),n);
p2=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
%Calculate the gradient of the tangent at the A point
ch=5.5*n;
A=[fc(ch) R_SCA(ch,1)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan1=m*(fc0'-A(1))+A(2);
Yhor1=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx2,fc_SCA]=polyxpoly(Yhor1,fc0',p2,fc');

%Capitate
p=polyfit(fc',R_CAP(:,1),n);
p2c=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
%Calculate the gradient of the tangent at the A point
ch=5.5*n;
A=[fc(ch) R_CAP(ch,1)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan1c=m*(fc0'-A(1))+A(2);
Yhor1c=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx2c,fc_CAP]=polyxpoly(Yhor1c,fc0',p2c,fc');

%Visualization of risidual analysis for translation 
figure(3); clf;
plot(fc',R_LUN(:,1),'ko');
hold on
plot(fc',R_SCA(:,1),'bo');
plot(fc',R_CAP(:,1),'go');

plot(fc',p1,'k-','LineWidth',2);
plot(fc',p2,'b-','LineWidth',2);
plot(fc',p2c,'g-','LineWidth',2);


plotV([fc_LUN xx1],'b*','MarkerSize',30);
plotV([fc_SCA xx2],'b*','MarkerSize',30);
plotV([fc_CAP xx2c],'b*','MarkerSize',30);

ylabel('Residual');
xlabel('fc');
legend('Translation lunate','Translation scaphoid','Translation capitate','location','northeast');
set(gca,'FontSize',18); grid on; box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tangent line at the point to fitted curve: rotation
%Lunate
n=8;
p=polyfit(fc',Rrot_LUN(:,1),n);
p3=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_LUN(ch,1)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan2=m*(fc0'-A(1))+A(2);
Yhor2=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx3,fc_rot1_LUN]=polyxpoly(Yhor2,fc0',p3,fc');

p=polyfit(fc',Rrot_LUN(:,2),n);
p4=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_LUN(ch,2)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan3=m*(fc0'-A(1))+A(2);
Yhor3=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx4,fc_rot2_LUN]=polyxpoly(Yhor3,fc0',p4,fc');

p=polyfit(fc',Rrot_LUN(:,3),n);
p5=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_LUN(ch,3)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan4=m*(fc0'-A(1))+A(2);
Yhor4=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx5,fc_rot3_LUN]=polyxpoly(Yhor4,fc0',p5,fc');

%Scaphoid
p=polyfit(fc',Rrot_SCA(:,1),n);
p6=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_SCA(ch,1)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan5=m*(fc0'-A(1))+A(2);
Yhor5=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx6,fc_rot1_SCA]=polyxpoly(Yhor5,fc0',p6,fc');

%
p=polyfit(fc',Rrot_SCA(:,2),n);
p7=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_SCA(ch,2)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan6=m*(fc0'-A(1))+A(2);
Yhor6=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx7,fc_rot2_SCA]=polyxpoly(Yhor6,fc0',p7,fc');

%
p=polyfit(fc',Rrot_SCA(:,3),n);
p8=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_SCA(ch,3)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan7=m*(fc0'-A(1))+A(2);
Yhor7=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx8,fc_rot3_SCA]=polyxpoly(Yhor7,fc0',p8,fc');

%
p=polyfit(fc',Rrot_CAP(:,1),n);
p9=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_CAP(ch,1)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan8=m*(fc0'-A(1))+A(2);
Yhor8=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx9,fc_rot1_CAP]=polyxpoly(Yhor8,fc0',p9,fc');

%
p=polyfit(fc',Rrot_CAP(:,2),n);
p10=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_CAP(ch,2)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan9=m*(fc0'-A(1))+A(2);
Yhor9=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx10,fc_rot2_CAP]=polyxpoly(Yhor9,fc0',p10,fc');

%
p=polyfit(fc',Rrot_CAP(:,3),n);
p11=p(1)*fc'.^n+p(2)*fc'.^(n-1)+p(3)*fc'.^(n-2)+p(4)*fc'.^(n-3)+p(5)*fc'.^(n-4)+p(6)*fc'.^(n-5)+p(7)*fc'.^(n-6)+p(8)*fc'.^(n-7)+p(9);
ch=7*n;
A=[fc(ch) Rrot_CAP(ch,3)];
m= n*p(1)*fc(ch).^(n-1)+...
    (n-1)*p(2)*fc(ch).^(n-2)+...
    (n-2)*p(3)*fc(ch).^(n-3)+...
    (n-3)*p(4)*fc(ch).^(n-4)+...
    (n-4)*p(5)*fc(ch).^(n-5)+...
    (n-5)*p(6)*fc(ch).^(n-6)+...
    (n-6)*p(7)*fc(ch).^(n-7)+...
    (n-7)*p(8)*fc(ch).^(n-8);
Ytan10=m*(fc0'-A(1))+A(2);
Yhor10=(-A(1)*m+A(2))*ones(size(fc0'));
%Find intersection point between horizontal line and polyfit
[xx11,fc_rot3_CAP]=polyxpoly(Yhor10,fc0',p11,fc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf;
plot(fc',Rrot_LUN(:,1),'ro');
hold on
plot(fc',Rrot_LUN(:,2),'go');
plot(fc',Rrot_LUN(:,3),'bo');

%po
plot(fc',p3,'r-','LineWidth',2);
plot(fc',p4,'g-','LineWidth',2);
plot(fc',p5,'b-','LineWidth',2);

plotV([fc_rot1_LUN xx3],'b*','MarkerSize',30);
plotV([fc_rot2_LUN xx4],'b*','MarkerSize',30);
plotV([fc_rot3_LUN xx5],'b*','MarkerSize',30);

ylabel('Residual Lunate');
xlabel('fc');
legend('Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
set(gca,'FontSize',18); grid on; box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); clf;
plot(fc',Rrot_SCA(:,1),'ro');
hold on
plot(fc',Rrot_SCA(:,2),'go');
plot(fc',Rrot_SCA(:,3),'bo');

%po
plot(fc',p6,'r-','LineWidth',2);
plot(fc',p7,'g-','LineWidth',2);
plot(fc',p8,'b-','LineWidth',2);


plotV([fc_rot1_SCA xx6],'b*','MarkerSize',30);
plotV([fc_rot2_SCA xx7],'b*','MarkerSize',30);
plotV([fc_rot3_SCA xx8],'b*','MarkerSize',30);

ylabel('Residual Scaphoid');
xlabel('fc');
legend('Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
set(gca,'FontSize',18); grid on; box on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf;
plot(fc',Rrot_CAP(:,1),'ro');
hold on
plot(fc',Rrot_CAP(:,2),'go');
plot(fc',Rrot_CAP(:,3),'bo');

%po
plot(fc',p9,'r-','LineWidth',2);
plot(fc',p10,'g-','LineWidth',2);
plot(fc',p11,'b-','LineWidth',2);

plot(fc0, Ytan8, 'r--','LineWidth',1);
 plot(fc0, Yhor9, 'g--','LineWidth',1);

plot(fc0, Ytan9, 'r--','LineWidth',1);
plot(fc0, Yhor9, 'g--','LineWidth',1)
% 
plot(fc0, Ytan10, 'r--','LineWidth',1);
 plot(fc0, Yhor10, 'g--','LineWidth',1)

plotV([fc_rot1_CAP xx9],'b*','MarkerSize',30);
plotV([fc_rot2_CAP xx10],'b*','MarkerSize',30);
plotV([fc_rot3_CAP xx11],'b*','MarkerSize',30);

ylabel('Residual Capitate');
xlabel('fc');
legend('Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation','location','northeast');
set(gca,'FontSize',18); grid on; box on;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Butterworth filter is applied based on optimised value of fc
%fc_LUN=0.6;fc_SCA=0.6;
% fc_rot2_LUN=0.25;fc_rot2_SCA=0.25;

LUN_translations_filt= matfilt(fs,fc_LUN, 2, LUN_translations ,'low');
SCA_translations_filt= matfilt(fs,fc_SCA, 2, SCA_translations ,'low');
CAP_translations_filt= matfilt(fs,fc_CAP, 2, CAP_translations ,'low');

%Slerp filter rotation based on optimized value of fc
fc_rot2_LUN=0.6;
fc_rot2_SCA=0.6;
fc_rot2_CAP=0.6;

[FilteredData] = Filter_Slerp(LUN_rotations,fc_rot2_LUN);
% Convert quaternions to rotations matrices
LUN_rotations_filt=rotmat(FilteredData, 'frame');

[FilteredData] = Filter_Slerp(SCA_rotations,fc_rot2_SCA);
% Convert quaternions to rotations matrices
SCA_rotations_filt=rotmat(FilteredData, 'frame');

[FilteredData] = Filter_Slerp(CAP_rotations,fc_rot2_CAP);
% Convert quaternions to rotations matrices
CAP_rotations_filt=rotmat(FilteredData, 'frame');

%% Compose matrices out of filtered translation and rotation
for i=1:1:length(LUN_rotations_filt)
    LUN_4x4_FILT(1:3,1:3)=LUN_rotations_filt(:,:,i);
    LUN_4x4_FILT(:,end)=LUN_translations_filt(i,:);
    LUN_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(LUN_4x4_FILT);
    
    SCA_4x4_FILT(1:3,1:3)=SCA_rotations_filt(:,:,i);
    SCA_4x4_FILT(:,end)=SCA_translations_filt(i,:);
    SCA_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(SCA_4x4_FILT);
    
    CAP_4x4_FILT(1:3,1:3)=CAP_rotations_filt(:,:,i);
    CAP_4x4_FILT(:,end)=CAP_translations_filt(i,:);
    CAP_6x1_TraRot_FILT(i,:) = GB_RT_to_eul_4x4(CAP_4x4_FILT);
    
    
    Trans_L=LUN_4x4_FILT;
    Trans_S=SCA_4x4_FILT;
    Trans_C=CAP_4x4_FILT;
  
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Lunate.mat',count)), 'Trans_L', '-mat');
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Scaphoid.mat',count)), 'Trans_S', '-mat');
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Capitate.mat',count)), 'Trans_C', '-mat');
    save(fullfile(SavePath,sprintf('Flex_Ext_%d_TriPlane_Radius.mat',count)), 'Trans_R', '-mat');
end
%% Showing the rotation of the bone relative to the reference position
time_tp = 0:seconds(1/200)*Step:seconds((length(LUN_rotations_filt)-1)/200*Step);
cFigure;
subplot(2,3,1); hold on;
plot(time_tp,CAP_6x1_TraRot_FILT(:,4),'k','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot_FILT(:,5),'k:','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot_FILT(:,6),'k--','LineWidth',2);

xline(seconds(0),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
xline(seconds(1.2),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',2);
xline(seconds(3.42),'-.','color',[0.9290, 0.6940, 0.1250],'LineWidth',2);
xline(seconds(3.67),'-.','color',[0.3010, 0.7450, 0.9330],'LineWidth',2);
xline(seconds(2.355),'-.','color',[0, 0.5, 0],'LineWidth',2);
xline(seconds(4.4850),'-.','color',[0, 0.5, 0],'LineWidth',2);

ylabel('Rotation ({\circ})');

ylim([-70,70]);
set(gca,'FontSize',24); 
grid on; box on;


subplot(2,3,2); hold on;
plot(time_tp,SCA_6x1_TraRot_FILT(:,4),'k','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot_FILT(:,5),'k:','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot_FILT(:,6),'k--','LineWidth',2);

ylim([-70,70]);
set(gca,'FontSize',24); grid on; box on;


subplot(2,3,3); hold on;
hl(1)=plot(time_tp,LUN_6x1_TraRot_FILT(:,4),'k','LineWidth',2);
hl(2)=plot(time_tp,LUN_6x1_TraRot_FILT(:,5),'k:','LineWidth',2);
hl(3)=plot(time_tp,LUN_6x1_TraRot_FILT(:,6),'k--','LineWidth',2);

ylim([-70,70]);

legend(hl,{'Pronation-Supination','Flexion-Extension','Radial-Ulnar Deviation',},'Position',[0.52 0.92 0.15 0.0869]);
legend('Orientation','horizontal');
legend('boxoff');

set(gca,'FontSize',24); 
grid on; box on;


subplot(2,3,4); hold on;
plot(time_tp,CAP_6x1_TraRot_FILT(:,1),'k','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot_FILT(:,2),'k:','LineWidth',2);
plot(time_tp,CAP_6x1_TraRot_FILT(:,3),'k--','LineWidth',2);

xline(seconds(0),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
xline(seconds(1.2),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',2);
xline(seconds(3.42),'-.','color',[0.9290, 0.6940, 0.1250],'LineWidth',2);
xline(seconds(3.67),'-.','color',[0.3010, 0.7450, 0.9330],'LineWidth',2);
xline(seconds(2.355),'-.','color',[0, 0.5, 0],'LineWidth',2);
xline(seconds(4.4850),'-.','color',[0, 0.5, 0],'LineWidth',2);

ylabel('Translation (mm)')

ylim([-8,13]);

set(gca,'FontSize',24); 
grid on; box on;
subplot(2,3,5); hold on;
plot(time_tp,SCA_6x1_TraRot_FILT(:,1),'k','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot_FILT(:,2),'k:','LineWidth',2);
plot(time_tp,SCA_6x1_TraRot_FILT(:,3),'k--','LineWidth',2);

ylim([-8,13]);
set(gca,'FontSize',24); grid on; box on;

subplot(2,3,6); hold on;
plot(time_tp,LUN_6x1_TraRot_FILT(:,1),'k','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot_FILT(:,2),'k:','LineWidth',2);
plot(time_tp,LUN_6x1_TraRot_FILT(:,3),'k--','LineWidth',2);

ylim([-8,13]);
set(gca,'FontSize',24); grid on; box on;

set(gcf,'Color',[1 1 1],'Position',[223 373 1970 682]);


