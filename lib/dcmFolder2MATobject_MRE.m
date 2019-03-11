function dcmFolder2MATobject_MRE(PathName,MaxVarSize)

% function dcmFolder2MATobject_MRE(PathName,MaxVarSize)
% ------------------------------------------------------------------------
% This function converts the DICOM files in the folder PathName to a MAT
% object which is stored as a file called IMDAT in a new subfolder called
% IMDAT.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/03/16
%-------------------------------------------------------------------------

%%
%Get/set maximum variable size (sets save steps in case the variable size
%is exceeded)
if isempty(MaxVarSize)
    MaxVarSize=1e9; %Default
end

%% Getting DICOM file names and image dimension parameters
files = dir(fullfile(PathName,'*.dcm'));
files={files(1:end).name};
files=sort(files(:));
NumberOfFiles=numel(files);

%%

firstWarning_TriggerTime=0;
firstWarning_Rotation=0;
if NumberOfFiles==0
    error(['Folder ',PathName,' skipped since it does not contain .dcm files']);
end

%% SPECIFY Tags to collect

collectTags = { 'Basename', ...
    'PatientID', ...
    'StudyInstanceUID', ...
    'SeriesInstanceUID', ...
    'StudyID', ...
    'InstanceNumber', ...
    'SeriesNumber', ...
    'SeriesDescription', ...
    'StudyDescription', ...
    'ImageType', ...
    'Rows', ...
    'Columns', ...
    'PixelSpacing',...
    'ImagePositionPatient', ...
    'ImageOrientationPatient',...
    'SliceLocation', ...
    'SliceThickness', ...
    'SpacingBetweenSlices', ...
    'SliceGap', ...
    'NumberOfTemporalPositions', ...
    'TemporalPositionIdentifier', ...
    'HeartRate', ...
    'TriggerTime', ...
    'RepetitionTime', ...
    'EchoTime', ...
    'FlipAngle', ...
    'RescaleSlope', ...
    'ScaleSlope',...
    'RescaleIntercept',...
    'RescaleType',...
    'Height',...
    'Width',...
    'MRImagePhaseNumber'...
    'SliceNumberMR'...
    'NumberOfPCDirections',...
    'NumberOfPhasesMR',...
    'NumberOfSlicesMR',...
    'MRSeriesNrOfSlices',...
    'MRSeriesNrOfPhases'...
    'MRSeriesNrOfDynamicScans',...
    'MRImagePhaseNumber',...
    'ImagePlaneNumber'...
    'MRScaleSlope'...
    'MRScaleIntercept'...
    };
% 'EchoNumbers' is private and just a echotime identifier or sequence number
% 'SliceNumberMR' is also private (2001,100a) and can be usefull for identifing slices that are on same position
% N.B. Filename is added by matlab including the path, Basename is one we add ourself

%% Creating output file name and folder
fName='IMDAT';
foldername_out=fullfile(PathName,fName,filesep);
if ~exist(foldername_out,'dir') %create output folder if it does not exist already
    mkdir(foldername_out);
end
savename_out1=fullfile(foldername_out,'IMDAT.mat');
if exist(savename_out1,'file') %create output folder if it does not exist already
    delete(savename_out1); %Recreate file if it already exists
end
matObj = matfile(savename_out1,'Writable',true);

%% SETTING DICOM DICTIONARY
fName=fullfile(PathName,files{1});
dcmInfo_full=dicominfo(fName);
try
    if ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'philips'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        dicomdict('set','gibbon_dict.txt');
        disp('DICOM dictionary set to: gibbon_dict.txt');
        dictSetting=1;
    elseif ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'pms'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        dicomdict('set','gibbon_dict.txt');
        disp('DICOM dictionary set to: gibbon_dict.txt');
        dictSetting=1;
    elseif ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'siemens'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        dicomdict('set','dicom-dict-siemens.txt');
        disp('DICOM dictionary set to: dicom-dict-siemens.txt');
        dictSetting=2;
    else
        dicomdict('factory');
        warning('Unknown vendor, using DICOM dictionary factory settings');
        dictSetting=3;
    end
catch %e.g. if the manufacturer field is missing
    dicomdict('factory');
    warning('Unknown vendor, using DICOM dictionary factory settings');
    dictSetting=3;
end

%% LOADING DICOM INFO
hw = waitbar(0,'Loading DICOM info...');
for c=1:1:numel(files)
    fName=fullfile(PathName,files{c});
    %         dicomdict('get');
    dcmInfo_full=dicominfo(fName);
    dcmInfo_full.Basename=fName; %Add custom field
    iFields = find(isfield(dcmInfo_full,collectTags)); %indices of existing fields
    for iField=1:numel(iFields)
        fieldName = collectTags{iFields(iField)};
        dcmInfo(c).(fieldName) = dcmInfo_full.(fieldName); %allocated after first iteration
    end
    if c==1
        dcmInfo=repmat(dcmInfo,1,numel(files)); %Allocate after size of first is known
    end
    waitbar(c/numel(files),hw,['Loading DICIM info...',num2str(round(100.*c/numel(files))),'%']);
end
waitbar(c/numel(files),hw,['Saving DICOM info to MAT-file...',num2str(round(100.*c/numel(files))),'%']);
% matObj.dcmInfo=dcmInfo;
close(hw);

%Geometry information based on first info file
[G]=dicom3Dpar(dcmInfo(1));
matObj.G=G; %Saving geometry information

%TODO this part is not memory protected yet, all DICOM info is loaded into one variable

%% LOADING IMAGE DATA

if isfield(dcmInfo,'EchoTime')
    EchoTimesAll=[dcmInfo(:).EchoTime];
else
    EchoTimesAll=nan;
end
EchoTimesUni=unique(EchoTimesAll);
matObj.EchoTimesUni=EchoTimesUni;
NumEchoTimes=numel(EchoTimesUni);

%Get image types
ImageTypesAll={dcmInfo.ImageType};
ImageTypesUni=unique(ImageTypesAll);
matObj.ImageTypesUni=ImageTypesUni;
NumImageTypes=numel(ImageTypesUni);

if isfield(dcmInfo,'NumberOfTemporalPositions') %Multiple dynamics
    NumberOfTemporalPositions=double(dcmInfo(1).NumberOfTemporalPositions);
else
    NumberOfTemporalPositions=1;
end

if isfield(dcmInfo,'TemporalPositionIdentifier') %Multiple dynamics
    TemporalPositionIdentifierAll=[dcmInfo(:).TemporalPositionIdentifier];
else
    TemporalPositionIdentifierAll=[];
end

if isfield(dcmInfo,'MRImagePhaseNumber') %Multiple dynamics
    MRImagePhaseNumbersAll=[dcmInfo(:).MRImagePhaseNumber];
    MRImagePhaseNumbersUni=unique(MRImagePhaseNumbersAll);
    NumberOfPhasesMR=numel(MRImagePhaseNumbersUni);%dcmInfo(1).NumberOfPhasesMR %MRSeriesNrOfSlices
else
    MRImagePhaseNumbersAll=[];
    NumberOfPhasesMR=1;
end


switch dictSetting
    case 1 %PHILIPS
        SliceNumberMRAll=[dcmInfo(:).ImagePlaneNumber];%[dcmInfo(:).SliceNumberMR];
        NumberOfSlices=numel(unique(SliceNumberMRAll));%double(dcmInfo(1).NumberOfSlicesMR);
    case 2 % SIEMENS
        SliceNumberMRAll=[dcmInfo(:).InstanceNumber];        
        NumberOfSlices=numel(unique(SliceNumberMRAll));%double(dcmInfo(1).NumberOfSlicesMR);
        NumberOfSlices=NumberOfSlices./NumEchoTimes;
        SliceNumberMRAll(SliceNumberMRAll>NumberOfSlices)=SliceNumberMRAll(SliceNumberMRAll>NumberOfSlices)-NumberOfSlices;        
end

NumberOfFilesPerType=NumberOfSlices*NumberOfTemporalPositions;

NumberOfRows=double(dcmInfo(1).Rows);
NumberOfColumns=double(dcmInfo(1).Columns);

%%
ImageSize=double([NumberOfRows NumberOfColumns NumberOfSlices NumberOfTemporalPositions NumberOfPhasesMR]);
matObj.ImageSize=ImageSize;

NumClass=['uint',num2str(dcmInfo_full.BitsAllocated)];

c=0;
hw = waitbar(0,'Loading DICOM image data...');
for iEcho=1:NumEchoTimes
    
    %Finding files for current EchoTime
    EchoTimeNow=EchoTimesUni(iEcho); %The current echo time
    logicEcho=EchoTimesAll==EchoTimeNow;
    
    %String to add to type spec
    if NumEchoTimes==1
        echoNameAppend=[];
    else
        echoNameAppend=['_EchoTime_',num2str(iEcho)];
    end
    
    for iType=1:NumImageTypes
        
        %Finding files for current type
        ImageType=ImageTypesUni(iType);
        logicType=strcmpi(ImageTypesAll,ImageType);
        
        logicNow=logicType&logicEcho; 
        
        TypeFiles=files(logicNow);        
        TypeName=['type_',num2str(iType),echoNameAppend];
        TypeNameDcmInfo=[TypeName,'_info'];
        
        SliceNumberMRAll_Type=SliceNumberMRAll(logicNow);
        if isempty(TemporalPositionIdentifierAll)
            TemporalPositionIdentifierAll_Type=ones(size(logicNow));
        else
            TemporalPositionIdentifierAll_Type=TemporalPositionIdentifierAll(logicNow);
        end
        
        if isempty(MRImagePhaseNumbersAll)
            MRImagePhaseNumbersAll_Type=ones(size(logicNow));
        else
            MRImagePhaseNumbersAll_Type=MRImagePhaseNumbersAll(logicNow);
        end
        
        dcmInfo_Type=dcmInfo(logicNow);
        
        VarSize=(prod(ImageSize))*dcmInfo_full.BitsAllocated/8;
        if VarSize>MaxVarSize %Save to MAT-file for each temporal step
            
            eval(['matObj.',TypeName,'=zeros(ImageSize,NumClass)']);
            slicesVisited=false(1,ImageSize(3));
            for qFile=1:1:numel(TypeFiles);
                
                m=dicomread(fullfile(PathName,TypeFiles{qFile})); %Current slice
                
                %Indices for current slice
                indTemp=TemporalPositionIdentifierAll_Type(qFile);
                indSlice=SliceNumberMRAll_Type(qFile);
                indPhase=MRImagePhaseNumbersAll_Type(qFile);
                
                waitbar(c/numel(files),hw,['Loading and saving DICOM image data to MAT-file...',num2str(round(100.*c/numel(files))),'%']);
                eval(['matObj.',TypeName,'(:,:,indSlice,indTemp,indPhase)=m;']);
                c=c+1;
                slicesVisited(indSlice)=1;
            end
            eval(['matObj.',TypeNameDcmInfo,'=dcmInfo(logicType);']);
            
        else %create 4D array and save whole array to MAT-file
            N_info=repmat(dcmInfo(1),ImageSize(3:end)); %Allocate memory for info structure
            N=zeros(ImageSize,NumClass);
            slicesVisited=false(1,ImageSize(3));
            for qFile=1:1:numel(TypeFiles);
                
                m=dicomread(fullfile(PathName,TypeFiles{qFile})); %Current slice
                
                %Indices for current slice
                indTemp=TemporalPositionIdentifierAll_Type(qFile);
                indSlice=SliceNumberMRAll_Type(qFile);
                indPhase=MRImagePhaseNumbersAll_Type(qFile);
                
                N(:,:,indSlice,indTemp,indPhase)=m;
                N_info(indSlice,indTemp,indPhase)=dcmInfo_Type(qFile);
                waitbar(c/numel(files),hw,['Loading DICIM image data...',num2str(round(100.*c/numel(files))),'%']);
                c=c+1;
                slicesVisited(indSlice)=1;
            end
            waitbar(c/numel(files),hw,['Saving DICOM image data to MAT-file...',num2str(round(100.*c/numel(files))),'%']);
            eval(['matObj.',TypeName,'=N;']);
            eval(['matObj.',TypeNameDcmInfo,'=N_info;']);
            
        end     
    end
    if any(~slicesVisited)
        error('Missing slices');
    end
end
close(hw);

dicomdict('factory');


 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
