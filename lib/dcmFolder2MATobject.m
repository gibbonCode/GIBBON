function dcmFolder2MATobject(varargin)

% function dcmFolder2MATobject(PathName,MaxVarSize,reOrderOpt,dicomDictType,fileExtension)
% ------------------------------------------------------------------------
% This function converts the DICOM files in the folder PathName to a MAT
% object which is stored as a file called IMDAT in a new subfolder called
% IMDAT.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log
% 2018/06/11 Suppressed/removed triggertime related warning (this appears
% to be a custom field which is nearly always missing).   
%------------------------------------------------------------------------

%% PARSE INPUT

switch nargin
    case 1
        PathName=varargin{1};
        MaxVarSize=[];
        reOrderOpt=0;
        dicomDictFactory=0;
        fileExtension='.dcm';
        optionStruct=[];
    case 2
        PathName=varargin{1};
        MaxVarSize=varargin{2};
        reOrderOpt=0;
        dicomDictFactory=0;
        fileExtension='.dcm';
        optionStruct=[];
    case 3
        PathName=varargin{1};
        MaxVarSize=varargin{2};
        reOrderOpt=varargin{3};
        dicomDictFactory=0;
        fileExtension='.dcm';
        optionStruct=[];
    case 4
        PathName=varargin{1};
        MaxVarSize=varargin{2};
        reOrderOpt=varargin{3};
        dicomDictFactory=varargin{4};
        fileExtension='.dcm';
    case 5
        PathName=varargin{1};
        MaxVarSize=varargin{2};
        reOrderOpt=varargin{3};
        dicomDictFactory=varargin{4};
        fileExtension=varargin{5};   
        optionStruct=[];
    case 6
        PathName=varargin{1};
        MaxVarSize=varargin{2};
        reOrderOpt=varargin{3};
        dicomDictFactory=varargin{4};
        fileExtension=varargin{5};
        optionStruct=varargin{6};
end

defaultOptionStruct.ignoreDynamic=0;

%Fix option structure, complete and remove empty values
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

%%
%Get/set maximum variable size (sets save steps in case the variable size
%is exceeded)
if isempty(MaxVarSize)
    MaxVarSize=1e9; %Default
end

if isempty(reOrderOpt)
    reOrderOpt=0; %Default
end

if isempty(dicomDictFactory)
    dicomDictFactory=0; %Default
end

if isempty(fileExtension)
    fileExtension='.dcm'; %Default
end

%% Getting DICOM file names and image dimension parameters
files = dir(fullfile(PathName,['*',fileExtension]));
files={files(1:end).name};
NumberOfFiles=numel(files);

%% Attempt to sort files
% This sorting assumes the files names make sense in terms of slice
% numbers. 

try %Try to get numbers from image file names
    nameNumbers=cellfun(@str2double,regexp(files,'\d+\d*','match'));
    [~,indSort]=sort(nameNumbers); %Get sort order
    files=files(indSort); %Sort file list based on numeric data in names
catch %Normal sort instead
    files=sort(files(:)); %Attempt sort (might cause dcm1,dcm10,...)
end

%%

firstWarning_TriggerTime=0;
firstWarning_Rotation=0;
if NumberOfFiles>0
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
        'AcquisitionMatrix'...
        'PixelBandwidth'...
        'ImagingFrequency'...
        'MagneticFieldStrength'...
        'InPlanePhaseEncodingDirection'...
        'PhaseEncodingDirection'...
        };
    % 'EchoNumbers' is private and just a echotime identifier or sequence number
    % 'SliceNumberMR' is also private (2001,100a) and can be usefull for identifing slices that are on same position
    % N.B. Filename is added by matlab including the path, Basename is one we add ourself
    
    %% Creating output file name and folder
    fName='IMDAT';
    foldername_out=fullfile(PathName,fName,filesep);
    if ~exist(foldername_out,'dir') %create output folder if it does not exist already
        mkdir(foldername_out); %Create new dir
    end
    savename_out1=fullfile(foldername_out,'IMDAT.mat');
    if exist(savename_out1,'file') %Check if file exist already
        delete(savename_out1); %Remove file
    end
    matObj = matfile(savename_out1,'Writable',true);
    
    %% SETTING DICOM DICTIONARY
    fName=fullfile(PathName,files{1});
    
    %First import using factory settings
    dicomdict('factory');
    dcmInfo_full=dicominfo(fName);  
    m=dicomread(fName); %Test upload first
    sizSlice=size(m);
    dictSetting=3;
    
    if dicomDictFactory==0
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
            elseif ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'GE MEDICAL SYSTEMS'))
                disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
                warning('No settings for this vendor using DICOM dictionary factory settings');
                dictSetting=3;
            else
                dicomdict('factory');
                warning('Unknown vendor, using DICOM dictionary factory settings');
                dictSetting=3;
            end
            %Test info import with new dictionary (if this fails we resume with
            % factory
            dcmInfo_full=dicominfo(fName);
        catch %e.g. if the manufacturer field is missing
            dicomdict('factory');
            warning('Unknown vendor, using DICOM dictionary factory settings');
            dictSetting=3;
            dcmInfo_full=dicominfo(fName);
        end
    elseif dicomDictFactory==1
%         dicomdict('factory');
%         dcmInfo_full=dicominfo(fName);
%         dictSetting=3;
    end
        
    %% LOADING DICOM INFO
    hw = waitbar(0,'Loading DICOM info...');    
    for c=1:1:numel(files)
        fName=fullfile(PathName,files{c});
        %dicomdict('get');
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
    
    %% Reorder files based on InstanceNumber (file naming may deviate from slice order) 
    if reOrderOpt
        [~,indSort]=sort([dcmInfo(:).InstanceNumber]);        
        files=files(indSort);
        dcmInfo=dcmInfo(indSort);
    end
    
    %%

    %Geometry information based on first info file
    [G]=dicom3Dpar(dcmInfo(1));
    matObj.G=G; %Saving geometry information
    
    %TODO this part is not memory protected yet, all DICOM info is loaded into
    %one variable
    
    % hw = waitbar(0,'Loading DICOM info...');
    % step_size=100;
    % q=1;
    % ind_save=[];
    % for c=1:1:numel(files)
    %     ind_save=[ind_save c];
    %     dcmInfo(q)=dicominfo(fullfile(foldername,files{c}));
    %     waitbar(c/numel(files),hw,['Loading DICIM info...',num2str(round(100.*c/numel(files))),'%']);
    %     if numel(ind_save)>=step_size
    %         waitbar(c/numel(files),hw,['Saving to MAT-file...',num2str(round(100.*c/numel(files))),'%']);
    %         matObj.dcmInfo(1,ind_save)=dcmInfo;
    %         ind_save=[];
    %         dcmInfo=[];
    %         q=1;
    %     end
    %     q=q+1;
    % end
    % matObj.dcmInfo(1,ind_save)=dcmInfo;
    % close(hw);
    
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
   
    if isfield(dcmInfo,'TriggerTime') %Multiple dynamics
        TriggerTimesAll=[dcmInfo(:).TriggerTime];
        TriggerTimesUni=unique(TriggerTimesAll);
        NumberOfTemporalPositions=numel(TriggerTimesUni); %NumberOfPhasesMR
    else
        TriggerTimesAll=zeros(1,NumberOfFiles);
        TriggerTimesUni=unique(TriggerTimesAll);
        if isfield(dcmInfo,'NumberOfTemporalPositions') %Multiple dynamics
            NumberOfTemporalPositions=double(dcmInfo(1).NumberOfTemporalPositions);
        else
            NumberOfTemporalPositions=1;
        end
    end
    
    if optionStruct.ignoreDynamic==1   
        NumberOfTemporalPositions=1;
    end
    
    NumberOfSlices=NumberOfFiles/NumImageTypes/NumEchoTimes/NumberOfTemporalPositions; %NumberOfSlicesMR           
%     NumberOfFilesPerType=NumberOfSlices*NumberOfTemporalPositions;
    
    switch dictSetting
        case 1 %PHILIPS
            %             NumberOfRows=double(dcmInfo(1).Width);
            %             NumberOfColumns=double(dcmInfo(1).Height);
            NumberOfRows=double(dcmInfo(1).Rows);
            NumberOfColumns=double(dcmInfo(1).Columns);
        case 2 % SIEMENS
            NumberOfRows=double(dcmInfo(1).Rows);
            NumberOfColumns=double(dcmInfo(1).Columns);
        case 3 %FACTORY
            NumberOfRows=double(dcmInfo(1).Rows);
            NumberOfColumns=double(dcmInfo(1).Columns);
            %             NumberOfRows=double(dcmInfo(1).Width);
            %             NumberOfColumns=double(dcmInfo(1).Height);
    end
    
    ImageSize=[NumberOfRows NumberOfColumns NumberOfSlices NumberOfTemporalPositions];
    
%     %Check if slice size matching DICOM parameters   
%     if sizSlice(1)~=ImageSize(1)
%         ImageSize(1)=sizSlice(1);
%         warning('Image size 2 does not match size specified in DICOM info. Using image size from image instead')
%     end
%     if sizSlice(2)~=ImageSize(2)
%         ImageSize(2)=sizSlice(2);
%         warning('Image size 2 does not match size specified in DICOM info. Using image size from image instead')
%     end
        
    matObj.ImageSize=ImageSize;
    
    NumClass=['uint',num2str(dcmInfo_full.BitsAllocated)];
    
    c=0;
    hw = waitbar(0,'Loading DICOM info...');
    
    for iEcho=1:NumEchoTimes
        %Finding files for current EchoTime
        EchoTimeNow=EchoTimesUni(iEcho); %The current echo time
        
        L_Echo=EchoTimesAll==EchoTimeNow;
        
        %String to add to type spec
        if NumEchoTimes==1
            echoNameAppend=[];
        else
            echoNameAppend=['_EchoTime_',num2str(iEcho)];
        end
        
        for iType=1:NumImageTypes
            
            %Finding files for current type
            ImageTypeNow=ImageTypesUni(iType); %The current image type
            L_Type=strcmpi(ImageTypesAll,ImageTypeNow);
            
            %Fix L_Echo in case EchoTIme is not defined
            if isnan(EchoTimesAll)
                L_Echo=true(size(L_Type));
            end
            
            L_now=L_Echo&L_Type;
            
            TypeFiles=files(L_now);
%             iStart=((iType-1)*NumberOfFilesPerType)+1; %     iEnd=iStart-1+NumberOfFilesPerType;
            
            TypeName=['type_',num2str(iType),echoNameAppend];
            TypeNameDcmInfo=[TypeName,'_info'];
            
            VarSize=(prod(ImageSize))*dcmInfo_full.BitsAllocated/8;
            if VarSize>MaxVarSize %Save to MAT-file for each temporal step
                if ImageSize(4)>1
                    matObj.(TypeName)=zeros([ImageSize(1:3) 2],NumClass);
                else
                    matObj.(TypeName)=zeros(ImageSize,NumClass);
                end
                for iTemp=1:1:NumberOfTemporalPositions
                    %Finding files for current temporal position
                    TriggerTimesType=TriggerTimesAll(L_now);%[dcmInfo(L_now).TriggerTime];
                    TriggerTime=TriggerTimesUni(iTemp);
                    L_Temp=(TriggerTimesType==TriggerTime);
                    %Creating volume by assembling slices
                    TempTypeFiles=TypeFiles(L_Temp);
                    M=zeros(ImageSize(1:3),NumClass);
                    for iSlice=1:1:numel(TempTypeFiles)
                        waitbar(c/numel(files),hw,['Loading DICIM image data...',num2str(round(100.*c/numel(files))),'%']);
                        load_name=fullfile(PathName,TempTypeFiles{iSlice});
                        m=dicomread(load_name);                        
                        if size(m,3)>1                            
                            m=mean(double(m),3);
                            warning('Multi-dimensional (e.g. RGB) slices where converted to 2D grayscale');
                        end                        
                        M(:,:,iSlice)=m;
                        
                        c=c+1;
                    end
                    waitbar(c/numel(files),hw,['Saving DICOM image data to MAT-file...',num2str(round(100.*c/numel(files))),'%']);
                    eval(['matObj.',TypeName,'(:,:,:,iTemp)=M;']);
                end
            else %create 4D array and save whole array to MAT-file
                N=zeros(ImageSize,NumClass);
                for iTemp=1:1:NumberOfTemporalPositions
                    %Finding files for current temporal position
                    TriggerTimesType=TriggerTimesAll(L_now);%[dcmInfo(L_now).TriggerTime];
                    TriggerTime=TriggerTimesUni(iTemp);
                    L_Temp=(TriggerTimesType==TriggerTime);
                    %Creating volume by assembling slices
                    TempTypeFiles=TypeFiles(L_Temp);
                    M=zeros(ImageSize(1:3),NumClass);
                    for iSlice=1:1:numel(TempTypeFiles)
                        waitbar(c/numel(files),hw,['Loading DICIM image data...',num2str(round(100.*c/numel(files))),'%']);
                        load_name=fullfile(PathName,TempTypeFiles{iSlice});
                        m=dicomread(load_name); %Current slice
                        if size(m,3)>1                            
                            m=mean(double(m),3);
                            warning('Multi-dimensional (e.g. RGB) slices where converted to 2D grayscale');
                        end
                        if size(M,1)~=size(m,1)
                            M(:,:,iSlice)=m'; %Try transpose                            
                            if firstWarning_Rotation==0
                                warning('Image data was rotated to match expected image size!');
                                firstWarning_Rotation=1;
                            end
                        else
                            M(:,:,iSlice)=m;
                        end
                        c=c+1;
                    end 
                    N(:,:,:,iTemp)=M;
                end
                waitbar(c/numel(files),hw,['Saving DICOM image data to MAT-file...',num2str(round(100.*c/numel(files))),'%']);
                eval(['matObj.',TypeName,'=N;']);
                eval(['matObj.',TypeNameDcmInfo,'=dcmInfo(L_now);']);
                
            end
        end
    end
    close(hw);
else    
    warning(['Folder ',PathName,' skipped since it does not contain .dcm files']);
end

dicomdict('factory');

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
