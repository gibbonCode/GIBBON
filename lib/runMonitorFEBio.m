function [runFlag]=runMonitorFEBio(FEBioRunStruct)

% function [runFlag]=runMonitorFEBio(FEBioRunStruct)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2013/04/01 Updated for GIBBON
% 2016/05/23 Fixed bug in relation to path names for FEBioRunStruct.cleanUpFileList
% 2020/11/24 Cleaned up and simplified code.
% 2020/11/24 Improved handling of internal/external mode
% 2020/11/24 Uses febio.m code to trigger FEBio
% 2020/11/25 Changed command window progress text
%------------------------------------------------------------------------

%% Parse input

%Create default structure
defaultFEBioRunStruct.FEBioPath=getFEBioPath; %Get FEBio path
defaultFEBioRunStruct.run_filename=[];
defaultFEBioRunStruct.run_logname=[];
defaultFEBioRunStruct.runMode='external';
defaultFEBioRunStruct.cleanUpFileList={};
defaultFEBioRunStruct.disp_on=1; %Display information on the command window
defaultFEBioRunStruct.t_check=0.5; %Time for checking log file (dont set too small)
defaultFEBioRunStruct.maxtpi=inf; %Max analysis time
defaultFEBioRunStruct.maxLogCheckTime=30; %Max log file checking time

%Construct default task kill command
[~,processName,exeExt]=fileparts(getFEBioPath);
if ispc
    defaultFEBioRunStruct.run_string_quit=['taskkill /F /IM ',processName,exeExt,' /T'];
else
    defaultFEBioRunStruct.run_string_quit=['killall ',processName];
end

%Complete struct if there are missing fields
FEBioRunStruct=structComplete(FEBioRunStruct,defaultFEBioRunStruct,0);


%Check if febio path is empty
if isempty(FEBioRunStruct.FEBioPath)
    FEBioRunStruct.FEBioPath='febio'; %Assume febio is known as an internal variable
end

%Check if logfile is empty and set it to default form if needed
if ~isempty(FEBioRunStruct.run_filename) && isempty(FEBioRunStruct.run_logname)
    [filePath,fileName,~]=fileparts(FEBioRunStruct.run_filename);
    FEBioRunStruct.run_logname=fullfile(filePath,[fileName,'.txt']);
end

%Check if log file path is empty and replace with feb file path if needed
[logFilePath,~,~]=fileparts(FEBioRunStruct.run_logname);
if isempty(logFilePath)
    [filePath,~,~]=fileparts(FEBioRunStruct.run_filename);
    FEBioRunStruct.run_logname=fullfile(filePath,FEBioRunStruct.run_logname);
end

%% Check for old FEBio version

if contains(lower(FEBioRunStruct.FEBioPath),'febio2')    
    warning('FEBio2 detected. FEBio2 support is depricated. Please upgrade to FEBio3'); 
end

%% Display start message

if FEBioRunStruct.disp_on==1
    startString=['-------->    RUNNING/MONITORING FEBIO JOB    <-------- ',datestr(now)];    
    stringLength=numel(startString);
    lineSep=repmat('%',1,stringLength); %Stripe of % symbols
    disp(' '); %Empty line
    disp(lineSep);
    disp(startString);
    disp(['FEBio path: ',FEBioRunStruct.FEBioPath])
end

%% Removing pre-existing files (e.g. from previous FEBio job with same name)

%Remove log file
if exist(FEBioRunStruct.run_logname,'file')==2
    if FEBioRunStruct.disp_on==1
        dispMessage('# Attempt removal of existing log files',stringLength);        
    end
    delete(FEBioRunStruct.run_logname);
    
    %Check if its gone
    if exist(FEBioRunStruct.run_logname,'file')==2
        error(['Deletion of ',FEBioRunStruct.run_logname,' not succesful, check user permissions']);
    else
        if FEBioRunStruct.disp_on==1
            dispMessage(' * Removal succesful',stringLength);            
        end
    end
end

% Remove .xplot file
[filePath,fileName,~]=fileparts(FEBioRunStruct.run_filename);

fileName_plot=fullfile(filePath,[fileName,'.xplt']);
if exist(fileName_plot,'file')==2
    if FEBioRunStruct.disp_on==1
        dispMessage('# Attempt removal of existing .xplt files',stringLength);
    end
    delete(fileName_plot);
    %Check if its gone
    if exist(fileName_plot,'file')==2
        error(['Deletion of ',fileName_plot,' not succesful, check user permissions']);
    else
        if FEBioRunStruct.disp_on==1
            dispMessage(' * Removal succesful',stringLength);             
        end
    end
end

%Remove other requested files (e.g. output files)
if ~isempty(FEBioRunStruct.cleanUpFileList)
    for q=1:1:numel(FEBioRunStruct.cleanUpFileList)
        
        fileToRemove=FEBioRunStruct.cleanUpFileList{q};
        
        [fileToRemovePath,~,~]=fileparts(fileToRemove);
        
        if isempty(fileToRemovePath) %Since FEBio 2.4.0 path is the same as .feb files path
            fileToRemove=fullfile(filePath,fileToRemove); %Create full path name
        end
        
        %Remove file if it exists
        if exist(fileToRemove,'file')==2
            if FEBioRunStruct.disp_on==1
                dispMessage(['# Attempt removal of user defined files'],stringLength);                 
            end
            delete(fileToRemove);
            
            %Check if its gone
            if exist(fileToRemove,'file')==2
                error(['Deletion of ',fileToRemove,' not succesful, check user permissions']);
            else
                if FEBioRunStruct.disp_on==1
                    dispMessage(' * Removal succesful',stringLength); 
                end
            end
        end
    end
end

%% Starting FEBio job

if FEBioRunStruct.disp_on==1
    dispMessage('# Starting FEBio... ',stringLength);     
    disp(['  Max. total analysis time is: ',num2str(FEBioRunStruct.maxtpi),' s']); 
end

FEBio_monitor=1; %i.e. currently monitoring
tPre=tic; %Initiate pre-start timer
febio(FEBioRunStruct); %Start FEBio
%tPost=tic; %Initiate post-run timer

%% Wait for log file

if any(strcmp(FEBioRunStruct.runMode,{'external_old','external'}))
    if FEBioRunStruct.disp_on==1
        dispMessage(' * Waiting for log file creation',stringLength);   
        disp(['   Max. wait time: ',num2str(FEBioRunStruct.maxLogCheckTime),' s']);
    end
    
    %Check for creation of log file
    fileFound=0;
    while fileFound==0
        %Check for existance of log file
        if exist(FEBioRunStruct.run_logname,'file')
            fileFound=1;
            if FEBioRunStruct.disp_on==1
                dispMessage(' * Log file found.',stringLength);
            end
        end
        
        if fileFound==0 %Check if it is taking too long
            t_fea=toc(tPre); %Current time since start
            
            if t_fea>=FEBioRunStruct.maxLogCheckTime 
                runFlag=0;
                FEBio_monitor=0;
                if FEBioRunStruct.disp_on==1
                    warning(['--- FAILED: Log file was not created in time. FEBio likely failed prior to logfile creation! --- ',datestr(now)]);
                    warning('Try increasing FEBioRunStruct.maxLogCheckTime in case FEBio is simply taking too long to start');
                    warning('Try setting FEBioRunStruct.runMode to "internal" to see potential errors prior to log file creation');
                end
                break
            end
            pause(FEBioRunStruct.t_check); %Wait a moment before checking again
        end
    end
else
    if exist(FEBioRunStruct.run_logname,'file')
        if FEBioRunStruct.disp_on==1
            dispMessage(' * Log file found.',stringLength);
        end
    else
        runFlag=0;
        FEBio_monitor=0;
        warning(['--- FAILED: Log file was not created in time. FEBio likely failed prior to logfile creation! --- ',datestr(now)]);
    end
end

%% Check log file

if FEBio_monitor==1
    
    if FEBioRunStruct.disp_on==1
        dispMessage('# Parsing log file...',stringLength);           
    end
    
    %Scan log file for the following targets
    targetsConvergence={'number of iterations   :','number of reformations :','------- converged at time :',' Elapsed time'};
    targetsEnd={' N O R M A L   T E R M I N A T I O N',....
                ' E R R O R   T E R M I N A T I O N',...
                '* Model initialization failed'};

    numLinesPrevious=0;    
    tStartCheck=datevec(now);
    while FEBio_monitor==1
        
        pause(FEBioRunStruct.t_check); %Wait
        
        %Import log file
        [T_log]=txtfile2cell(FEBioRunStruct.run_logname);
        
        %Check log file content
        if ~isempty(T_log) && numel(T_log)>numLinesPrevious %Change in number of lines
            tStartCheck=datevec(now); %reset change checking clock
            
            % Show convergence/elapsed time data
            T_check=T_log(numLinesPrevious+1:end);
            logicConverged=gcontains(T_check,targetsConvergence);
            if any(logicConverged)
                T_show=T_check(logicConverged);
                if FEBioRunStruct.disp_on==1
                    for q=1:numel(T_show)
                        dispMessage(T_show{q},stringLength);
                    end
                end
            end
            
            %Show termination data
            logicTermination=gcontains(T_check,targetsEnd); 
            if any(logicTermination)
                T_show=T_check(logicTermination);
                if FEBioRunStruct.disp_on==1
                    for q=1:numel(T_show)
                        disp(T_show{q});
                    end
                end
                FEBio_monitor=0; %Stop on normal or error termination
                
                if any(gcontains(T_show,targetsEnd{1}))
                    runFlag=1; %Signal for normal termination
                end
                if any(gcontains(T_show,targetsEnd(2:end)))
                    runFlag=0; %Signal for error termination
                end
            end
        else
            tNow=datevec(now);
            dt=etime(tNow,tStartCheck);
            if dt>FEBioRunStruct.maxLogCheckTime
                runFlag=0;
                FEBio_monitor=0;
                if FEBioRunStruct.disp_on==1
                    warning(['--- FAILED: Log file was not changed within maxLogCheckTime. ',datestr(now)]);
                    try
                        system(FEBioRunStruct.run_string_quit);
                    catch
                    end
                end
                break
            end
        end
        numLinesPrevious=numel(T_log); %Update number of lines
                
        if any(strcmp(FEBioRunStruct.runMode,{'external_old','external'})) && FEBio_monitor==1
            %Police timing
            t_fea=toc(tPre); %Current time since start
            [runFlag,FEBio_monitor]=policeFEBioJob(FEBioRunStruct,t_fea);
        end
    end
end

if FEBioRunStruct.disp_on==1
    if runFlag==1
        dispMessage('# Done ',stringLength);   
        disp(lineSep);
    else
        dispMessage('# Failed or terminated! ',stringLength);        
        disp(lineSep);
    end
end

end


function [runFlag,FEBio_monitor]=policeFEBioJob(FEBioRunStruct,t_fea)

if t_fea>=FEBioRunStruct.maxtpi
    if FEBioRunStruct.disp_on==1
        disp('* Maximum analysis time exceeded');
        disp(['FEBio job will be terminated! ',datestr(now)]);
    end
    runFlag=0;
    FEBio_monitor=0;
    system(FEBioRunStruct.run_string_quit);
else
    runFlag=0;
    FEBio_monitor=1;
end

end

function dispMessage(m,stringLength)
d=datestr(now);
nRep=stringLength-numel(d)-numel(m);
if nRep>0
    s=repmat(' ',1,nRep);
else
    s=' ';
end
disp([m,s,d]);
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
