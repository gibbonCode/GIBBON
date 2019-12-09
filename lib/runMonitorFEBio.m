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
%------------------------------------------------------------------------

%% 

lineSep=repmat('%',1,45);
if FEBioRunStruct.disp_on==1
    disp(' '); %Empty line
    disp(lineSep); %Stripe
    disp(['--- STARTING FEBIO JOB --- ',datestr(now)]);
end

%% Removing pre-existing files (e.g. from previsous FEBio job with same name) 

%Remove log file
if exist(FEBioRunStruct.run_logname,'file')==2
    delete(FEBioRunStruct.run_logname);
    
    %Check if its gone
    if exist(FEBioRunStruct.run_logname,'file')==2
        error(['Deletion of ',FEBioRunStruct.run_logname,' not succesful, check user permissions']);
    end
end

% Remove .xplot file
[filePath,fileName,~]=fileparts(FEBioRunStruct.run_filename);

fileName_plot=fullfile(filePath,[fileName,'.xplt']);
if exist(fileName_plot,'file')==2
    delete(fileName_plot);
    %Check if its gone
    if exist(fileName_plot,'file')==2
        error(['Deletion of ',fileName_plot,' not succesful, check user permissions']);
    end
end

%Remove other requested files (e.g. output files)
if ~isfield(FEBioRunStruct,'cleanUpFileList')
    FEBioRunStruct.cleanUpFileList={}; 
end

if ~isempty(FEBioRunStruct.cleanUpFileList)
    for q=1:1:numel(FEBioRunStruct.cleanUpFileList)
        
        fileToRemove=FEBioRunStruct.cleanUpFileList{q};
        
        [fileToRemovePath,~,~]=fileparts(fileToRemove);
        
        if isempty(fileToRemovePath) %Since FEBio 2.4.0 path is the same as .feb files path
            fileToRemove=fullfile(filePath,fileToRemove); %Create full path name
        end
        
        %Remove file if it exists
        if exist(fileToRemove,'file')==2
            delete(fileToRemove);
            %Check if its gone
            if exist(fileToRemove,'file')==2
                error(['Deletion of ',fileToRemove,' not succesful, check user permissions']);
            end
        end        
    end
end

%%

if ~isfield(FEBioRunStruct,'runMode')
    FEBioRunStruct.runMode='external'; %Default behaviour runs FEBio "externally"
end

if ~isfield(FEBioRunStruct,'FEBioPath')
    FEBioRunStruct.FEBioPath=getFEBioPath;
end

if isempty(FEBioRunStruct.FEBioPath)    
    FEBioRunStruct.FEBioPath='febio'; %Assume febio is known as an internal variable
end

if ~isfield(FEBioRunStruct,'run_string')
    switch FEBioRunStruct.runMode
        case 'external_old'
            FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
            runExternal=1;
        case 'external'
            if ismac % Code to run on Mac plaform
                %TO DO IMPROVE THIS LINE FOR MAC
                FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
            elseif isunix && ~ismac % Code to run on Linux plaform                               
%                 FEBioRunStruct.run_string=['gnome-terminal -e ""',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"" &'];                
%                 FEBioRunStruct.run_string=['xterm -e ""',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"" &'];                
%                 FEBioRunStruct.run_string=['konsole -e ""',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"" &'];                
                FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];                
            elseif ispc % Code to run on Windows platform
                FEBioRunStruct.run_string=['start /min "GIBBON - FEBio" "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
            else
                warning('Unknown operational system, run command might be inappropriate');
                FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
            end
            
             runExternal=1;
        case 'internal'
            if ispc
                FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
            elseif isunix
                FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"']; 
            else
                FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
            end
            runExternal=0;
    end
end

if ~isfield(FEBioRunStruct,'run_string_quit')
    FEBioRunStruct.run_string_quit='taskkill /F /IM FEBio2.exe /T';
end

%% Starting FEBio job

FEBio_go=0; %i.e. not stopped
tPre=tic; 
% FEBioRunStruct.run_string
system(FEBioRunStruct.run_string); %START FEBio NOW!!!!!!!!
tPost=tic;

%% Wait for log file

if FEBioRunStruct.disp_on==1
    disp('Waiting for log file...');
end

fileFound=0;
while fileFound==0
    pause(FEBioRunStruct.t_check); %Wait
    if exist(FEBioRunStruct.run_logname,'file')
        [T_log]=txtfile2cell(FEBioRunStruct.run_logname);
        fileFound=1;
    end
    t_fea=toc(tPre); %Current time since start

    switch FEBioRunStruct.runMode
        case 'internal'
            
        otherwise                        
            if t_fea>=FEBioRunStruct.maxLogCheckTime %&& logTimingError==1
                runFlag=0;
                FEBio_go=1;
                if FEBioRunStruct.disp_on==1
                    warning(['--- FAILED: Log file was not created in time. FEBio likely failed prior to logfile creation! --- ',datestr(now)]);
                    warning('Try increasing FEBioRunStruct.maxLogCheckTime in case FEBio is simply taking too long to start')
                    if runExternal==1
                        warning('Try setting FEBioRunStruct.runMode to "internal" to see potential errors prior to log file creation');
                    end
                end
                break
            end
    end
end

%% Check log file

if FEBio_go==0
    
    if FEBioRunStruct.disp_on==1
        disp(['Proceeding to check log file...',datestr(now)]);
    end
    
    %Scan log file for the following targets
    targets={'------- converged at time :',' N O R M A L   T E R M I N A T I O N',' E R R O R   T E R M I N A T I O N'};
    line_count=1;
    numLinesPrevious=-1;
    
    tStartCheck=datevec(now); 
    while FEBio_go==0
        
        pause(FEBioRunStruct.t_check); %Wait
        
        %Import log file
        [T_log]=txtfile2cell(FEBioRunStruct.run_logname);
        
        if ~isempty(T_log) && numel(T_log)~=numLinesPrevious %Change in number of lines 
            tStartCheck=datevec(now); %reset change checking clock
            
            %Check log file content
            while 1
                if line_count>numel(T_log)
                    break;
                end
                
                %Get line
                l=T_log{line_count};
                
                %Display line
                if (strfind(l,targets{1}))
                    if FEBioRunStruct.disp_on==1 && FEBioRunStruct.disp_log_on==1
                        disp(l); %display line
                    end
                end
                
                %Check for normal termination
                if (strfind(l,targets{2})) % Found: N O R M A L  T E R M I N A T I O N
                    runFlag=1;
                    FEBio_go=1;
                    break
                end
                
                %Check for error termination
                if (strfind(l,targets{3})) % Found: E R R O R   T E R M I N A T I O N
                    if FEBioRunStruct.disp_on==1
                        disp(l); %display line
                    end
                    runFlag=0;
                    FEBio_go=1;
                    if FEBioRunStruct.disp_on==1
                        disp(['--- FAILED: FEBio error! Check log-file --- ',datestr(now)]);
                    end
                    break
                end
                line_count=line_count+1;
            end
        else
            tNow=datevec(now); 
            dt=etime(tNow,tStartCheck);            
            if dt>FEBioRunStruct.maxLogCheckTime
                runFlag=0;
                FEBio_go=1;
                if FEBioRunStruct.disp_on==1
                    warning(['--- FAILED: Log file was not changed within maxLogCheckTime. ',datestr(now)]);
%                     try
%                         system(FEBioRunStruct.run_string_quit);
%                     end
                end
                break
            end
        end
        numLinesPrevious=numel(T_log);
        
        switch FEBioRunStruct.runMode
            case 'internal'
                
            otherwise
                if FEBio_go==0
                    %Police timing
                    t_fea=toc(tPre); %Current time since start
                    [runFlag,FEBio_go]=policeFEBioJob(FEBioRunStruct,t_fea);
                end
        end
    end    
end

if FEBioRunStruct.disp_on==1
    if runFlag==1
        disp(['--- Done --- ',datestr(now)]);
    else
        disp(['--- Failed or ended! --- ',datestr(now)]);
    end
end

end

function [runFlag,FEBio_go]=policeFEBioJob(FEBioRunStruct,t_fea)

if t_fea>=FEBioRunStruct.maxtpi 
    if FEBioRunStruct.disp_on==1
        disp(['--- FAILED: Maximum analysis time exceeded, FEBio job will be terminated! --- ',datestr(now)]);
    end
    runFlag=0;
    FEBio_go=1;
    system(FEBioRunStruct.run_string_quit);        
else
    runFlag=0;
    FEBio_go=0;
end

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
