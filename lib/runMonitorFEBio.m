function [runFlag]=runMonitorFEBio(FEBioRunStruct)

%%
lineSep=repmat('%',1,45);
if FEBioRunStruct.disp_on==1;
    disp(' '); %Empty line
    disp(lineSep); %Stripe
    disp(['--- STARTING FEBIO JOB --- ',datestr(now)]);
end

%% Removing existing log-file
if exist(FEBioRunStruct.run_logname,'file')
    delete(FEBioRunStruct.run_logname);
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
            FEBioRunStruct.run_string=['start /min "GIBBON - FEBio" "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
            runExternal=1;
        case 'internal'
            FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
            runExternal=0;
    end
end

if ~isfield(FEBioRunStruct,'run_string_quit')
    FEBioRunStruct.run_string_quit='taskkill /F /IM FEBio.exe /T';
end

%% Starting FEBio job

FEBio_go=0; %i.e. not stopped

system(FEBioRunStruct.run_string); %START FEBio NOW!!!!!!!!
tic; %Start time (after starting FEBio to cope with "internal" mode)

%% Wait for log file

if FEBioRunStruct.disp_on==1;
    disp('Waiting for log file...');
end

fileFound=0;
while fileFound==0
    pause(FEBioRunStruct.t_check); %Wait
    if exist(FEBioRunStruct.run_logname,'file')
        [T_log]=txtfile2cell(FEBioRunStruct.run_logname);
        fileFound=1;
    end

    t_fea=toc; %Current time since start
    if t_fea>=FEBioRunStruct.maxLogCheckTime %&& logTimingError==1
        runFlag=0;
        FEBio_go=1;
        if FEBioRunStruct.disp_on==1; 
            disp(['--- FAILED: Log file was not created in time. FEBio likely failed proir to logfile creation! --- ',datestr(now)]);
            if runExternal==1
                disp('Try setting FEBioRunStruct.runMode to "internal" to see potential errors prior to log file creation');
            end
        end
        break
    end
end

%% Check log file

if FEBio_go==0
    
    if FEBioRunStruct.disp_on==1;
        disp(['Proceeding to check log file...',datestr(now)]);
    end
    
    %Scan log file for the following targets
    targets={'------- converged at time :',' N O R M A L   T E R M I N A T I O N',' E R R O R   T E R M I N A T I O N'};
    line_count=1;
    while FEBio_go==0
        
        pause(FEBioRunStruct.t_check); %Wait
        
        %Import log file
        [T_log]=txtfile2cell(FEBioRunStruct.run_logname);
        
        %Check log file content
        while 1
            if line_count>numel(T_log);
                break;
            end
            
            %Get line
            l=T_log{line_count};
            
            %Display line
            if (strfind(l,targets{1}))
                if FEBioRunStruct.disp_on==1 && FEBioRunStruct.disp_log_on==1;
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
                if FEBioRunStruct.disp_on==1;
                    disp(l); %display line
                end
                runFlag=0;
                FEBio_go=1;
                if FEBioRunStruct.disp_on==1;
                    disp(['--- FAILED: FEBio error! Check log-file --- ',datestr(now)]);
                end
                break
            end
            line_count=line_count+1;
        end
        
        if FEBio_go==0
            %Police timing
            t_fea=toc; %Current time since start
            [runFlag,FEBio_go]=policeFEBioJob(FEBioRunStruct,t_fea);
        end
        
    end    
end

if FEBioRunStruct.disp_on==1;
    if runFlag==1
        disp(['--- Done --- ',datestr(now)]);
    else
        disp(['--- Failed or ended! --- ',datestr(now)]);
    end
end

end

function [runFlag,FEBio_go]=policeFEBioJob(FEBioRunStruct,t_fea)

if t_fea>=FEBioRunStruct.maxtpi 
    if FEBioRunStruct.disp_on==1;
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

