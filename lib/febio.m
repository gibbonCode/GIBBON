function febio(varargin)

% function febio(febFileName)
% ------------------------------------------------------------------------
% This function calls FEBio from MATLAB and runs it in the command window.
% If no input is provided this commend simply triggers FEBio to run without
% an input file. 
% If the input is a valid FEBio model input file path then FEBio will run
% the analysis for that file. 
% The input may also be a structure with the following default fields: 
%   defaultOptionStruct.FEBioPath=getFEBioPath; %The FEBio path
%   defaultOptionStruct.run_filename=[]; %The .feb file name
%   defaultOptionStruct.run_logname=[]; %The log file name
%   defaultOptionStruct.runMode='internal'; %Running internally
%
%
% 2020/11/24 Created
% ------------------------------------------------------------------------

%% Parse input

%Create default structure
defaultOptionStruct.FEBioPath=getFEBioPath; %Get FEBio path
defaultOptionStruct.run_filename=[];
defaultOptionStruct.run_logname=[];
defaultOptionStruct.runMode='internal';

switch nargin
    case 0
        inputVar=[];
    case 1        
        inputVar=varargin{1};
end

if ~isempty(inputVar)
    if ~isstruct(inputVar) %Assume a feb file name was provided
        FEBioRunStruct.run_filename=inputVar;
        [filePath,fileName,~]=fileparts(inputVar);
        logFileName=fullfile(filePath,[fileName,'.txt']);
        FEBioRunStruct.run_logname=logFileName;    
    else
        FEBioRunStruct=inputVar;
    end    
    %Complete struct if there are missing fields
    FEBioRunStruct=structComplete(FEBioRunStruct,defaultOptionStruct,0);
else
    FEBioRunStruct=[];
end

%% Run FEBio

if isstruct(FEBioRunStruct) && ~isfield(FEBioRunStruct,'run_string')
    switch FEBioRunStruct.runMode
        case 'external_old'
            FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
        case 'external'
            if ismac % Code to run on Mac plaform
                %TO DO IMPROVE THIS LINE FOR MAC
                FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'" &'];
            elseif isunix && ~ismac % Code to run on Linux plaform                               
                FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'" &'];                
            elseif ispc % Code to run on Windows platform
                FEBioRunStruct.run_string=['start /min "GIBBON - FEBio" "',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'"'];
            else
                warning('Unknown operational system, run command might be inappropriate');
                FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'" &'];
            end            
        case 'internal'
            if ispc
                FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'"'];
            elseif isunix
                FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'"']; 
            else
                FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" "',FEBioRunStruct.run_filename,'"'];
            end
    end
elseif isempty(FEBioRunStruct)
    FEBioRunStruct.run_string=['"',defaultOptionStruct.FEBioPath,'"'];
end

% if isstruct(FEBioRunStruct) && ~isfield(FEBioRunStruct,'run_string')
%     switch FEBioRunStruct.runMode
%         case 'external_old'
%             FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
%         case 'external'
%             if ismac % Code to run on Mac plaform
%                 %TO DO IMPROVE THIS LINE FOR MAC
%                 FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
%             elseif isunix && ~ismac % Code to run on Linux plaform                               
%                 FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];                
%             elseif ispc % Code to run on Windows platform
%                 FEBioRunStruct.run_string=['start /min "GIBBON - FEBio" "',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
%             else
%                 warning('Unknown operational system, run command might be inappropriate');
%                 FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'"',' -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'" &'];
%             end            
%         case 'internal'
%             if ispc
%                 FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
%             elseif isunix
%                 FEBioRunStruct.run_string=['nice "',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"']; 
%             else
%                 FEBioRunStruct.run_string=['"',FEBioRunStruct.FEBioPath,'" -i "',FEBioRunStruct.run_filename,'" -o "',FEBioRunStruct.run_logname,'"'];
%             end
%     end
% elseif isempty(FEBioRunStruct)
%     FEBioRunStruct.run_string=['"',defaultOptionStruct.FEBioPath,'"'];
% end

% FEBioRunStruct.run_string

system(FEBioRunStruct.run_string); %START FEBio NOW!!!!!!!!

end