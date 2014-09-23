function [mFileName_FOR,mFileName_PARFOR]=parsefor(mFileNameInput,cleanOpt,parOn_target,parOff_target,parEx_target)

% function parsefor(mFileNameInput,parOn_target,parOff_target)
% ------------------------------------------------------------------------
%
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/06/2013
%------------------------------------------------------------------------


%% Save settings
mFileName_FOR=[mFileNameInput(1:end-2),'_FOR.m'];
mFileName_PARFOR=[mFileNameInput(1:end-2),'_PARFOR.m'];

%% LOAD M-file
fid=fopen(mFileNameInput);
T=textscan(fid,'%s','delimiter', '\n','Whitespace','','bufsize',1e6);
T=T{1,1};
fclose(fid);

%% Initialize output parameters

%Output cells
T_FOR=T;
T_PARFOR=T;

%Logic variables for output
logicFor=true(numel(T),1);
logicParFor=true(numel(T),1);

%% Parsing M-file

%Form of wrap target strings
if nargin==2 %use default target types
    parOn_target='%% //parOn';
    parOff_target='%% //parOff';
    parEx_target='%% //parEx';
end

%Initialize actions
commentify=-1; %-1 is off, 1 is on
unCommentify=-1; %-1 is off, 1 is on
for lineIndex=1:1:numel(T)
    
    currentLine=T{lineIndex}; %The current line

    %Commenting code out
    if commentify==1 
        switch cleanOpt
            case 'c'
                logicFor(lineIndex)=0; %Do not include this line in FOR file
            case 'd'
                T_FOR(lineIndex)={['%',currentLine]}; %add percentage sign
        end
    end
    
    %Uncommenting by assuming first character is % sign
    if unCommentify==1        
        indCommentSigns=regexp(currentLine,'\w*%');
        if ~isempty(indCommentSigns)
            currentLineFix=currentLine;
            currentLineFix(indCommentSigns(1))='';
            T_FOR(lineIndex)={currentLineFix};%{currentLine(2:end)};
        end
        switch cleanOpt
            case 'c'
                logicParFor(lineIndex)=0; %Do not include this line in PARFOR file
            case 'd'                
                T_PARFOR(lineIndex)={['%',currentLine]}; %add percentage sign
        end      
    end
    
    %Testing current line for the parOn target 
    if strfind(currentLine,parOn_target) 
        commentify=-commentify; %Flip sign of commentify to turn on/off
        switch cleanOpt
            case 'c'
                logicFor(lineIndex)=0; %Do not include this line in FOR file
                logicParFor(lineIndex)=0; %Do not include this line in PARFOR file
            case 'd'
        end
    end
    
    %Testing current line for the parOff target
    if strfind(currentLine,parOff_target)
        unCommentify=-unCommentify; %Flip sign of commentify to turn on/off
        switch cleanOpt
            case 'c'
                logicFor(lineIndex)=0; %Do not include this line in FOR file
                logicParFor(lineIndex)=0; %Do not include this line in PARFOR file
            case 'd'
        end
    end
    
    %Testing current line for the parEx target 
    if strfind(currentLine,parEx_target)
        switch cleanOpt
            case 'c'
                logicFor(lineIndex)=0; %Do not include this line in FOR file
                logicParFor(lineIndex)=0; %Do not include this line in PARFOR file
            case 'd'
                T_FOR(lineIndex)={['%',currentLine]}; %add percentage sign
                T_PARFOR(lineIndex)={['%',currentLine]}; %add percentage sign
        end        
    end
    
end

%% Saving FOR M-file
fid=fopen(mFileName_FOR,'wt');
for lineIndexSave=1:size(T,1); 
    if logicFor(lineIndexSave)
        fprintf(fid,'%s\n',T_FOR{lineIndexSave});
    end
end
%Add last lines for information
fprintf(fid,'%s\n','%% --- PARSEFOR for-loop compatible file ---');
fprintf(fid,'%s\n',['%Created on: ',datestr(clock)]);
fprintf(fid,'%s\n',['%Original file: ',mFileNameInput]);
fclose(fid);

%% Saving PARFOR M-file
fid=fopen(mFileName_PARFOR,'wt');
for lineIndexSave=1:size(T,1); 
    if logicParFor(lineIndexSave)
        fprintf(fid,'%s\n',T_PARFOR{lineIndexSave});
    end
end
%Add last lines for information
fprintf(fid,'%s\n','%% --- PARSEFOR parfor-loop compatible file ---');
fprintf(fid,'%s\n',['%Created on: ',datestr(clock)]);
fprintf(fid,'%s\n',['%The original file can be found here: ',mFileNameInput]);
fclose(fid);

