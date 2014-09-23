clear all; close all; clc;

testCase=0; %% //parEx
switch testCase %% //parEx
    case 1 %% //parEx
        
        %% //parOn
        try
            matlabpool close;
        catch exception
            disp('--> No existing pools to close');
        end
        
        %Define profile and number of workers
        defaultProfile = parallel.defaultClusterProfile;
        myCluster = parcluster(defaultProfile);
        myCluster.NumWorkers=2;
        
        %Open matlab pool
        matlabpool(myCluster, 'open');
        
        %% //parOn
        n=1e6;
        tic;
        A=zeros(n,1);
        
        %% //parOn
        parfor q=1:n;
            
            %% //parOn
            
            %% //parOff
        %for q=1:n;
            
            %% //parOff
            A(q)=q;
        end
        t=toc
        %% //parOn
        try
            matlabpool close;
        catch exception
            disp('--> No existing pools to close');
        end
    case 0 %% //parEx
        
end %% //parEx
%% //parOn

mFileNameInput=[mfilename('fullpath'),'.m']; %% //parEx
[mFileName_FOR,mFileName_PARFOR]=parsefor(mFileNameInput,'d'); %% //parEx
run(mFileName_FOR); %% //parEx

%% 
%
% <<gibbVerySmall.gif>>
% 
% GIBBON 
% 
% Kevin M. Moerman (kevinmoerman@hotmail.com)