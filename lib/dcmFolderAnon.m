function dcmFolderAnon(pathName,varargin)

% function dcmFolderAnon(pathName,varargin)
% ------------------------------------------------------------------------
%
% This function anonomizes the dicom files in the folder pathName using the
% dicomanon function. The first input is the folder path name, the
% following inputs are possible inputs of the dicomanon function.
%
% See also: dicomanon
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/02/09: Created
%------------------------------------------------------------------------

%%

files = dir(fullfile(pathName,'*.dcm'));
files={files(1:end).name};
files=sort(files(:));

%% SEt dicom dictionary
try
    
    fName=fullfile(pathName,files{1});
    
    %First import using factory settings
    dicomdict('factory');
    dcmInfo_full=dicominfo(fName);
    
    if ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'philips'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        dicomdict('set','gibbon_dict.txt');
        disp('DICOM dictionary set to: gibbon_dict.txt');        
    elseif ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'pms'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        dicomdict('set','gibbon_dict.txt');
        disp('DICOM dictionary set to: gibbon_dict.txt');        
    elseif ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'siemens'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        dicomdict('set','dicom-dict-siemens.txt');
        disp('DICOM dictionary set to: dicom-dict-siemens.txt');        
    elseif ~isempty(strfind(lower(dcmInfo_full.Manufacturer),'GE MEDICAL SYSTEMS'))
        disp(['Detected ',dcmInfo_full.Manufacturer,' files']);
        warning('No settings for this vendor using DICOM dictionary factory settings');        
    else
        dicomdict('factory');
        warning('Unknown vendor, using DICOM dictionary factory settings');        
    end
    %Test info import with new dictionary (if this fails we resume with
    % factory
    dcmInfo_full=dicominfo(fName);
catch %e.g. if the manufacturer field is missing
    dicomdict('factory');
    warning('Unknown vendor, using DICOM dictionary factory settings');    
end
  
%%
hw = waitbar(0,'Anonomizing DICOM info...');
try
    for q=1:1:numel(files)
        fileName=fullfile(pathName,files{q});
        dicomanon(fileName,fileName,varargin{:});
        waitbar(q/numel(files),hw,['Anonomizing DICOM info...',num2str(round(100.*q/numel(files))),'%']);
    end
catch ME
    close(hw);
    rethrow(ME);
    dicomdict('factory');
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
