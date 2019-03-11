function [varargout]=febioDoc(varargin)

% function [hFig]=febioDoc(docType)
%------------------------------------------------------------------------
% View febio documentation PDF files using a MATLAB browser. Use docType
% 1,'UM','um' for the user manual, or use 2,'TM','tm' for the theory
% manual.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/02/02 : Created
%------------------------------------------------------------------------

%% Parse input
switch nargin
    case 0
        docType=1;
    case 1
        docType=varargin{1};
end

switch docType
    case {'UM','um',1} %User manual
        fileTarget='FEBio_um';        
    case {'TM','tm',2} %Theory manual
        fileTarget='FEBio_tm';        
end

%% Get pdf file location

febioPath=getFEBioPath;
documentationPath=fullfile(fileparts(fileparts(febioPath)),'doc');

%Get list of pdf files in folder
files = dir(fullfile(documentationPath,'*.pdf'));
files={files(1:end).name};
files=sort(files(:));
NumberOfFiles=numel(files);

%% View pdf file

if NumberOfFiles>1
    for q=1:1:NumberOfFiles
        if gcontains(files{q},fileTarget) %If it contains the search pattern
            fileName=fullfile(documentationPath,files{q}); %The full pdf file name
            [hFig]=pdfView(fileName); %View pdf file using pdfView            
            if nargout==1
                varargout{1}=hFig;
            end
            break
        end
    end
else
    error('No documentation file found')
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
