function [Fn,Vn]=ggremesh(varargin)

% function [Fn,Vn]=ggremesh(F,V,optionStruct)
%-------------------------------------------------------------------------
% This function uses the external library Geogram to remesh the input
% triangulation defined by the faces F and the vertices V. In particular
% the code "vorpalite" is used. An additional option structure may be
% provided where users can set particular parameters for Geogram. 
% 
% Below the options and defaults are provided: 
% optionStruct.nb_pts=size(V,1); %number of points
% optionStruct.anisotropy=0; %Use anisotropy (~=0) to capture geometry or favour isotropic triangles (=0)
% optionStruct.pre.max_hole_area=100; %Max hole area for pre-processing step
% optionStruct.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
% optionStruct.post.max_hole_area=100; %Max hole area for post-processing step
% optionStruct.post.max_hole_edges=0; %Max number of hole edges for post-processing step
% optionStruct.disp_on=1; %Turn on/off displaying of Geogram text
%
% Instead of nb_pts users can also specify a pointSpacing to be used
% instead of nb_pts. This is not a Geogram feature but a GIBBON option
% which is translated to the number of points for Geogram remeshing. This
% is and example for a desired point spacing of 4:  
% optionStruct.pointSpacing=4
%
%
% Geogram GitHub repository: https://github.com/BrunoLevy/geogram
% 
%
% Lévy B., Bonneel N. (2013) Variational Anisotropic Surface Meshing with
% Voronoi Parallel Linear Enumeration. In: Jiao X., Weill JC. (eds)
% Proceedings of the 21st International Meshing Roundtable. Springer,
% Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-33573-0_21 
% 
% See also: 
% http://alice.loria.fr/publications/papers/2012/Vorpaline_IMR/vorpaline.pdf
% https://www.ljll.math.upmc.fr/hecht/ftp/ff++days/2013/BrunoLevy.pdf
%
%
% 2020/12/04 KMM Created
% 2021/07/16 KMM Switched default to disp_on=0
%-------------------------------------------------------------------------


%% Parse input

switch nargin    
    case 2
        F=varargin{1};
        V=varargin{2};
        optionStruct=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        optionStruct=varargin{3};
end

%Set default structure
defaultOptionStruct.nb_pts=size(V,1); %resample with same number of points
defaultOptionStruct.anisotropy=0; %Use anisotropy (~=0) to capture geometry or favour isotropic triangles (=0)
defaultOptionStruct.pre.max_hole_area=100; %Max hole area for pre-processing step
defaultOptionStruct.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
defaultOptionStruct.post.max_hole_area=100; %Max hole area for post-processing step
defaultOptionStruct.post.max_hole_edges=0; %Max number of hole edges for post-processing step
defaultOptionStruct.disp_on=0; %Turn on/off displaying of Geogram text

%Complement input with default if missing
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

%Check for point spacing field
if isfield(optionStruct,'pointSpacing')
    optionStruct.nb_pts=spacing2numVertices(F,V,optionStruct.pointSpacing); 
    optionStruct=rmfield(optionStruct,'pointSpacing');
end

%Get display setting
disp_on=optionStruct.disp_on;
optionStruct=rmfield(optionStruct,'disp_on');

%%

%Setting pathnames
pathNameLib=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','geogram');
pathNameTempFiles=fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','temp');

%Create temporary file name
inputFileName=fullfile(pathNameTempFiles,'temp.obj');
outputFileName=fullfile(pathNameTempFiles,'temp_out.obj');

%Create runName for binary
compString=computer; 
switch compString
    case 'PCWIN64' %Windows 64-bit
        pathNameTetGenFile=fullfile(pathNameLib,'win64','bin');
        runName=fullfile(pathNameTetGenFile,'vorpalite.exe');
    case 'GLNXA64' %Linux 64-bit       
        pathNameTetGenFile=fullfile(pathNameLib,'lin64','bin');
        runName=fullfile(pathNameTetGenFile,'vorpalite');
    case 'MACI64'  %MAC 64-bit      
        pathNameTetGenFile=fullfile(pathNameLib,'mac64','bin');
        runName=fullfile(pathNameTetGenFile,'vorpalite');
    otherwise
        error('Your platform does not seem to be supported');
end

%% 

if disp_on==1
    startString=['------>  Geogram/vorpalite for resmeshing  <------ ',datestr(now)];
    stringLength=numel(startString);
    lineSep=repmat('%',1,stringLength); %Stripe of % symbols
    disp(' '); %Empty line
    disp(lineSep);
    disp(startString);
end

%%  Export to obj

if disp_on==1
    dispMessage('# Export mesh input file.',stringLength);
end

patch2obj(inputFileName,F,V);

%% Start Geogram -> vorpalite

if disp_on==1
    dispMessage('# Run Geomgram/vorpalite.',stringLength);
end

%Compose basic run string
runString=['"',runName,'" "',inputFileName,'" "',outputFileName,'"']; 

%Grow run string with additional options
fieldNameSet=fieldnames(optionStruct);
for q=1:1:numel(fieldNameSet)
    fieldNameNow=fieldNameSet{q};    
    structVar=optionStruct.(fieldNameNow);
    
    if isstruct(structVar)        
        preFix=fieldNameNow; 
        fieldNameSetPre=fieldnames(structVar);
        for qs=1:1:numel(fieldNameSetPre)
            fieldNameNow=fieldNameSetPre{qs};
            structVarSub=structVar.(fieldNameNow);
            [runString]=appendOption(runString,preFix,fieldNameNow,structVarSub);
        end
    else
        [runString]=appendOption(runString,[],fieldNameNow,structVar);
    end
end

if disp_on==1
    [runStatus,runOut]=system(runString,'-echo');
else
    [runStatus,runOut]=system(runString);
end

%% Import mesh

if disp_on==1
    dispMessage('# Importing remeshed geometry.',stringLength);
end

try
    objImportOptionStruct.fullMode=0;
    objStruct=import_obj(outputFileName,objImportOptionStruct);
    Fn=objStruct.F;
    Vn=objStruct.V;
catch ME
    warning('import of OBJ file not successful. Set disp_on=1 to see any Geogram error messages.');    
    rethrow(ME);
end

%%

if disp_on==1
    dispMessage('# Removing temporary files.',stringLength);
end

%Delete the file
delete(inputFileName); 
delete(outputFileName); 

%%
if disp_on==1
    dispMessage('# Done!',stringLength);
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

function [runString]=appendOption(runString,preFix,varName,varValue)

if isnumeric(varValue) || islogical(varValue)
    varValue=sprintf('%.16g',varValue);
end

if ~isempty(preFix)
    runString=[runString,' ',preFix,':',varName,'=',varValue];
else
    runString=[runString,' ',varName,'=',varValue];
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
