function patch2obj(varargin)

% function patch2obj(objFileName,F,V,C,cMap,cLim,mtlStruct)
% ------------------------------------------------------------------------
% This function exports the patch data defined by the faces (F), vertices
% (V) and the color data (C) to the OBJ (Wavefront .obj) format. The
% function generates a .obj file, a .mtl file, and a .jpg image file. The
% .obj file contains the geometry information and texture/color coordinates
% to use. The .mtl file contains the material information and refers to the
% image to use to look up the colors based on the texture coordinates in
% the .obj file.
% The color data C should ideally define either the vertex or face colors
% in the form of an nx1 array. If face colors are provided these are
% re-sampled (averaged) to vertex colors which is the required format for
% OBJ files. Colors are obtained from the input color map as well as the
% color limits. The input structure mtlStruct defines the MTL file
% components. With the default entries:
%
% mtlStruct.Ka=[1 1 1]; %Ambient color
% mtlStruct.Kd=[1 1 1]; %Diffuse color
% mtlStruct.Ks=[0 0 0]; %Specular color, black=off
% mtlStruct.Ns=0; %Specular component [0-1000]
% mtlStruct.Ni=1.45; %Optical density/index of refraction
% mtlStruct.d=1; %"dissolved"/transparancy [0-1]
% mtlStruct.Tr=0; %1 - d, used instead of d by some software
% mtlStruct.illum=1; %Illumination model
%
% Illumination models:
% 0. Color on and Ambient off
% 1. Color on and Ambient on
% 2. Highlight on
% 3. Reflection on and Ray trace on
% 4. Transparency: Glass on, Reflection: Ray trace on
% 5. Reflection: Fresnel on and Ray trace on
% 6. Transparency: Refraction on, Reflection: Fresnel off and Ray trace on
% 7. Transparency: Refraction on, Reflection: Fresnel on and Ray trace on
% 8. Reflection on and Ray trace off
% 9. Transparency: Glass on, Reflection: Ray trace off
% 10. Casts shadows onto invisible surfaces
%
%
% For more information on the OBJ file format see:
% https://en.wikipedia.org/wiki/Wavefront_.obj_file
% http://paulbourke.net/dataformats/obj/minobj.html
%
% Kevin Mattheus Moerman
%
% Change log:
% 2020/05/26 Created
% ------------------------------------------------------------------------

%% parse input

switch nargin
    case 3
        objFileName=varargin{1};
        F=varargin{2};
        V=varargin{3};
        C=[];
        cMap=gray(250);
        cLim=[0 1];
        mtlStruct=[];
    case 4
        objFileName=varargin{1};
        F=varargin{2};
        V=varargin{3};
        C=varargin{4};
        cMap=viridis(250);
        cLim=[];
        mtlStruct=[];
    case 5
        objFileName=varargin{1};
        F=varargin{2};
        V=varargin{3};
        C=varargin{4};
        cMap=varargin{5};
        cLim=[];
        mtlStruct=[];
    case 6
        objFileName=varargin{1};
        F=varargin{2};
        V=varargin{3};
        C=varargin{4};
        cMap=varargin{5};
        cLim=varargin{6};
        mtlStruct=[];
    case 7
        objFileName=varargin{1};
        F=varargin{2};
        V=varargin{3};
        C=varargin{4};
        cMap=varargin{5};
        cLim=varargin{6};
        mtlStruct=varargin{7};
end

%Parse mtlStruct (make complete is fields are missing)
mtlStructDefault.Ka=[1 1 1];
mtlStructDefault.Kd=[1 1 1];
mtlStructDefault.Ks=[0 0 0];
mtlStructDefault.illum=1;
mtlStructDefault.Ns=0;
mtlStructDefault.Ni=1.45;
mtlStruct=structComplete(mtlStruct,mtlStructDefault,1); %Complement provided with default if missing or empty

% Define Tr and d and add if missing
if ~isfield(mtlStruct,'d') && isfield(mtlStruct,'Tr')
    mtlStruct.d=1-mtlStruct.Tr;
end

if ~isfield(mtlStruct,'Tr') && isfield(mtlStruct,'d')
    mtlStruct.Tr=1-mtlStruct.d;
end

if ~isfield(mtlStruct,'Tr') && ~isfield(mtlStruct,'d')
    mtlStruct.d=1;
    mtlStruct.Tr=0;
end

if isempty(C)
    colorTextureOutput=false;
else
    colorTextureOutput=true;
    
    %Check if colors are vertex colors
    if size(C,1)==size(F,1) %if face colors are provided
        C=faceToVertexMeasure(F,V,C); %Convert to vertex data
    end
end

%If coordinates are 2D add zero valued Z-coordinates
if size(V,2)==2
    V(:,3)=0;
end

%Settings
formatDouble='%6.7e';
formatInteger='%d';

%% Open file for writing
fid = fopen(objFileName,'w');
fprintf(fid,'%s\n',['# Created using GIBBON, ',datestr(now)]);

%%
if colorTextureOutput
    
    %% Set up folder and file name variables
    % Get path names
    [pathName,fileName,~]=fileparts(objFileName);
    
    %Create mtl file names
    mtlName=[fileName,'.mtl'];
    mtlFileName=fullfile(pathName,mtlName);
    
    %Create jpg file name
    jpgName=[fileName,'.png'];
    jpgFileName=fullfile(pathName,jpgName);
    
    %% Save mtl file
    
    T={'newmtl material0'};
    fieldNames=fieldnames(mtlStruct);
    for q=1:1:numel(fieldNames)
        v=mtlStruct.(fieldNames{q});
        if all(isrounded(v))
            t_form=repmat(['%d',' '],1,size(v,2));
        else
            t_form=repmat(['%f',' '],1,size(v,2));
        end
        t_form=t_form(1:end-1); %Crop trailing delimiter
        t=sprintf(t_form,v);
        T{q+1,1}=[fieldNames{q},' ',t];
    end
    T{end+1,1}=['map_Kd ',jpgName];
    
    cell2txtfile(mtlFileName,T);
    
    %% Save texture image
    
    m=1; %required pixel replication factor (should be uneven)
    if size(C,2)==1
        [~,cmapIndex]=cmaperise(C,cMap,cLim); %texture linear indices
        ind=repmat(1:1:size(cMap,1),m,1);
        cMap=cMap(ind,:);
        cmapIndex=cmapIndex*m-floor(m/2);
        cmapImage=uint16(permute(cMap*65535,[1 3 2]));
    elseif size(C,2)==3
        %WIP interpolation inaccuraries appear in visualized results if colors
        %are not properly interpolatable based on texture coordinates.
        cmapImage=uint16(C*65535);
        [cmapImage,~,cmapIndex]=unique(cmapImage,'rows');
        ind=repmat(1:1:size(cmapImage,1),m,1);
        cmapImage=cmapImage(ind,:);
        cmapIndex=cmapIndex*m-floor(m/2);
        cmapImage=permute(cmapImage,[1 3 2]);
    else
        error('Color data should be nx1 or nx3, where n is either the number of faces or the number of vertices');
    end
    
    imwrite(cmapImage,jpgFileName,'png','BitDepth',16);
    
    %% Create texture coordinates and indices
    pixelSize=1/size(cmapImage,1);
    T=[0.5*ones(size(cmapImage,1),1) linspace(1-pixelSize/2,pixelSize/2,size(cmapImage,1))'];
    CF_ind=cmapIndex(F); %Texture indices for all face vertices
    
    %% Add object line
    % fprintf(fid,'%s\n','o patchObject');
    
    %% Add mtl file description
    fprintf(fid,'%s\n','# ---------------------------------------------------');
    fprintf(fid,'%s\n','# Specify MTL file to use');
    fprintf(fid,['%s\n','mtllib ',mtlName]);
    
    %% Add vertices
    [fid]=addVertices(fid,V,formatDouble);
    
    %% Add normals
    [fid]=addNormals(fid,F,V,formatDouble);
    
    %% Add texture coordinates
    
    %Create text
    t_form=repmat([formatDouble,' '],1,size(T,2));
    t_form=t_form(1:end-1); %Crop trailing delimiter
    t_form=['vt ',t_form,'\n']; %Append end of line character
    t=sprintf(t_form,T');      %Convert to string
    t=t(1:end-1); %Take away last end of line character
    
    %Write text to file
    fprintf(fid,'%s\n','# ---------------------------------------------------');
    fprintf(fid,'%s\n','# Texture coordinates');
    fprintf(fid,'%s\r\n',t);
    
    %% Add material description for face set
    fprintf(fid,'%s\n','# ---------------------------------------------------');
    fprintf(fid,'%s\n','# Specify material from MTL file for face set below');
    fprintf(fid,'%s\n','usemtl material0');
    
    %% Add faces
    [fid]=addFaces(fid,F,formatInteger,CF_ind);
    
else %Simple minimal OBJ file for geometry only
    
    %% Add vertices
    [fid]=addVertices(fid,V,formatDouble);
    
    %% Add normals
    %     [fid]=addNormals(fid,F,V,formatDouble);
    
    %% Add faces
    [fid]=addFaces(fid,F,formatInteger,[]);
    
end

%% Close file
fclose(fid);

end

function [fid]=addVertices(fid,V,formatDouble)

%Create text
t_form=repmat([formatDouble,' '],1,size(V,2));
t_form=t_form(1:end-1); %Crop trailing delimiter
t_form=['v ',t_form,'\n']; %Append end of line character
t=sprintf(t_form,V');      %Convert to string
t=t(1:end-1); %Take away last end of line character

%Write text to file
fprintf(fid,'%s\n','# ---------------------------------------------------');
fprintf(fid,'%s\n','# Vertices');
fprintf(fid,'%s\r\n',t);

end

function [fid]=addNormals(fid,F,V,formatDouble)
[~,~,N]=patchNormal(F,V); %Vertex normals

%Create text
t_form=repmat([formatDouble,' '],1,size(N,2));
t_form=t_form(1:end-1); %Crop trailing delimiter
t_form=['vn ',t_form,'\n']; %Append end of line character
t=sprintf(t_form,N');      %Convert to string
t=t(1:end-1); %Take away last end of line character

%Write text to file
fprintf(fid,'%s\n','# ---------------------------------------------------');
fprintf(fid,'%s\n','# Vertex normal vectors');
fprintf(fid,'%s\r\n',t);
end

function [fid]=addFaces(fid,F,formatInteger,CF_ind)
% Form:
% f v1 v2 v3 # A face with no texture indices
% f f v1/vt1 v2/vt2 v3/vt3 # A face with texture indices
% f f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 # A face with texture and normal indices

if isempty(CF_ind)
    FF=F;
    t_form=repmat([formatInteger,' '],1,size(F,2));
    t_form=t_form(1:end-1); %Crop trailing delimiter
    
    %     else %Normal direction indices
    %         nRep=2;
    %         d='//';
    %
    %         FF=zeros(size(F,1),nRep*size(F,2));
    %
    %         columnIndFaces=1:2:size(FF,2);
    %         FF(:,columnIndFaces)=F;
    %
    %         columnIndNormal=2:2:size(FF,2);
    %         FF(:,columnIndNormal)=F;
    %
    %         t_form_p=repmat([formatInteger,d],1,nRep);
    %         t_form_p=t_form_p(1:end-numel(d)); %Crop trailing delimiter
    %
    %         t_form=repmat([t_form_p,' '],1,size(F,2));
    %         t_form=t_form(1:end-1); %Crop trailing delimiter
    
else
    nRep=3;
    d='/';
    
    FF=zeros(size(F,1),nRep*size(F,2));
    
    columnIndFaces=1:3:size(FF,2);
    FF(:,columnIndFaces)=F;
    
    columnIndTexture=2:3:size(FF,2);
    FF(:,columnIndTexture)=CF_ind;
    
    columnIndNormal=3:3:size(FF,2);
    FF(:,columnIndNormal)=F;
    
    t_form_p=repmat([formatInteger,d],1,nRep);
    t_form_p=t_form_p(1:end-numel(d)); %Crop trailing delimiter
    
    t_form=repmat([t_form_p,' '],1,size(F,2));
    t_form=t_form(1:end-1); %Crop trailing delimiter
end

t_form=['f ',t_form,'\n']; %Append end of line character
t=sprintf(t_form,FF');      %Convert to string
t=t(1:end-1); %Take away last end of line character

%Write text to file
fprintf(fid,'%s\n','# ---------------------------------------------------');
fprintf(fid,'%s\n','# Faces');
fprintf(fid,'%s\r\n',t);

end