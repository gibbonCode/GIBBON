%% DEMO_import_FEB_export_INP
% Below is a demonstration of how import a FEB file and subsequently export
% the geometry into an INP file. 

%%

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=0.5;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=50; 

%% Importing .feb file

%Set main folders
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName_FEB=fullfile(defaultFolder,'data','FEB'); %Where to load the FEB file
pathName_INP=fullfile(defaultFolder,'data','INP'); %Where to export the INP file

febFileNamePart='example_HEX_QUAD.feb';
febFileName=fullfile(pathName_FEB,febFileNamePart);
[febXML,nodeStruct,elementCell]=import_FEB(febFileName);

%% Plotting model

% Plotting the example model surfaces
hf1=cFigure;
title('Visualizing model','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

uniqueMaterialIndices=[];
for q=1:1:numel(elementCell)
    uniqueMaterialIndices=unique([uniqueMaterialIndices(:); elementCell{q}.E_mat(:)]);
    switch elementCell{q}.E_type
        case {'tri3', 'quad4'}
            F=elementCell{q}.E;
            V=nodeStruct.N;
            C=elementCell{q}.E_mat;            
       case {'hex8', 'tet4'}
            [F,C]=element2faces(elementCell{q}.E,elementCell{q}.E_mat); %Creates faces and colors (e.g. stress) for patch based plotting
    end
    hp=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','flat','Cdata',C,'FaceAlpha',0.8);   
end

colormap(jet(numel(uniqueMaterialIndices))); hc=colorbar; caxis([min(uniqueMaterialIndices)-0.5 max(uniqueMaterialIndices)+0.5]);
axis equal; view(3); axis tight;  grid on; set(gca,'FontSize',fontSize);
camlight('headlight');
drawnow;

%% EXPORTING INP FILES FOR EACH ELEMENT TYPE

%You can change this example to do this for material type instead. Just use
%the material indices to select the elements from the lists. However the
%export_INP function can only handle 1 element type at a time at the moment

for q=1:1:numel(elementCell)
    
    inpFileNamepart=[febFileNamePart(1:end-4),'_',num2str(q),'.inp']; %filename for inp file
    inpFileName=fullfile(pathName_INP,inpFileNamepart);

    elementStruct=elementCell{q};    
    %Setting appropriate element type line for ABAQUS. CHECK THESE!
    switch elementStruct.E_type
        case 'tri3'
            elementStruct.E_type='*ELEMENT, TYPE=STRI3, ELSET=PART-DEFAULT_1_EB1';
        case 'quad4'
            elementStruct.E_type='*ELEMENT, TYPE=S4R, ELSET=PART-DEFAULT_1_EB1';
        case 'tet4'
            elementStruct.E_type='*ELEMENT, TYPE=C3D4, ELSET=PART-DEFAULT_1_EB1';            
        case 'hex8'            
            elementStruct.E_type='*ELEMENT, TYPE=C3D8R, ELSET=PART-DEFAULT_1_EB1';
    end
    export_INP(elementStruct,nodeStruct,inpFileName);
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
