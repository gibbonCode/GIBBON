function varargout=set_mat_par_FEBIO(febXML,saveName,mat_cell)

% function FEB_XML=set_mat_par_FEBIO(filename,savename,mat_cell)
% ------------------------------------------------------------------------
%
% Example:
% mat_cell={};
% mat_struct.id=1;
% mat_struct.par_names={'density'};
% mat_struct.par_values={pi};
% mat_cell{1}=mat_struct;
%
% mat_struct.id=3;
% mat_struct.par_names={'c1','k'};
% mat_struct.par_values={[66,5,5],77};
% mat_cell{2}=mat_struct;
%
% Leads to:
%
% <Material>
%       <material id="1" name="rigid_body" type="rigid body">
%          <density>3.1415927e+000</density>
%          <center_of_mass>0,0,0</center_of_mass>
%       </material>
%       <material id="3" name="fibre_component" type="trans iso Mooney-Rivlin">
%          <c1>6.6000000e+001,5.0000000e+000,5.0000000e+000</c1>
%          <c2>0</c2>
%          <c3>0.01</c3>
%          <c4>1</c4>
%          <c5>0</c5>
%          <k>7.7000000e+001</k>
%       </material>
% </Material>
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/10/10
%------------------------------------------------------------------------

%%
if exist(febXML,'file')==2
    
    FEB_XML = xmlread(febXML); %Read in XML object
    MAT_FEB_XML = FEB_XML.getElementsByTagName('Material'); %Material level
    mat_FEB_XML = MAT_FEB_XML.item(0).getElementsByTagName('material'); %material level
    numMats=mat_FEB_XML.getLength; %Number of material entries found
    
    %Retrieving material ID's
    mat_ID=ones(1,numMats); %Initialise ID's
    for q=1:1:numMats
        mat_ID(q)=str2double(mat_FEB_XML.item(q-1).getAttribute('id').toCharArray()'); %Read in ID
    end
    
    %Altering material entries
    for q=1:1:numel(mat_cell)
        mat_id=mat_cell{q}.id; %The current material id
        mat_id_ind=find(mat_ID==mat_id)-1; %material index in XML
        
        for qp=1:1:numel(mat_cell{q}.par_names)
            par_data_text=sprintf('%6.7e,',mat_cell{q}.par_values{qp}); %formatted as e.g. 4.4408921e-016
            par_data_text=par_data_text(1:end-1); %because an extra comma is added to the end
            if iscell(mat_cell{q}.par_names{qp}); %uncoupled solid mixture type entry
                %here it is assumed that mat_cell{i}.par_names{j}{1:3} refers to
                %solid, type and parameter.
                solid_element=mat_FEB_XML.item(mat_id_ind).getElementsByTagName(mat_cell{q}.par_names{qp}{1});
                for k=0:1:solid_element.getLength-1
                    if strcmp(solid_element.item(k).getAttribute('type').toCharArray',mat_cell{q}.par_names{qp}{2})
                        solid_element.item(k).getElementsByTagName(mat_cell{q}.par_names{qp}{3}).item(0).getFirstChild.setData(par_data_text);%Setting parameter
                        break
                    end
                end
            else
                mat_FEB_XML.item(mat_id_ind).getElementsByTagName(mat_cell{q}.par_names{qp}).item(0).getFirstChild.setData(par_data_text); %Setting parameter
            end
            
        end
    end
    
    %Saving XML file
    % xmlwrite(savename,FEB_XML);
    if ~isempty(saveName)
        write_XML_no_extra_lines(saveName,FEB_XML); %Only save if a save name is given
    end
    
    %Optional output
    varargout{1}=FEB_XML;
    
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
