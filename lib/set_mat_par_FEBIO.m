function FEB_XML=set_mat_par_FEBIO(filename,savename,mat_cell)

% Alters material parameters in XML formatted file
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
% 25/01/2012, Kevin Mattheus Moerman
% kevinmoerman@gmail.com

%%
FEB_XML = xmlread(filename);
MAT_FEB_XML = FEB_XML.getElementsByTagName('Material'); %Material level
mat_FEB_XML = MAT_FEB_XML.item(0).getElementsByTagName('material'); %material level
no_mats=mat_FEB_XML.getLength;

%Retrieving material ID's
mat_ID=ones(1,no_mats);
for i=0:1:no_mats-1
    mat_ID(i+1)=str2double(mat_FEB_XML.item(i).getAttribute('id').toCharArray()');
end

for i=1:1:numel(mat_cell)
   mat_id=mat_cell{i}.id; %material id
   mat_ind=find(mat_ID==mat_id)-1; %material index in XML
   for j=1:1:numel(mat_cell{i}.par_names)       
       par_data_text=sprintf('%6.7e,',mat_cell{i}.par_values{j}); %formatted as e.g. 4.4408921e-016       
       par_data_text=par_data_text(1:end-1); %because an extra comma is added to the end
       if iscell(mat_cell{i}.par_names{j}); %uncoupled solid mixture entry
           %here it is assumed that mat_cell{i}.par_names{j}{1:3} refers to
           %solid, type and parameter.
           solid_element=mat_FEB_XML.item(mat_ind).getElementsByTagName(mat_cell{i}.par_names{j}{1});
           for k=0:1:solid_element.getLength-1
               if strcmp(solid_element.item(k).getAttribute('type').toCharArray',mat_cell{i}.par_names{j}{2})
                   solid_element.item(k).getElementsByTagName(mat_cell{i}.par_names{j}{3}).item(0).getFirstChild.setData(par_data_text);%Setting parameter
                   break
               end
           end
       else 
           mat_FEB_XML.item(mat_ind).getElementsByTagName(mat_cell{i}.par_names{j}).item(0).getFirstChild.setData(par_data_text); %Setting parameter
       end
       
   end
end

%Saving XML file
% xmlwrite(savename,FEB_XML);
write_XML_no_extra_lines(savename,FEB_XML);

end