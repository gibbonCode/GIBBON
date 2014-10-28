function docNode=addMaterialLevel_FEB_v1p2(docNode,FEB_struct)


%%
febio_spec = docNode.getDocumentElement;

ElementAddNode = docNode.createElement('Material');
febio_spec.appendChild(ElementAddNode);
MATNode = docNode.getElementsByTagName('Material').item(0);

% Adding material fields
disp('Adding Material level')
A=[FEB_struct.Geometry.ElementMat{:}];
mat_ind=unique(A);
parent_node =MATNode;

if FEB_struct.disp_opt==1;
    hw = waitbar(0,'Adding material level...');
end
for i=1:1:numel(mat_ind)
    
    %Adding material node
    material_node = docNode.createElement('material'); %create material entry
    material_node = parent_node.appendChild(material_node); %add material entry
    
    %Checking material type
    Mat_i=FEB_struct.Materials{i};
    mat_type=Mat_i.type;
    switch mat_type
        case {'uncoupled solid mixture','solid mixture'}
            
            %%
            
            %Mat ID
            attr = docNode.createAttribute('id'); %Create attribute
            attr.setNodeValue(num2str(mat_ind(i))); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Name
            attr = docNode.createAttribute('name'); %Create attribute
            attr.setNodeValue(['mat_',num2str(i)]); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Type
            attr = docNode.createAttribute('type'); %Create attribute
            attr.setNodeValue(mat_type); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat properties and attributes
            if isfield(Mat_i,'props')
                mat_props=Mat_i.props;
                mat_prop_vals=Mat_i.vals;
                
                %Access property attributes
                if isfield(Mat_i,'PropAttrName')
                    mat_prop_attr_name=Mat_i.PropAttrName;
                    mat_prop_attr_val=Mat_i.PropAttrVal;
                else
                    mat_prop_attr_name=[];
                    mat_prop_attr_val=[];
                end
                
                %Material properties
                for q=1:1:numel(mat_props)
                    mat_prop=mat_props{q};
                    mat_prop_val=mat_prop_vals{q};
                    prop_node = docNode.createElement(mat_prop); %create entry
                    prop_node = material_node.appendChild(prop_node); %add entry
                    t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                    
                    %Add potential property attributes (e.g. property load curves)
                    if ~isempty(mat_prop_attr_name)
                        if q<=numel(mat_prop_attr_name)
                            currentPropAttrName=mat_prop_attr_name{q};
                            if ~isempty(currentPropAttrName)
                                currentPropAttrVal=mat_prop_attr_val{q};
                                
                                t_form=repmat('%6.7e, ',1,size(currentPropAttrVal,2));
                                t_form=t_form(1:end-2);
                                
                                attr = docNode.createAttribute(currentPropAttrName); %Create attribute
                                attr.setNodeValue(sprintf(t_form,currentPropAttrVal)); %Set text
                                prop_node.setAttributeNode(attr); %Add attribute
                            end
                        end
                    end
                    
                end
            end
            
            %Mat aniso_type
            if isfield(Mat_i,'aniso_type')
                mat_aniso_type=Mat_i.aniso_type;
                switch mat_aniso_type
                    case 'none'
                        
                    case 'fiber'
                        prop_node = docNode.createElement('fiber'); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        attr = docNode.createAttribute('type'); %Create attribute
                        attr.setNodeValue('user'); %Set text
                        prop_node.setAttributeNode(attr); %Add attribute
                    case 'mat_axis'
                        prop_node = docNode.createElement('mat_axis'); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        attr = docNode.createAttribute('type'); %Create attribute
                        attr.setNodeValue('vector'); %Set text
                        prop_node.setAttributeNode(attr); %Add attribute
                end
            end
            %%
            
            for j=1:1:numel(FEB_struct.Materials{i}.Mats)
                
                solid_node = docNode.createElement('solid'); %create entry
                solid_node = material_node.appendChild(solid_node); %add entry
                
                %Mat Type
                attr = docNode.createAttribute('type'); %Create attribute
                attr.setNodeValue(Mat_i.Mats{j}.type); %Set text
                solid_node.setAttributeNode(attr); %Add attribute
                
                %Mat ID
                attr = docNode.createAttribute('id'); %Create attribute
                attr.setNodeValue(num2str(j)); %Set text
                solid_node.setAttributeNode(attr); %Add attribute
                
                %%
                %                 mat_props=Mat_i.Mats{j}.props;
                %                 mat_prop_vals=Mat_i.Mats{j}.vals;
                %                 %Material properties
                %                 for q=1:1:numel(mat_props)
                %                     mat_prop=mat_props{q};
                %                     mat_prop_val=mat_prop_vals{q};
                %                     prop_node = docNode.createElement(mat_prop); %create entry
                %                     prop_node = solid_node.appendChild(prop_node); %add entry
                %                     t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                %                     prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                %                 end
                
                %Mat properties and attributes
                if isfield(Mat_i.Mats{j},'props')
                    mat_props=Mat_i.Mats{j}.props;
                    mat_prop_vals=Mat_i.Mats{j}.vals;
                    
                    %Access property attributes
                    if isfield(Mat_i.Mats{j},'PropAttrName')
                        mat_prop_attr_name=Mat_i.Mats{j}.PropAttrName;
                        mat_prop_attr_val=Mat_i.Mats{j}.PropAttrVal;
                    else
                        mat_prop_attr_name=[];
                        mat_prop_attr_val=[];
                    end
                    
                    %Material properties
                    for q=1:1:numel(mat_props)
                        mat_prop=mat_props{q};
                        mat_prop_val=mat_prop_vals{q};
                        prop_node = docNode.createElement(mat_prop); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                        prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                        
                        %Add potential property attributes (e.g. property load curves)
                        if ~isempty(mat_prop_attr_name)
                            if q<=numel(mat_prop_attr_name)
                                currentPropAttrName=mat_prop_attr_name{q};
                                if ~isempty(currentPropAttrName)
                                    currentPropAttrVal=mat_prop_attr_val{q};
                                    
                                    t_form=repmat('%6.7e, ',1,size(currentPropAttrVal,2));
                                    t_form=t_form(1:end-2);
                                    
                                    attr = docNode.createAttribute(currentPropAttrName); %Create attribute
                                    attr.setNodeValue(sprintf(t_form,currentPropAttrVal)); %Set text
                                    prop_node.setAttributeNode(attr); %Add attribute
                                end
                            end
                        end
                    end
                end
                
                %%
                %Fiber spec
                mat_aniso_type=Mat_i.Mats{j}.aniso_type;
                switch mat_aniso_type
                    case 'none'
                        
                    case 'fiber'
                        prop_node = docNode.createElement('fiber'); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        attr = docNode.createAttribute('type'); %Create attribute
                        attr.setNodeValue('user'); %Set text
                        prop_node.setAttributeNode(attr); %Add attribute
                    case 'mat_axis'
                        prop_node = docNode.createElement('mat_axis'); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        attr = docNode.createAttribute('type'); %Create attribute
                        attr.setNodeValue('vector'); %Set text
                        prop_node.setAttributeNode(attr); %Add attribute
                end
                
            end
        case 'multigeneration'
            
            if ~isfield(FEB_struct.Materials{i},'start_times')
                
                %% BEFORE FEBIO 2.0
                %Mat ID
                attr = docNode.createAttribute('id'); %Create attribute
                attr.setNodeValue(num2str(mat_ind(i))); %Set text
                material_node.setAttributeNode(attr); %Add attribute
                
                %Mat Name
                attr = docNode.createAttribute('name'); %Create attribute
                attr.setNodeValue(['mat_',num2str(i)]); %Set text
                material_node.setAttributeNode(attr); %Add attribute
                
                %Mat Type
                attr = docNode.createAttribute('type'); %Create attribute
                attr.setNodeValue(mat_type); %Set text
                material_node.setAttributeNode(attr); %Add attribute
                
                for j=1:1:numel(FEB_struct.Materials{i}.Mats)
                    
                    solid_node = docNode.createElement('solid'); %create entry
                    solid_node = material_node.appendChild(solid_node); %add entry
                    
                    %Mat Type
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue(Mat_i.Mats{j}.type); %Set text
                    solid_node.setAttributeNode(attr); %Add attribute
                    
                    %Mat Generation
                    attr = docNode.createAttribute('gen'); %Create attribute
                    attr.setNodeValue(num2str(j)); %Set text
                    solid_node.setAttributeNode(attr); %Add attribute
                    
                    %Mat ID
                    attr = docNode.createAttribute('id'); %Create attribute
                    attr.setNodeValue(num2str(j)); %Set text
                    solid_node.setAttributeNode(attr); %Add attribute
                    
                    %%
                    %                     mat_props=Mat_i.Mats{j}.props;
                    %                     mat_prop_vals=Mat_i.Mats{j}.vals;
                    %                     %Material properties
                    %                     for q=1:1:numel(mat_props)
                    %                         mat_prop=mat_props{q};
                    %                         mat_prop_val=mat_prop_vals{q};
                    %                         prop_node = docNode.createElement(mat_prop); %create entry
                    %                         prop_node = solid_node.appendChild(prop_node); %add entry
                    %                         t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                    %                         prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                    %                     end
                    
                    %Mat properties and attributes
                    if isfield(Mat_i.Mats{j},'props')
                        mat_props=Mat_i.Mats{j}.props;
                        mat_prop_vals=Mat_i.Mats{j}.vals;
                        
                        %Access property attributes
                        if isfield(Mat_i.Mats{j},'PropAttrName')
                            mat_prop_attr_name=Mat_i.Mats{j}.PropAttrName;
                            mat_prop_attr_val=Mat_i.Mats{j}.PropAttrVal;
                        else
                            mat_prop_attr_name=[];
                            mat_prop_attr_val=[];
                        end
                        
                        %Material properties
                        for q=1:1:numel(mat_props)
                            mat_prop=mat_props{q};
                            mat_prop_val=mat_prop_vals{q};
                            prop_node = docNode.createElement(mat_prop); %create entry
                            prop_node = solid_node.appendChild(prop_node); %add entry
                            t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                            prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                            
                            %Add potential property attributes (e.g. property load curves)
                            if ~isempty(mat_prop_attr_name)
                                if q<=numel(mat_prop_attr_name)
                                    currentPropAttrName=mat_prop_attr_name{q};
                                    if ~isempty(currentPropAttrName)
                                        currentPropAttrVal=mat_prop_attr_val{q};
                                        
                                        t_form=repmat('%6.7e, ',1,size(currentPropAttrVal,2));
                                        t_form=t_form(1:end-2);
                                        
                                        attr = docNode.createAttribute(currentPropAttrName); %Create attribute
                                        attr.setNodeValue(sprintf(t_form,currentPropAttrVal)); %Set text
                                        prop_node.setAttributeNode(attr); %Add attribute
                                    end
                                end
                            end
                        end
                    end
                    
                    %%
                    
                    %Fiber spec
                    mat_aniso_type=Mat_i.Mats{j}.aniso_type;
                    switch mat_aniso_type
                        case 'none'
                            
                        case 'fiber'
                            prop_node = docNode.createElement('fiber'); %create entry
                            prop_node = solid_node.appendChild(prop_node); %add entry
                            attr = docNode.createAttribute('type'); %Create attribute
                            attr.setNodeValue('user'); %Set text
                            prop_node.setAttributeNode(attr); %Add attribute
                        case 'mat_axis'
                            prop_node = docNode.createElement('mat_axis'); %create entry
                            prop_node = solid_node.appendChild(prop_node); %add entry
                            attr = docNode.createAttribute('type'); %Create attribute
                            attr.setNodeValue('vector'); %Set text
                            prop_node.setAttributeNode(attr); %Add attribute
                    end
                    
                end
            else
                
                %% AS OF FEBIO 2.0
                
                %Mat ID
                attr = docNode.createAttribute('id'); %Create attribute
                attr.setNodeValue(num2str(mat_ind(i))); %Set text
                material_node.setAttributeNode(attr); %Add attribute
                
                %Mat Name
                attr = docNode.createAttribute('name'); %Create attribute
                attr.setNodeValue(['mat_',num2str(i)]); %Set text
                material_node.setAttributeNode(attr); %Add attribute
                
                %Mat Type
                attr = docNode.createAttribute('type'); %Create attribute
                attr.setNodeValue(mat_type); %Set text
                material_node.setAttributeNode(attr); %Add attribute
                
                for j=1:1:numel(FEB_struct.Materials{i}.Mats)
                    
                    generation_node = docNode.createElement('generation'); %create entry
                    generation_node = material_node.appendChild(generation_node); %add entry
                    
                    %Generation ID
                    attr = docNode.createAttribute('id'); %Create attribute
                    attr.setNodeValue(num2str(j)); %Set text
                    generation_node.setAttributeNode(attr); %Add attribute
                    
                    %Define start time
                    start_time_node = docNode.createElement('start_time'); %create entry
                    start_time_node = generation_node.appendChild(start_time_node); %add entry
                    start_time_val=Mat_i.start_times(j);
                    t_form=repmat('%6.7e, ',1,size(start_time_val,2)); t_form=t_form(1:end-2);
                    start_time_node.appendChild(docNode.createTextNode(sprintf(t_form,start_time_val))); %append data text child
                    
                    %Define solid
                    solid_node = docNode.createElement('solid'); %create entry
                    solid_node = generation_node.appendChild(solid_node); %add entry
                    
                    %Mat Type
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue(Mat_i.Mats{j}.type); %Set text
                    solid_node.setAttributeNode(attr); %Add attribute
                    
                    %%
                    %                     mat_props=Mat_i.Mats{j}.props;
                    %                     mat_prop_vals=Mat_i.Mats{j}.vals;
                    %
                    %                     %Material properties
                    %                     for q=1:1:numel(mat_props)
                    %                         mat_prop=mat_props{q};
                    %                         mat_prop_val=mat_prop_vals{q};
                    %                         prop_node = docNode.createElement(mat_prop); %create entry
                    %                         prop_node = solid_node.appendChild(prop_node); %add entry
                    %                         t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                    %                         prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                    %                     end
                    
                    %Mat properties and attributes
                    if isfield(Mat_i.Mats{j},'props')
                        mat_props=Mat_i.Mats{j}.props;
                        mat_prop_vals=Mat_i.Mats{j}.vals;
                        
                        %Access property attributes
                        if isfield(Mat_i.Mats{j},'PropAttrName')
                            mat_prop_attr_name=Mat_i.Mats{j}.PropAttrName;
                            mat_prop_attr_val=Mat_i.Mats{j}.PropAttrVal;
                        else
                            mat_prop_attr_name=[];
                            mat_prop_attr_val=[];
                        end
                        
                        %Material properties
                        for q=1:1:numel(mat_props)
                            mat_prop=mat_props{q};
                            mat_prop_val=mat_prop_vals{q};
                            prop_node = docNode.createElement(mat_prop); %create entry
                            prop_node = solid_node.appendChild(prop_node); %add entry
                            t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                            prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                            
                            %Add potential property attributes (e.g. property load curves)
                            if ~isempty(mat_prop_attr_name)
                                if q<=numel(mat_prop_attr_name)
                                    currentPropAttrName=mat_prop_attr_name{q};
                                    if ~isempty(currentPropAttrName)
                                        currentPropAttrVal=mat_prop_attr_val{q};
                                        
                                        t_form=repmat('%6.7e, ',1,size(currentPropAttrVal,2));
                                        t_form=t_form(1:end-2);
                                        
                                        attr = docNode.createAttribute(currentPropAttrName); %Create attribute
                                        attr.setNodeValue(sprintf(t_form,currentPropAttrVal)); %Set text
                                        prop_node.setAttributeNode(attr); %Add attribute
                                    end
                                end
                            end
                        end
                    end
                    
                    %%
                    
                    %Fiber spec
                    mat_aniso_type=Mat_i.Mats{j}.aniso_type;
                    switch mat_aniso_type
                        case 'none'
                            
                        case 'fiber'
                            prop_node = docNode.createElement('fiber'); %create entry
                            prop_node = solid_node.appendChild(prop_node); %add entry
                            attr = docNode.createAttribute('type'); %Create attribute
                            attr.setNodeValue('user'); %Set text
                            prop_node.setAttributeNode(attr); %Add attribute
                        case 'mat_axis'
                            prop_node = docNode.createElement('mat_axis'); %create entry
                            prop_node = solid_node.appendChild(prop_node); %add entry
                            attr = docNode.createAttribute('type'); %Create attribute
                            attr.setNodeValue('vector'); %Set text
                            prop_node.setAttributeNode(attr); %Add attribute
                    end
                    
                end
            end
            
            %%
        otherwise
            
            %Mat ID
            attr = docNode.createAttribute('id'); %Create attribute
            attr.setNodeValue(num2str(mat_ind(i))); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Name
            attr = docNode.createAttribute('name'); %Create attribute
            attr.setNodeValue(['mat_',num2str(i)]); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Type
            attr = docNode.createAttribute('type'); %Create attribute
            attr.setNodeValue(mat_type); %Set text
            material_node.setAttributeNode(attr); %Add attribute
           
            mat_aniso_type=Mat_i.aniso_type;
  
            %%
            
            %Mat properties and attributes
            if isfield(Mat_i,'props')
                mat_props=Mat_i.props;
                mat_prop_vals=Mat_i.vals;
                
                %Access property attributes
                if isfield(Mat_i,'PropAttrName')
                    mat_prop_attr_name=Mat_i.PropAttrName;
                    mat_prop_attr_val=Mat_i.PropAttrVal;
                else
                    mat_prop_attr_name=[];
                    mat_prop_attr_val=[];
                end
                
                %Material properties
                for q=1:1:numel(mat_props)
                    mat_prop=mat_props{q};
                    mat_prop_val=mat_prop_vals{q};
                    prop_node = docNode.createElement(mat_prop); %create entry
                    prop_node = material_node.appendChild(prop_node); %add entry
                    t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                    
                    %Add potential property attributes (e.g. property load curves)
                    if ~isempty(mat_prop_attr_name)
                        if q<=numel(mat_prop_attr_name)
                            currentPropAttrName=mat_prop_attr_name{q};
                            if ~isempty(currentPropAttrName)
                                currentPropAttrVal=mat_prop_attr_val{q};
                                
                                t_form=repmat('%6.7e, ',1,size(currentPropAttrVal,2));
                                t_form=t_form(1:end-2);
                                
                                attr = docNode.createAttribute(currentPropAttrName); %Create attribute
                                attr.setNodeValue(sprintf(t_form,currentPropAttrVal)); %Set text
                                prop_node.setAttributeNode(attr); %Add attribute
                            end
                        end
                    end
                    
                end
            end
            
            %%
            %Fiber spec
            switch mat_aniso_type
                case 'none' %No fibre type specified, do nothing
                    
                case 'fiber'
                    prop_node = docNode.createElement('fiber'); %create entry
                    prop_node = material_node.appendChild(prop_node); %add entry
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue('user'); %Set text
                    prop_node.setAttributeNode(attr); %Add attribute
                case 'mat_axis'
                    prop_node = docNode.createElement('mat_axis'); %create entry
                    prop_node = material_node.appendChild(prop_node); %add entry
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue('vector'); %Set text
                    prop_node.setAttributeNode(attr); %Add attribute
            end
    end
    %     prop_node.appendChild(docNode.createTextNode(' ')); %append data text child
    if FEB_struct.disp_opt==1;
        waitbar(i/numel(mat_ind),hw,'Adding material level...');
    end
end
if FEB_struct.disp_opt==1;
    close(hw);
    drawnow;
end