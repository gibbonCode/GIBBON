function [docNode]=addBoundaryLevel_FEB(docNode,FEB_struct)

% if  isfield(FEB_struct,'Boundary')
%     
%     disp('Adding Boundary level')
%     febio_spec = docNode.getDocumentElement;
%     
%     %Check for/create Boundary level
%     if docNode.getElementsByTagName('Boundary').getLength==0; %level does not exist yet
%         ElementAddNode = docNode.createElement('Boundary');
%         febio_spec.appendChild(ElementAddNode);
%         BoundaryNode = docNode.getElementsByTagName('Boundary').item(0);
%         IsBoundaryNodeNew=1;
%     else
%         BoundaryNode = docNode.getElementsByTagName('Boundary').item(0);
%         IsBoundaryNodeNew=0;
%     end
%     
%     %% Adding boundary components
%     
%     [docNode]=addBoundaryComponents_FEB(docNode,BoundaryNode,FEB_struct);       
%     
% end

disp('Adding Boundary level')

rootNode = docNode.getDocumentElement;

boundaryNode = docNode.createElement('Boundary');
boundaryNode = rootNode.appendChild(boundaryNode);

%% Adding boundary components

[docNode]=addBoundaryComponents_FEB(docNode,boundaryNode,FEB_struct);

