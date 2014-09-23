function export_STL_txt(fileName,stlStruct)

% ------------------------------------------------------------------------
% function export_STL_txt(fileName,stlStruct)
%
% This function generates an STL file using the input structure stlStruct.
% Multiple solids are supported as cell entries within the structure.
%
% Example if the structure array for a single solid:
%
% stlStruct =
%    solidNames: {'standford_bunny'} %Cell containing solid names
% solidVertices: {[8745x3 double]} %Cell constaining vertices for each solid
%    solidFaces: {[2915x3 double]} %Cell containing faces for each solid
%  solidNormals: {[2915x3 double]} %Cell containing normals for each solid
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/04/08
%
% CHANGE LOG:
% 2014/04/08
% * Created
%
% ------------------------------------------------------------------------

%% Write STL file in txt format

%Remove file if it exists
if exist(fileName,'file')==2; 
    delete(fileName);
end

%Append STLs to file
writeOpt='a';
for q=1:1:numel(stlStruct.solidNames);    
    patch2STL(fileName,stlStruct.solidVertices{q},stlStruct.solidFaces{q},stlStruct.solidNormals{q},stlStruct.solidNames{q},writeOpt);
end

