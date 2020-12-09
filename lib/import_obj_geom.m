function [F,V]=import_obj_geom(fileName)

T=txtfile2cell(fileName);
V=cell2mat(cellfun(@(x) (sscanf(x,'v %f %f %f')'),T,'UniformOutput',0));
F=cell2mat(cellfun(@(x) (sscanf(x,'f %d %d %d')'),T,'UniformOutput',0));