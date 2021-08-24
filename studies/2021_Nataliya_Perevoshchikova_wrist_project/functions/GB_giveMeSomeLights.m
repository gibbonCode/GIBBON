function [] = GB_giveMeSomeLights(thanks,parentName)
% Bardiya Akhbari
% This is an add-on for plotboneVRML function

if thanks == 1
    if ~exist('parentName','var')
        colormap(jet);
        light('Position',[150 0 0],'Style','local');
        light('Position',[0 150 0],'Style','local');
        light('Position',[0 0 150],'Style','local');
        light('Position',[200 200 100],'Style','local');
    else
        light('Position',[150 0 0],'Style','local','Parent',parentName);
        light('Position',[0 150 0],'Style','local','Parent',parentName);
        light('Position',[0 0 150],'Style','local','Parent',parentName);
        light('Position',[200 200 100],'Style','local','Parent',parentName);
    end
elseif thanks == 66
    light('Position',[0 0 -50],'Style','local');
    light('Position',[0 -150 50],'Style','local');
    light('Position',[0 150 50],'Style','local');
%     light('Position',[200 50 100],'Style','local');    
elseif thanks == 77
    light('Position',[150 0 0],'Style','local');
    light('Position',[0 -150 0],'Style','local');
    light('Position',[0 0 150],'Style','local');
    light('Position',[200 50 100],'Style','local');
elseif thanks == 88
    light('Position',[-150 0 0],'Style','local');
    light('Position',[0 50 0],'Style','local');
    light('Position',[0 0 150],'Style','local');
    light('Position',[-200 50 -100],'Style','local');   
elseif thanks == 99
    light('Position',[-150 0 0],'Style','local');
    light('Position',[0 150 0],'Style','local');
    light('Position',[0 0 150],'Style','local');
    light('Position',[200 200 100],'Style','local');
end