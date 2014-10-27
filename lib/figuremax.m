function [hf]=figuremax(varargin)

% function [fig]=figuremax(Cbg,Cdef)
% ------------------------------------------------------------------------
%
% This function opens a maximized figure window using the color settings
% specified.
% 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 30/09/2014
%------------------------------------------------------------------------

%%

scr_offset_default= round(max(get(0,'ScreenSize'))/15);

switch nargin
    case 0
        Cbg='w';
        Cdef='white';
        scr_offset=scr_offset_default; %e.g. quick and dirty way to cope with bottom taskbar in windows
    case 1
        Cbg=varargin{1};
        Cdef='white';
        scr_offset=scr_offset_default; %e.g. quick and dirty way to cope with bottom taskbar in windows
    case 2
        Cbg=varargin{1};
        Cdef=varargin{2};
        scr_offset=scr_offset_default; %e.g. quick and dirty way to cope with bottom taskbar in windows        
    case 3
        Cbg=varargin{1};
        Cdef=varargin{2};
        scr_offset=varargin{3}; %e.g. quick and dirty way to cope with bottom taskbar in windows        
end

%Open figure with handle
hf=figure; 
set(hf,'name','GIBBON');
set(hf,'renderer','OpenGL'); %Default renderer changed, options: painters | zbuffer | OpenGL

%Specify color scheme
colordef(hf,Cdef); 
set(hf,'Color',Cbg);

%Setting renderer. For RGB and colormap driven There are some bugs for the hardware option, hence changed here
if ispc
%     opengl software;
    %On UNIX systems, start MATLAB with the command, matlab softwareopengl
end

%Set renderer, set background color, maximize figure size
if scr_offset==0
    try 
        maxFig(hf);
    catch
        set(hf,'units','normalized','outerposition',[0 0 1 1]);
    end
else
    scrsz = get(0,'ScreenSize');
    set(hf,'units','pixels','outerposition',[scr_offset scr_offset scrsz(3)-2*scr_offset scrsz(4)-2*scr_offset]);
end

hold off;
 
end
