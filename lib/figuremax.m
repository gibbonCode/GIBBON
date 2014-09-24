function [fig]=figuremax(varargin)

% function [fig]=figuremax(Cbg,Cdef)
% ------------------------------------------------------------------------
%
% This function opens a maximized figure window using the color settings
% specified.
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% Last update 14/03/2012
%------------------------------------------------------------------------
%%

switch nargin
    case 0
        Cbg='w';
        Cdef='white';
        scr_offset=60; %e.g. quick and dirty way to cope with bottom taskbar in windows
    case 1
        Cbg=varargin{1};
        Cdef='white';
        scr_offset=60; %e.g. quick and dirty way to cope with bottom taskbar in windows
    case 2
        Cbg=varargin{1};
        Cdef=varargin{2};
        scr_offset=60; %e.g. quick and dirty way to cope with bottom taskbar in windows        
    case 3
        Cbg=varargin{1};
        Cdef=varargin{2};
        scr_offset=varargin{3}; %e.g. quick and dirty way to cope with bottom taskbar in windows        
end
%Open figure with handle
fig=figure; 
set(fig,'renderer','OpenGL'); %Default renderer changed, options: painters | zbuffer | OpenGL

%Setting renderer. For RGB and colormap driven There are some bugs for the hardware option, hence changed here
opengl software; 

%Specify color scheme
colordef (fig,Cdef); 

%Set renderer, set background color, maximize figure size
% set(fig,'Color',Cbg,'units','normalized','outerposition',[0 0 1 1]);

scrsz = get(0,'ScreenSize'); 
set(fig,'Color',Cbg,'units','pixels','outerposition',[1+scr_offset scr_offset scrsz(3)-2*scr_offset scrsz(4)-2*scr_offset]);
hold off;
 
end
