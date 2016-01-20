function maximizeFigureWindow(varargin)

%% Parse input

switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
end

if isempty(hf)
    hf=gcf; %Take the current figure window handle if the handle is empty
end

%% Turn off warning (JavaFrame will become obsolete, this will lead to error in the future)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

%% Maximize using JaveFrame
try
    pause(0.01);     
    jFrame=get(hf,'JavaFrame');
%     set(jFrame,'Maximized',1);
    jFrame.setMaximized(1);
catch exception    
    rethrow(exception)
end

